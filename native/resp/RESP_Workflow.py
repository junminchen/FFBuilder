import os
import sys
import numpy as np

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from utils.FFutils import read_pdb_skeleton, generate_ff_xml

# Optional dependencies
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

try:
    import pyscf
    from pyscf import gto, dft, lib, scf
    import gpu4pyscf
    from gpu4pyscf.dft import rks
    PYSCF_AVAILABLE = True
except ImportError:
    PYSCF_AVAILABLE = False

# SCALE_FACTOR = 1.0 (Original RESP)
SCALE_FACTOR = 1.0

def get_symmetry_groups(symbols, coords):
    """Simple symmetry grouping based on symbol and distance to center for small ions like BF4/PF6."""
    natm = len(symbols)
    if natm == 1: return [[0]]
    
    # Calculate distance to center of mass (or geometric center)
    center = np.mean(coords, axis=0)
    dists = np.linalg.norm(coords - center, axis=1)
    
    groups = {}
    for i in range(natm):
        # Key is (Symbol, Rounded Distance)
        key = (symbols[i], round(dists[i], 2))
        if key not in groups: groups[key] = []
        groups[key].append(i)
    
    return list(groups.values())

def calculate_real_resp(symbols, coords, total_charge=0):
    """Robust CHelpG-style fitting with symmetry and PySCF ESP calculation."""
    if not PYSCF_AVAILABLE:
        print("PySCF not available, returning dummy charges.")
        return np.random.uniform(-0.1, 0.1, len(symbols))

    np.random.seed(42)
    ang2bohr = 1.8897259886
    coords_bohr = coords * ang2bohr
    
    mol = gto.Mole()
    mol.atom = "\n".join([f"{s} {c[0]} {c[1]} {c[2]}" for s, c in zip(symbols, coords_bohr)])
    mol.unit = 'Bohr'; mol.basis = '6-31G*'; mol.charge = total_charge
    mol.build()
    
    try:
        print(f"Running DFT calculation for {len(symbols)} atoms...")
        mf = rks.RKS(mol, xc='b3lyp').to_gpu()
        mf.kernel()
        dm = mf.make_rdm1().get()
    except Exception as e:
        print(f"GPU failed or not found, using CPU: {e}")
        mf = dft.RKS(mol, xc='b3lyp')
        mf.kernel()
        dm = mf.make_rdm1()

    # CHelpG Grid Generation (Cubic grid)
    spacing = 0.3 * ang2bohr
    # Box around the molecule plus 3.0 Angstrom buffer
    min_c = np.min(mol.atom_coords(), axis=0) - 3.0 * ang2bohr
    max_c = np.max(mol.atom_coords(), axis=0) + 3.0 * ang2bohr
    
    x = np.arange(min_c[0], max_c[0], spacing)
    y = np.arange(min_c[1], max_c[1], spacing)
    z = np.arange(min_c[2], max_c[2], spacing)
    gx, gy, gz = np.meshgrid(x, y, z)
    grid_pts = np.c_[gx.ravel(), gy.ravel(), gz.ravel()]
    
    # Standard CHelpG VdW Radii (Bohr)
    vdw_map = {'H': 1.2*ang2bohr, 'C': 1.7*ang2bohr, 'N': 1.55*ang2bohr, 'O': 1.5*ang2bohr, 
               'F': 1.47*ang2bohr, 'P': 1.8*ang2bohr, 'S': 1.8*ang2bohr, 'B': 1.9*ang2bohr, 'Li': 1.82*ang2bohr}
    
    valid_mask = np.zeros(len(grid_pts), dtype=bool)
    # Point must be within 2.8A of at least one atom's VdW surface
    for i in range(mol.natm):
        r = np.linalg.norm(grid_pts - mol.atom_coords()[i], axis=1)
        r_vdw = vdw_map.get(mol.atom_symbol(i), 1.5*ang2bohr)
        valid_mask |= (r > r_vdw) & (r < r_vdw + 2.8 * ang2bohr)

    # Secondary check: point must not be inside ANY atom's VdW
    for i in range(mol.natm):
        r = np.linalg.norm(grid_pts - mol.atom_coords()[i], axis=1)
        r_vdw = vdw_map.get(mol.atom_symbol(i), 1.5*ang2bohr)
        valid_mask &= (r > r_vdw)

    grid_pts = grid_pts[valid_mask]
    print(f"Generated {len(grid_pts)} CHelpG grid points.")

    # Target ESP Calculation
    v_nuc = np.zeros(len(grid_pts))
    for i in range(mol.natm):
        v_nuc += mol.atom_charge(i) / np.linalg.norm(grid_pts - mol.atom_coords()[i], axis=1)
    
    v_elec = np.zeros(len(grid_pts))
    for i, pt in enumerate(grid_pts):
        mol.set_rinv_origin(pt)
        v_elec[i] = -np.einsum('ij,ij', dm, mol.intor('int1e_rinv'))
    target_esp = v_nuc + v_elec

    # Symmetry Grouping
    sym_groups = get_symmetry_groups(symbols, coords)
    num_groups = len(sym_groups)
    
    from scipy.optimize import minimize
    def objective(x):
        q = np.zeros(mol.natm)
        for g_idx, indices in enumerate(sym_groups):
            for i in indices: q[i] = x[g_idx]
        
        v_calc = np.zeros(len(grid_pts))
        for i in range(mol.natm):
            v_calc += q[i] / np.linalg.norm(grid_pts - mol.atom_coords()[i], axis=1)
        return np.sum((v_calc - target_esp)**2)

    cons = ({'type': 'eq', 'fun': lambda x: np.sum([x[g_idx] * len(sym_groups[g_idx]) for g_idx in range(num_groups)]) - total_charge})
    
    # Use zero as initial guess for stability
    res = minimize(objective, np.zeros(num_groups), constraints=cons, method='SLSQP', options={'ftol': 1e-12})
    
    final_q = np.zeros(mol.natm)
    for g_idx, indices in enumerate(sym_groups):
        for i in indices: final_q[i] = res.x[g_idx]
            
    print(f"CHelpG Charges calculated. Total: {np.sum(final_q):.4f}, RMS: {np.sqrt(res.fun/len(grid_pts)):.6e}")
    return final_q

def process_molecule(name, smiles, total_charge=0):
    print(f"Processing {name} ({smiles})...")
    pdb_path = f"molecules/{name}/{name}.pdb"
    if not os.path.exists(pdb_path):
        pdb_path = f"molecules/{name}.pdb"
    symbols, coords = read_pdb_skeleton(pdb_path)
    if symbols is None:
        if RDKIT_AVAILABLE:
            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            symbols = [a.GetSymbol() for a in mol.GetAtoms()]
            coords = mol.GetConformer().GetPositions()
        else:
            raise ValueError(f"No geometry for {name}")

    raw_charges = calculate_real_resp(symbols, coords, total_charge)
    return symbols, coords, raw_charges * SCALE_FACTOR, None
