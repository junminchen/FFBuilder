import os
import numpy as np
from FFutils import read_pdb_skeleton
from pyscf import gto, dft, scf
try:
    import gpu4pyscf
    from gpu4pyscf.dft import rks
    HAS_GPU = True
except ImportError:
    from pyscf import dft
    HAS_GPU = False
from scipy.optimize import minimize

def run_chelpg_bf4():
    name = "BF4"
    total_charge = -1.0
    
    # 1. Setup Molecule (from your PDB)
    symbols, coords = read_pdb_skeleton(f"molecules/{name}/{name}.pdb")
    # Convert to Bohr for PySCF
    ang2bohr = 1.8897259886
    coords_bohr = coords * ang2bohr
    
    mol = gto.Mole()
    mol.atom = "\n".join([f"{s} {c[0]} {c[1]} {c[2]}" for s, c in zip(symbols, coords_bohr)])
    mol.unit = 'Bohr'
    mol.basis = '6-31G*' # Standard OPLS/RESP basis
    mol.charge = total_charge
    mol.build()
    
    print(f"--- CHelpG Calculation for {name} ---")
    print(f"Basis: {mol.basis}, Charge: {mol.charge}")
    
    # 2. DFT Calculation
    if HAS_GPU:
        print("Using GPU4PySCF...")
        mf = rks.RKS(mol, xc='b3lyp').to_gpu()
    else:
        print("Using CPU PySCF...")
        mf = dft.RKS(mol, xc='b3lyp')
    
    mf.kernel()
    dm = mf.make_rdm1()
    if HAS_GPU: dm = dm.get()
    
    # 3. Generate CHelpG Grid
    # Cubic grid with 0.3 Angstrom spacing
    spacing = 0.3 * ang2bohr
    box_size = 4.0 * ang2bohr # 4A from center
    
    x = np.arange(-box_size, box_size, spacing)
    y = np.arange(-box_size, box_size, spacing)
    z = np.arange(-box_size, box_size, spacing)
    gx, gy, gz = np.meshgrid(x, y, z)
    grid_pts = np.c_[gx.ravel(), gy.ravel(), gz.ravel()]
    
    # Filter grid points: 
    # Must be outside VdW radii and within 2.8A of any atom (CHelpG standard)
    # VdW radii (Bohr): B=3.59, F=2.78 (approximate)
    vdw = {'B': 1.9 * ang2bohr, 'F': 1.47 * ang2bohr}
    
    valid_mask = np.zeros(len(grid_pts), dtype=bool)
    for i in range(mol.natm):
        r = np.linalg.norm(grid_pts - mol.atom_coords()[i], axis=1)
        # Outside VdW
        if i == 0: # B
            valid_mask |= (r > vdw['B']) & (r < 2.8 * ang2bohr + vdw['B'])
        else: # F
            valid_mask |= (r > vdw['F']) & (r < 2.8 * ang2bohr + vdw['F'])

    # Final cleanup: ensure no point is inside ANY VdW sphere
    for i in range(mol.natm):
        r = np.linalg.norm(grid_pts - mol.atom_coords()[i], axis=1)
        valid_mask &= (r > vdw.get(mol.atom_symbol(i), 1.5*ang2bohr))

    grid_pts = grid_pts[valid_mask]
    print(f"Number of CHelpG grid points: {len(grid_pts)}")
    
    # 4. Calculate Target ESP
    # V_total = V_nuc + V_elec
    v_nuc = np.zeros(len(grid_pts))
    for i in range(mol.natm):
        v_nuc += mol.atom_charge(i) / np.linalg.norm(grid_pts - mol.atom_coords()[i], axis=1)
    
    v_elec = np.zeros(len(grid_pts))
    # Batch ESP calculation for efficiency
    for i, pt in enumerate(grid_pts):
        mol.set_rinv_origin(pt)
        v_elec[i] = -np.einsum('ij,ij', dm, mol.intor('int1e_rinv'))
    target_esp = v_nuc + v_elec

    # 5. Symmetry-Constrained Fitting
    # x[0] = B charge, x[1] = F charge (shared by all 4)
    def objective(x):
        q = np.array([x[0], x[1], x[1], x[1], x[1]])
        v_calc = np.zeros(len(grid_pts))
        for i in range(mol.natm):
            v_calc += q[i] / np.linalg.norm(grid_pts - mol.atom_coords()[i], axis=1)
        
        chi_sq = np.sum((v_calc - target_esp)**2)
        return chi_sq

    cons = ({'type': 'eq', 'fun': lambda x: x[0] + 4*x[1] - total_charge})
    
    # Initial guess: B=0.6, F=-0.4 (Mulliken-like)
    res = minimize(objective, [0.6, -0.4], constraints=cons, method='SLSQP')
    
    print("\n--- Final CHelpG Results for BF4 ---")
    print(f"B Charge: {res.x[0]:.4f}")
    print(f"F Charge: {res.x[1]:.4f}")
    print(f"Total: {res.x[0] + 4*res.x[1]:.4f}")
    print(f"Fitting RMS Error: {np.sqrt(res.fun/len(grid_pts)):.6f} Hartree/e")

if __name__ == "__main__":
    run_chelpg_bf4()
