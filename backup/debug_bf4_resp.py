import os
import numpy as np
from FFutils import read_pdb_skeleton
from pyscf import gto, dft, scf
import gpu4pyscf
from gpu4pyscf.dft import rks
from scipy.optimize import minimize

def debug_bf4():
    name = "BF4"
    smiles = "[B-:1]([F:2])([F:3])([F:4])([F:5])"
    total_charge = -1
    
    symbols, coords = read_pdb_skeleton(f"molecules/{name}/{name}.pdb")
    print(f"Symbols: {symbols}")
    
    # Unit: Bohr
    coords_bohr = coords * 1.8897259886
    mol_str = "\n".join([f"{s} {c[0]} {c[1]} {c[2]}" for s, c in zip(symbols, coords_bohr)])
    
    mol = gto.Mole()
    mol.atom = mol_str
    mol.unit = 'Bohr'
    mol.basis = '6-31G*'
    mol.charge = total_charge
    mol.build()
    
    print("Running DFT...")
    mf = rks.RKS(mol, xc='b3lyp').to_gpu()
    mf.kernel()
    dm_cpu = mf.make_rdm1().get()
    
    # Grid
    vdw_radii_bohr = {'B': 1.9*1.89, 'F': 1.47*1.89}
    grid_pts = []
    np.random.seed(42)
    for i in range(mol.natm):
        for factor in [1.4, 1.6, 1.8, 2.0]:
            r = vdw_radii_bohr.get(mol.atom_symbol(i), 1.5 * 1.89) * factor
            pts = np.random.normal(0, 1, (200, 3))
            pts /= np.linalg.norm(pts, axis=1)[:, None]
            grid_pts.append(mol.atom_coords()[i] + pts * r)
    grid_pts = np.vstack(grid_pts)
    
    print(f"Number of Grid points: {len(grid_pts)}")
    
    # ESP calculation
    v_nuc = np.zeros(len(grid_pts))
    for i in range(mol.natm):
        dist = np.linalg.norm(grid_pts - mol.atom_coords()[i], axis=1)
        v_nuc += mol.atom_charge(i) / dist
    
    v_elec = np.zeros(len(grid_pts))
    for i, pt in enumerate(grid_pts):
        mol.set_rinv_origin(pt)
        rinv_mat = mol.intor('int1e_rinv')
        v_elec[i] = -np.einsum('ij,ij', dm_cpu, rinv_mat)
        
    target_esp = v_nuc + v_elec
    
    print(f"Avg V_nuc: {np.mean(v_nuc):.4f}")
    print(f"Avg V_elec: {np.mean(v_elec):.4f}")
    print(f"Avg Target ESP: {np.mean(target_esp):.4f}")
    
    # Fitting
    a = 0.0005
    b = 0.1
    def objective(q):
        v_calc = np.zeros(len(grid_pts))
        for i in range(mol.natm):
            dist = np.linalg.norm(grid_pts - mol.atom_coords()[i], axis=1)
            v_calc += q[i] / dist
        chi_sq = np.sum((v_calc - target_esp)**2)
        restraint = 0
        for i in range(mol.natm):
            if mol.atom_symbol(i) != 'H':
                restraint += a * (np.sqrt(q[i]**2 + b**2) - b)
        return chi_sq + restraint
    
    cons = ({'type': 'eq', 'fun': lambda q: np.sum(q) - total_charge})
    res = minimize(objective, np.zeros(mol.natm), constraints=cons, method='SLSQP')
    
    q_final = res.x
    print("\n--- Final RESP Charges ---")
    for i, s in enumerate(symbols):
        print(f"Atom {i+1} {s}: {q_final[i]:.4f}")
    print(f"Total Charge: {np.sum(q_final):.4f}")

if __name__ == "__main__":
    debug_bf4()
