import os
import numpy as np
from FFutils import read_pdb_skeleton
from pyscf import gto, dft, scf
import gpu4pyscf
from gpu4pyscf.dft import rks
from scipy.optimize import minimize

def debug_bf4_v3():
    name = "BF4"
    total_charge = -1
    
    symbols, coords = read_pdb_skeleton(f"molecules/{name}/{name}.pdb")
    coords_bohr = coords * 1.8897259886
    
    mol = gto.Mole()
    mol.atom = "\n".join([f"{s} {c[0]} {c[1]} {c[2]}" for s, c in zip(symbols, coords_bohr)])
    mol.unit = 'Bohr'; mol.basis = '6-31G*'; mol.charge = total_charge
    mol.build()
    
    print("Running DFT...")
    mf = rks.RKS(mol, xc='b3lyp').to_gpu()
    mf.kernel()
    dm_cpu = mf.make_rdm1().get()
    
    # Mulliken charges for reference
    mulliken = mf.to_cpu().mulliken_meta()[1]
    print(f"Mulliken charges: B={mulliken[0]:.4f}, F={np.mean(mulliken[1:]):.4f}")
    
    # Grid
    grid_pts = []
    np.random.seed(42)
    vdw = {'B': 1.9*1.89, 'F': 1.47*1.89}
    for i in range(mol.natm):
        for f in [1.2, 1.4, 1.6, 1.8]: # Use closer shells to reduce shielding
            r = vdw.get(mol.atom_symbol(i), 1.5*1.89) * f
            pts = np.random.normal(0, 1, (400, 3))
            pts /= np.linalg.norm(pts, axis=1)[:, None]
            grid_pts.append(mol.atom_coords()[i] + pts * r)
    grid_pts = np.vstack(grid_pts)
    
    # Target ESP
    v_nuc = np.zeros(len(grid_pts))
    for i in range(mol.natm):
        v_nuc += mol.atom_charge(i) / np.linalg.norm(grid_pts - mol.atom_coords()[i], axis=1)
    
    v_elec = np.zeros(len(grid_pts))
    for i, pt in enumerate(grid_pts):
        mol.set_rinv_origin(pt)
        v_elec[i] = -np.einsum('ij,ij', dm_cpu, mol.intor('int1e_rinv'))
    target_esp = v_nuc + v_elec

    # Fitting with Symmetry and Normalized Chi-sq
    def objective(x):
        q = np.array([x[0], x[1], x[1], x[1], x[1]])
        v_calc = np.zeros(len(grid_pts))
        for i in range(mol.natm):
            v_calc += q[i] / np.linalg.norm(grid_pts - mol.atom_coords()[i], axis=1)
        
        # NORMALIZE chi_sq by number of points
        chi_sq = np.sum((v_calc - target_esp)**2) / len(grid_pts)
        
        # Hyperbolic restraint (Very strong for buried B)
        a_buried = 0.05 
        a_f = 0.001
        b = 0.1
        restraint = a_buried * (np.sqrt(x[0]**2 + b**2) - b) + 4 * a_f * (np.sqrt(x[1]**2 + b**2) - b)
        return chi_sq + restraint

    cons = ({'type': 'eq', 'fun': lambda x: x[0] + 4*x[1] - total_charge})
    # Start from Mulliken
    res = minimize(objective, [mulliken[0], np.mean(mulliken[1:])], constraints=cons, method='SLSQP')
    
    print("\n--- Final V3 RESP Charges ---")
    print(f"B Charge: {res.x[0]:.4f}")
    print(f"F Charge: {res.x[1]:.4f}")
    print(f"Total: {res.x[0] + 4*res.x[1]:.4f}")

if __name__ == "__main__":
    debug_bf4_v3()
