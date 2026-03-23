import os
import numpy as np
from FFutils import read_pdb_skeleton
from pyscf import gto, dft, lib
import gpu4pyscf
from gpu4pyscf.dft import rks
import resp

def verify_bf4():
    name = "BF4"
    total_charge = -1
    
    symbols, coords = read_pdb_skeleton(f"molecules/{name}/{name}.pdb")
    # Unit: Bohr
    coords_bohr = coords * 1.8897259886
    mol_str = "\n".join([f"{s} {c[0]} {c[1]} {c[2]}" for s, c in zip(symbols, coords_bohr)])
    
    mol = gto.Mole()
    mol.atom = mol_str
    mol.unit = 'Bohr'; mol.basis = '6-31G*'; mol.charge = total_charge
    mol.build()
    
    print("Running DFT for BF4...")
    mf = rks.RKS(mol, xc='b3lyp').to_gpu()
    mf.kernel()
    dm = mf.make_rdm1().get()
    
    # Grid Generation (using the same logic as our workflow)
    vdw = {'B': 1.9*1.89, 'F': 1.47*1.89}
    grid_pts = []
    np.random.seed(42)
    for i in range(mol.natm):
        for factor in [1.4, 1.6, 1.8, 2.0]:
            r = vdw.get(mol.atom_symbol(i), 1.5*1.89) * factor
            pts = np.random.normal(0, 1, (200, 3))
            pts /= np.linalg.norm(pts, axis=1)[:, None]
            grid_pts.append(mol.atom_coords()[i] + pts * r)
    grid_pts = np.vstack(grid_pts)
    
    # ESP calculation
    v_nuc = np.zeros(len(grid_pts))
    for i in range(mol.natm):
        v_nuc += mol.atom_charge(i) / np.linalg.norm(grid_pts - mol.atom_coords()[i], axis=1)
    
    v_elec = np.zeros(len(grid_pts))
    for i, pt in enumerate(grid_pts):
        mol.set_rinv_origin(pt)
        v_elec[i] = -np.einsum('ij,ij', dm, mol.intor('int1e_rinv'))
    target_esp = v_nuc + v_elec

    # Call the installed 'resp' library
    options = resp.RESPOptions(restraint_height=0.0005, restraint_width=0.1)
    
    print("Fitting charges with 'resp' library...")
    # The 'resp' package expects: coords (Natm, 3), target_esp (Ngrid), grid_pts (Ngrid, 3), total_charge
    # All in atomic units
    charges = resp.fit_charges(mol.atom_coords(), target_esp, grid_pts, total_charge, options)
    
    print("\n--- BF4 RESP Results (using lib) ---")
    for i, s in enumerate(symbols):
        print(f"Atom {i+1} {s}: {charges[i]:.4f}")
    print(f"Total: {np.sum(charges):.4f}")

if __name__ == "__main__":
    verify_bf4()
