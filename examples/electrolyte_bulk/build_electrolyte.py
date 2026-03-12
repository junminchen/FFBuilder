import os
import sys
import numpy as np
import openmm
import openmm.app as app
import openmm.unit as unit

# ================= 配置部分 =================
SOLVENTS = {"PC": 20, "DMC": 20} 
FF_FILES = ["forcefields/opls_solvent.xml", "forcefields/opls_salt.xml"]
# ===========================================

def build_simple_pdb():
    print("--- Building PDB (PC/DMC/LiPF6) ---")
    lines = []
    atom_idx = 1
    res_idx = 1
    box_size = 30.0

    def add_mol(name, coords, names, res_name):
        nonlocal atom_idx, res_idx
        center = np.random.uniform(2, box_size-2, 3)
        for n, c in zip(names, coords):
            lines.append(f"ATOM  {atom_idx:5d} {n:4s} {res_name:3s}  {res_idx:4d}    {c[0]+center[0]:8.3f}{c[1]+center[1]:8.3f}{c[2]+center[2]:8.3f}  1.00  0.00")
            atom_idx += 1
        res_idx += 1

    pc_names = ["C00", "C01", "C02", "O03", "C04", "O05", "O06", "H07", "H08", "H09", "H10", "H11", "H12"]
    pc_coords = np.random.normal(0, 1, (13, 3))
    for _ in range(SOLVENTS["PC"]): add_mol("PC", pc_coords, pc_names, "PC")

    dmc_names = ["C00", "O01", "C02", "O03", "O04", "C05", "H06", "H07", "H08", "H09", "H10", "H11"]
    dmc_coords = np.random.normal(0, 1, (12, 3))
    for _ in range(SOLVENTS["DMC"]): add_mol("DMC", dmc_coords, dmc_names, "DMC")

    pf6_names = ["P01", "F02", "F03", "F04", "F05", "F06", "F07"]
    pf6_coords = np.array([[0,0,0], [0,-1.6,0], [1.6,0,0], [-1.6,0,0], [0,0,1.6], [0,0,-1.6], [0,1.6,0]])
    for _ in range(5):
        add_mol("LiA", [[0,0,0]], ["Li01"], "LiA")
        add_mol("PF6", pf6_coords, pf6_names, "PF6")

    with open("init_system.pdb", "w") as f:
        f.write("\n".join(lines) + "\nEND")
    return "init_system.pdb"

def run_simulation(pdb_file):
    print("--- Running OpenMM Simulation ---")
    pdb = app.PDBFile(pdb_file)
    topology = pdb.topology
    
    for res in topology.residues():
        at = {a.name.strip(): a for a in res.atoms()}
        if res.name == "PF6":
            p = at["P01"]
            for i in range(2, 8): topology.addBond(p, at[f"F{i:02d}"])
        elif res.name == "PC":
            c00, c01, c02, o03, c04, o05, o06 = at["C00"], at["C01"], at["C02"], at["O03"], at["C04"], at["O05"], at["O06"]
            topology.addBond(c00, c01); topology.addBond(c01, c02); topology.addBond(c02, o03)
            topology.addBond(o03, c04); topology.addBond(c04, o05); topology.addBond(c04, o06)
            topology.addBond(c01, o06)
            for h in ["H07", "H08", "H09"]: topology.addBond(c00, at[h])
            topology.addBond(c01, at["H10"])
            for h in ["H11", "H12"]: topology.addBond(c02, at[h])
        elif res.name == "DMC":
            c00, o01, c02, o03, o04, c05 = at["C00"], at["O01"], at["C02"], at["O03"], at["O04"], at["C05"]
            topology.addBond(c00, o01); topology.addBond(o01, c02); topology.addBond(c02, o03)
            topology.addBond(c02, o04); topology.addBond(o04, c05)
            for h in ["H06", "H07", "H08"]: topology.addBond(c00, at[h])
            for h in ["H09", "H10", "H11"]: topology.addBond(c05, at[h])

    forcefield = app.ForceField(*FF_FILES)
    system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff)
    integrator = openmm.VerletIntegrator(0.001*unit.picoseconds)

    # 真正的初始化测试
    platform = None
    sim = None
    try:
        platform = openmm.Platform.getPlatformByName("CUDA")
        sim = app.Simulation(topology, system, integrator, platform)
        # 尝试触发 Context 创建以检查 CUDA 驱动兼容性
        sim.context.setPositions(pdb.positions)
        print("Platform: CUDA (Enabled & Working)")
    except Exception as e:
        print(f"CUDA failed to initialize ({e}). Falling back to Reference.")
        platform = openmm.Platform.getPlatformByName("Reference")
        integrator = openmm.VerletIntegrator(0.001*unit.picoseconds) # 重置
        sim = app.Simulation(topology, system, integrator, platform)
        sim.context.setPositions(pdb.positions)

    sim.minimizeEnergy()
    print(f"Success! Final Potential Energy: {sim.context.getState(getEnergy=True).getPotentialEnergy()}")

if __name__ == "__main__":
    run_simulation(build_simple_pdb())
