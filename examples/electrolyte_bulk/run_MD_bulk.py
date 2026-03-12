import os
import sys
import openmm
import openmm.app as app
import openmm.unit as unit

# ================= MD 仿真配置 =================
PDB_FILE = 'init_system.pdb'
FF_FILES = ['forcefields/opls_solvent.xml', 'forcefields/opls_salt.xml']
TEMPERATURE = 298.15 * unit.kelvin
PRESSURE = 1.0 * unit.bar
TIMESTEP = 0.002 * unit.picoseconds
TOTAL_STEPS = 5000000  # 10 ns
REPORT_STEPS = 5000    # 10 ps
# ===============================================

def run_npt_simulation():
    print('--- Starting Bulk MD (NPT) Simulation ---')
    if not os.path.exists(PDB_FILE):
        print(f'Error: {PDB_FILE} not found.')
        return

    pdb = app.PDBFile(PDB_FILE)
    topology = pdb.topology
    topology.setUnitCellDimensions((3.0, 3.0, 3.0)*unit.nanometer)

    for res in topology.residues():
        at = {a.name.strip(): a for a in res.atoms()}
        if res.name == 'PF6':
            p = at['P01']
            for i in range(2, 8): topology.addBond(p, at[f'F{i:02d}'])
        elif res.name == 'PC':
            c00, c01, c02, o03, c04, o05, o06 = at['C00'], at['C01'], at['C02'], at['O03'], at['C04'], at['O05'], at['O06']
            topology.addBond(c00, c01); topology.addBond(c01, c02); topology.addBond(c02, o03)
            topology.addBond(o03, c04); topology.addBond(c04, o05); topology.addBond(c04, o06)
            topology.addBond(c01, o06)
            for h in ['H07', 'H08', 'H09']: topology.addBond(c00, at[h])
            topology.addBond(c01, at['H10'])
            for h in ['H11', 'H12']: topology.addBond(c02, at[h])
        elif res.name == 'DMC':
            c00, o01, c02, o03, o04, c05 = at['C00'], at['O01'], at['C02'], at['O03'], at['O04'], at['C05']
            topology.addBond(c00, o01); topology.addBond(o01, c02); topology.addBond(c02, o03)
            topology.addBond(c02, o04); topology.addBond(o04, c05)
            for h in ['H06', 'H07', 'H08']: topology.addBond(c00, at[h])
            for h in ['H09', 'H10', 'H11']: topology.addBond(c05, at[h])

    forcefield = app.ForceField(*FF_FILES)
    system = forcefield.createSystem(topology, nonbondedMethod=app.PME, 
                                    nonbondedCutoff=1.0*unit.nanometer, constraints=app.HBonds)
    system.addForce(openmm.MonteCarloBarostat(PRESSURE, TEMPERATURE))
    
    integrator = openmm.LangevinMiddleIntegrator(TEMPERATURE, 1/unit.picosecond, TIMESTEP)
    
    try:
        platform = openmm.Platform.getPlatformByName('CUDA')
        properties = {'Precision': 'mixed'}
        simulation = app.Simulation(topology, system, integrator, platform, properties)
        simulation.context.setPositions(pdb.positions)
        print('Platform: CUDA')
    except:
        platform = openmm.Platform.getPlatformByName('Reference')
        simulation = app.Simulation(topology, system, integrator, platform)
        simulation.context.setPositions(pdb.positions)
        print('Platform: Reference')

    print('Minimizing...')
    simulation.minimizeEnergy()
    
    print(f'Running {TOTAL_STEPS} steps...')
    simulation.reporters.append(app.StateDataReporter(sys.stdout, 100, step=True, potentialEnergy=True, density=True))
    simulation.step(TOTAL_STEPS)
    print('Success!')

if __name__ == "__main__":
    run_npt_simulation()
