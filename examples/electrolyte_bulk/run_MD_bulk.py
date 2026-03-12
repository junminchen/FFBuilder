import os
import sys
import openmm
import openmm.app as app
import openmm.unit as unit

# ================= MD 仿真配置 =================
# 1. 输入文件
PDB_FILE = "init_system.pdb"  # 由建模脚本生成的初始构型
FF_FILES = ["../../forcefields/opls_solvent.xml", "../../forcefields/opls_salt.xml"]

# 2. 模拟参数
TEMPERATURE = 298.15 * unit.kelvin
PRESSURE = 1.0 * unit.bar
TIMESTEP = 0.002 * unit.picoseconds
TOTAL_STEPS = 5000000  # 10 ns (5M * 2fs)
REPORT_STEPS = 5000    # 每 10 ps 记录一次数据

# 3. 输出文件
OUTPUT_TRAJ = "trajectory.dcd"
OUTPUT_LOG = "simulation.log"
CHECKPOINT = "checkpoint.chk"
# ===============================================

def run_npt_simulation():
    print("--- Starting Bulk MD (NPT) Simulation ---")
    
    # 1. 加载构型与力场
    if not os.path.exists(PDB_FILE):
        print(f"Error: {PDB_FILE} not found. Please build the system first.")
        return

    pdb = app.PDBFile(PDB_FILE)
    forcefield = app.ForceField(*FF_FILES)
    
    # 2. 创建 OpenMM 系统
    # 注意：电解液系统必须使用 PME 处理离子间的长程静电
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,  # 固定氢键以允许 2fs 步长
        rigidWater=True
    )
    
    # 添加压力计 (Barostat) 以实现 NPT 模拟
    system.addForce(openmm.MonteCarloBarostat(PRESSURE, TEMPERATURE))
    
    # 3. 设置积分器 (Langevin 热浴)
    integrator = openmm.LangevinMiddleIntegrator(TEMPERATURE, 1/unit.picosecond, TIMESTEP)
    
    # 4. 选择计算平台 (强制 CUDA)
    try:
        platform = openmm.Platform.getPlatformByName("CUDA")
        properties = {'Precision': 'mixed'}
        print(f"Using platform: CUDA (Acceleration Enabled)")
    except Exception:
        platform = openmm.Platform.getPlatformByName("Reference")
        properties = {}
        print("Warning: CUDA not found, falling back to CPU (Slow!)")

    # 5. 初始化仿真对象
    simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)
    simulation.context.setPositions(pdb.positions)
    
    # 6. 能量最小化
    print("Step 1: Minimizing energy...")
    simulation.minimizeEnergy()
    
    # 7. 预平衡 (NVT -> NPT)
    print("Step 2: Equilibration (100 ps)...")
    simulation.step(50000) 
    
    # 8. 设置数据记录器
    simulation.reporters.append(app.DCDReporter(OUTPUT_TRAJ, REPORT_STEPS))
    simulation.reporters.append(app.StateDataReporter(
        OUTPUT_LOG, REPORT_STEPS, step=True, potentialEnergy=True, 
        temperature=True, density=True, speed=True, volume=True
    ))
    # 屏幕打印
    simulation.reporters.append(app.StateDataReporter(
        sys.stdout, REPORT_STEPS, step=True, potentialEnergy=True, density=True, speed=True
    ))
    
    # 9. 正式生产运行
    print(f"Step 3: Production Run ({TOTAL_STEPS} steps)...")
    simulation.step(TOTAL_STEPS)
    
    # 10. 保存最终状态
    simulation.saveCheckpoint(CHECKPOINT)
    print("--- Simulation Completed Successfully ---")

if __name__ == "__main__":
    run_npt_simulation()
