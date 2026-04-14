# FFBuilder

> Automated Force Field Builder for OPLS-AA — with support for both **LigParGen/BOSS** (high-accuracy CM1A charges) and **native RESP** (BOSS-free, requires your own QM calculations).

---

## Project Overview

FFBuilder generates OpenMM-compatible OPLS-AA force field XML files for organic molecules (electrolytes, solvents, salts). It has **two independent pipelines**:

| Pipeline | Charge Method | BOSS Required | Requires |
|---------|--------------|---------------|----------|
| **LigParGen** | CM1A / CM1A-LBCC | Yes | BOSS software |
| **Web LigParGen** | CM1A / CM1A-LBCC | No (web service) | Network access |
| **RESP** | RESP / AM1-BCC | No | QMD calculation + external tool |

The two pipelines **output compatible XML formats** — both can be merged into the same global `forcefields/*.xml`.

---

## Quick Start

```bash
# 1. Clone
git clone https://github.com/junminchen/FFBuilder.git
cd FFBuilder

# 2. Create conda environment
conda create -n ffbuilder -c rdkit -c conda-forge rdkit openbabel pandas numpy
conda activate ffbuilder

# 3. Install LigParGen 2.3 (needed for both LigParGen and RESP paths)
cd LigParGen_2.3 && pip install -e . && cd ..

# 4. Install FFBuilder
pip install -e .
```

---

## Pipeline 1: Web LigParGen — No BOSS Required

Uses the Yale LigParGen web server (`traken.chem.yale.edu`). CM1A charges generated server-side.

**Requirements:** Network access, `requests`, `beautifulsoup4`, `MDAnalysis`, `rdkit`

```bash
pip install requests beautifulsoup4 MDAnalysis
```

### Single Molecule

```python
from MD_OPLS_Workflow import OPLSWorkflow

wf = OPLSWorkflow({'FAN': 'FCC#N'})
wf.mol_root = 'molecules'
wf.run()  # downloads → fixes PDB → consolidates FF → merges to forcefields/opls_solvent_update.xml
```

### Batch from mol_db.py

```python
from mol_db import mol_db
from MD_OPLS_Workflow import OPLSWorkflow

wf = OPLSWorkflow(mol_db)
wf.run()
```

### CLI

```bash
python MD_OPLS_Workflow.py
```

---

## Pipeline 2: LigParGen Local — BOSS Required

Runs LigParGen 2.3 locally. Requires **BOSS software** installed and `$BOSSdir` set.

### BOSS Installation

BOSS is a Fortran/C quantum chemistry program from the Jorgensen lab at Yale.

```bash
# 1. Obtain boss0824.tar.gz from Yale (or your existing copy)

# 2. Extract to ~/BOSS/
mkdir -p ~/BOSS
tar -xzf boss0824.tar.gz -C ~/BOSS/

# 3. Verify
export BOSSdir=~/BOSS/boss
ls $BOSSdir/scripts/xSP   # should exist

# 4. Add to shell profile (macOS/Linux)
echo 'export BOSSdir=~/BOSS/boss' >> ~/.zshrc
echo 'export PATH=$BOSSdir/scripts:$PATH' >> ~/.zshrc
source ~/.zshrc
```

**Note:** BOSS binaries are 32-bit Linux ELF. They require:
- Linux with 32-bit compatibility libraries, **or**
- A Docker/Linux VM, **or**
- Remote HPC execution via SSH

#### Docker Alternative

```dockerfile
FROM i386/ubuntu:18.04
RUN apt-get update && apt-get install -y csh openjdk python3 python3-pip
COPY boss0824.tar.gz /tmp/
RUN mkdir -p /root/BOSS && tar -xzf /tmp/boss0824.tar.gz -C /root/BOSS/
ENV BOSSdir=/root/BOSS/boss
```

#### Remote SSH Alternative

If BOSS is installed on a remote Linux machine (`your-hpc`):

```python
import subprocess

def run_ligpargen_remote(smiles, name, resname, work_dir):
    cmd = f"""
    export BOSSdir=~/BOSS/boss
    cd {work_dir}
    LigParGen -s '{smiles}' -r {resname} -c 0 -o 0 -l
    """
    result = subprocess.run(
        ['ssh', 'your-hpc', cmd],
        capture_output=True, text=True
    )
    return result.returncode == 0
```

### Usage

```bash
# Single molecule
python -c "
from ht_screening.ligpargen_local import run_ligpargen_local, extract_parameters_from_lpg_xml
ok = run_ligpargen_local('FCC#N', 'FAN', res_name='FAN', charge=0)
print('Success:', ok)
"

# Batch from mol_db.py
python AutoFF_LPG_Builder.py --limit 10

# With Active Learning workflow
python run_ht_screening.py -c candidates.txt -n 5
```

---

## Pipeline 3: RESP — BOSS-Free (Custom QM Charges)

Uses RESP/AM1-BCC charges from external QM calculations. No BOSS needed.

**Requirements:** `rdkit`, `numpy`; optionally `Psi4` or `ORCA` for QM calculations.

### Setup

```bash
# Install Psi4 for quantum calculations (optional, if not using pre-computed charges)
pip install psi4
```

### Usage

```bash
# Process molecules with pre-computed charges
python AutoFF_Builder.py --limit 10

# Or use the RESP workflow directly
python -c "
from RESP_Workflow import process_molecule
from FFutils import generate_ff_xml

symbols, coords, charges, rd_mol = process_molecule('FAN', 'FCC#N', charge=0)
generate_ff_xml('FAN', symbols, coords, charges, rd_mol, 'FAN.xml')
"
```

### Custom QM → RESP Workflow

If you have ORCA or Psi4 charge files:

```python
from RESP_Workflow import load_qm_charges, assign_charges_to_mol

# Load charges from external QM output
charges = load_qm_charges('FAN_orca.chg', method='RESP')
symbols, coords, _, rd_mol = process_molecule('FAN', 'FCC#N', charge=0)
atom_charges = assign_charges_to_mol(rd_mol, charges)
generate_ff_xml('FAN', symbols, coords, atom_charges, rd_mol, 'FAN.xml')
```

---

## Directory Structure

```
FFBuilder/
├── mol_db.py              # Unified molecule database (salts + solvents)
├── FFutils.py             # Shared utilities (charge, XML generation)
│
├── MD_OPLS_Workflow.py    # Web LigParGen pipeline (Stage 1-4)
├── AutoFF_LPG_Builder.py  # Local LigParGen batch builder
├── AutoFF_Builder.py     # RESP-based FF builder
│
├── ht_screening/         # Active learning screening module
│   ├── screening.py       # GPR + Expected Improvement (PyTorch)
│   ├── ligpargen_local.py # Local LigParGen wrapper
│   └── ht_workflow.py     # Batch + merge orchestrator
│
├── LigParGen_2.3/        # LigParGen 2.3 Python package (BOSS-dependent)
├── RESP_Workflow.py       # RESP charge workflow
│
├── forcefields/           # Output directory
│   ├── opls_solvent.xml   # Base OPLS-AA force field
│   ├── opls_ht_update.xml # HT/Active Learning updates
│   └── opls_lpg_update.xml # LigParGen local updates
│
└── molecules/             # Per-molecule working directories
```

---

## Comparing Web vs Local LigParGen

Both produce OPLS-AA XML files. Key differences:

| Aspect | Web LigParGen | Local LigParGen |
|--------|--------------|-----------------|
| Charge method | CM1A (server-side) | CM1A / CM1A-LBCC (BOSS) |
| CM5 charges | ❌ Not available | ✅ Via ORCA (`-q` flag) |
| Network | Required | Not after setup |
| Geometry optimization | Server BOSS | Local BOSS |
| Processing speed | ~5-10s/molecule (network + queue) | ~5s/molecule (local) |
| Batch size | Rate-limited | Unlimited |

---

## Molecular Database

`mol_db.py` contains ~107 pre-defined electrolyte molecules:

- **9 salts**: Li/Na BF4, FSI, TFSI, PF6, BOB, DFOB, DFP, etc.
- **~98 solvents**: EC, PC, DMC, EMC, FEC, VC, AN, DOL, DEC, and fluorinated analogues

SMILES use RDKit atom indices (`[C:1]`, `[O:2]`) for charge assignment.

```python
from mol_db import mol_db, salts, solvents
print(f"Total: {len(mol_db)} molecules ({len(salts)} salts + {len(solvents)} solvents)")
```

---

## Force Field Output

All pipelines produce OpenMM-compatible XML:

```xml
<ForceField>
  <AtomTypes>
    <Type name="FAN_1" class="FAN_1" element="F" mass="19.00"/>
    ...
  </AtomTypes>
  <Residues>
    <Residue name="FAN">
      <Atom name="F1" type="FAN_1" charge="-0.1325"/>
      ...
    </Residue>
  </Residues>
  <HarmonicBondForce>...</HarmonicBondForce>
  <HarmonicAngleForce>...</HarmonicAngleForce>
  <PeriodicTorsionForce>...</PeriodicTorsionForce>
  <NonbondedForce>...</NonbondedForce>
</ForceField>
```

---

## Active Learning Workflow (BOSS + LigParGen)

```bash
# 1. Prepare candidate SMILES
echo "FCC#N\nC1CCCCC1\nCC(=O)O" > candidates.txt

# 2. (Optional) Provide training data
echo "smiles,property" > training_data.csv
echo "CCO,0.5" >> training_data.csv

# 3. Run one iteration
python run_ht_screening.py -c candidates.txt -t training_data.csv -n 3

# 4. Check suggested molecules
cat batches/batch0/suggested_smiles.txt
```

Active learning uses GPR with Expected Improvement to suggest the most promising molecules for the next round of MD simulation.

---

## Citation

If you use FFBuilder, cite:

- **LigParGen**: Dodda et al., *Nucleic Acids Research* 2017, W331-W336
- **CM1A-LBCC**: Dodda et al., *J. Phys. Chem. B* 2017, 121, 3864-3870
- **OPLS-AA**: Jorgensen et al., *J. Am. Chem. Soc.* 1996, 118, 11225-11236
