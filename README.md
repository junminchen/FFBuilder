# FFBuilder

> Automated OPLS-AA force field generation for molecular dynamics — with **BOSS-free** (web LigParGen) and **BOSS-dependent** (local LigParGen, RESP) pipelines.

---

## Project Structure

```
FFBuilder/
├── README.md
├── .gitignore
│
├── data/                      ← All data: molecules, force fields, examples
│   ├── mol_db.py              ← 107 pre-defined electrolyte molecules
│   ├── molecules/             ← Per-molecule working dirs (see Output below)
│   ├── forcefields/           ← Base OPLS-AA XML files
│   ├── compare_FAN/           ← FAN(FCC#N) comparison results
│   ├── examples/              ← Electrolyte bulk system example
│   └── ETH.mol
│
├── web/                       ← BOSS-free pipeline (network required)
│   └── MD_OPLS_Workflow.py    ← traken.chem.yale.edu LigParGen web scraper
│
├── native/                    ← BOSS-dependent pipelines
│   ├── ligpargen_local/
│   │   ├── LigParGen_2.3/     ← LigParGen v2.3 Python package
│   │   ├── boss0824.tar.gz     ← BOSS 5.1 distribution (extract to ~/BOSS/)
│   │   └── AutoFF_LPG_Builder.py
│   └── resp/
│       ├── RESP_Workflow.py   ← RESP/AM1-BCC charge workflow
│       └── AutoFF_Builder.py  ← RESP batch builder
│
├── ht_workflow/               ← Active learning screening (GPR + EI)
│   ├── screening.py            ← Morgan FP + PCA + GPR + Expected Improvement
│   ├── ht_workflow.py         ← Batch orchestration
│   ├── ht_utils.py            ← Simulation result checker (template)
│   ├── ligpargen_local.py     ← Local LigParGen wrapper
│   └── run_ht_screening.py    ← CLI entry point
│
└── utils/                      ← Shared utilities
    ├── FFutils.py             ← Charge detection, XML generation, PDB skeleton
    ├── forcefield.py          ← ForceFieldManager (OpenMM XML loading/merging)
    ├── build_electrolyte.py   ← Electrolyte system builder
    ├── run_MD_bulk.py         ← MD batch runner
    ├── compare_charges.py     ← Charge comparison tool
    ├── generate_parity_svg.py
    ├── update_opls_resp.py
    └── mol_list.py
```

---

## Quick Start

```bash
# Clone
git clone https://github.com/junminchen/FFBuilder.git
cd FFBuilder

# Environment
conda create -n ffbuilder -c rdkit -c conda-forge rdkit openbabel pandas numpy requests
conda activate ffbuilder

# Install
pip install -e ./native/ligpargen_local/LigParGen_2.3
```

---

## Three Pipelines

| Pipeline | Charge | BOSS | Network | Output |
|----------|--------|------|---------|--------|
| **Web LigParGen** | CM1A | ❌ | ✅ | Per-molecule XML |
| **Local LigParGen** | CM1A / CM1A-LBCC / CM5 | ✅ | After setup | Per-molecule XML |
| **RESP** | RESP / AM1-BCC | ❌ | ❌ | Per-molecule XML |

---

## Pipeline 1: Web LigParGen (BOSS-free)

> Network access required. CM1A charges generated server-side at Yale.

```python
from web.MD_OPLS_Workflow import OPLSWorkflow

wf = OPLSWorkflow({'FAN': 'FCC#N', 'FEC': 'FCC1OC1O1'})
wf.run()
```

Or from `mol_db.py`:

```python
from web.MD_OPLS_Workflow import OPLSWorkflow
from data.mol_db import mol_db

wf = OPLSWorkflow(mol_db)
wf.run()
```

Or CLI:

```bash
# Edit utils/mol_list.py first to set your molecules
python web/MD_OPLS_Workflow.py
```

### Output

Each molecule produces a self-contained OpenMM force field XML:

```
data/molecules/<NAME>/
├── <NAME>.xml              ← Raw LigParGen XML (server output)
├── <NAME>.pdb              ← Raw LigParGen PDB
├── monomer.pdb               ← PDB with element column + canonical atom names
└── monomer_<NAME>/
    └── ff.xml                 ← ✅ OpenMM-compatible OPLS-AA XML (use this)
```

**To load in OpenMM:**

```python
from openmm.app import ForceField
ff = ForceField('data/molecules/FAN/monomer_FAN/ff.xml')
```

**To merge multiple molecules in OpenMM at runtime:**

```python
from openmm.app import ForceField
ff1 = ForceField('data/molecules/FAN/monomer_FAN/ff.xml')
ff2 = ForceField('data/molecules/FEC/monomer_FEC/ff.xml')
# OpenMM does not natively support merging — use utils/forcefield.py:
from utils.forcefield import ForceFieldManager
mgr = ForceFieldManager()
mgr.load_ff('data/molecules/FAN/monomer_FAN/ff.xml')
mgr.load_ff('data/molecules/FEC/monomer_FEC/ff.xml')
mgr.save('data/forcefields/custom_combined.xml')
```

### Input Format

`mol_dict` is a plain `dict`:

```python
{'NAME': 'SMILES', ...}
# e.g. {'FAN': 'FCC#N', 'FEC': 'FCC1OC1O1', 'Li': '[Li+]'}
```

For salt molecules, use `[Li+]` (cation) or `[Fluorine-]`/`[BF4-]` etc. (anion) with net charge automatically detected.

---

## Pipeline 2: Local LigParGen (BOSS required)

### BOSS Setup

BOSS is a Fortran/C quantum chemistry program (Yale Jorgensen lab). Binaries are **32-bit Linux ELF** — requires Linux or Docker.

**Option A: Docker**

```bash
cd native/ligpargen_local
docker build -t ligpargen .
docker run -v $(pwd)/../../../data:/data ligpargen bash -c "export BOSSdir=/root/BOSS/boss && python native/ligpargen_local/AutoFF_LPG_Builder.py"
```

**Option B: Remote SSH**

Point to an HPC with BOSS installed:

```python
import subprocess

def run_ligpargen_remote(smiles, name, res_name='UNK', charge=0):
    cmd = f"export BOSSdir=~/BOSS/boss && LigParGen -s '{smiles}' -r {res_name} -c {charge} -o 0 -l"
    result = subprocess.run(['ssh', 'your-hpc', cmd], capture_output=True)
    return result.returncode == 0
```

### Usage

```bash
export BOSSdir=~/BOSS/boss

# Single molecule
python -c "
import sys; sys.path.insert(0, '.')
from ht_workflow.ligpargen_local import run_ligpargen_local
ok = run_ligpargen_local('FCC#N', 'FAN', res_name='FAN', charge=0)
print('Success:', ok)
"

# Batch
python native/ligpargen_local/AutoFF_LPG_Builder.py --limit 10
```

---

## Pipeline 3: RESP (BOSS-free, custom QM charges)

Uses RESP or AM1-BCC charges from your own QM calculations. No BOSS needed.

**Requirements:** `rdkit`, `psi4` (optional, for on-the-fly QM), OR pre-computed charge files.

```python
from native.resp.RESP_Workflow import load_qm_charges, process_molecule, assign_charges_to_mol
from utils.FFutils import generate_ff_xml

# Load charges from ORCA/Psi4 output
charges = load_qm_charges('FAN_orca.chg', method='RESP')
symbols, coords, _, rd_mol = process_molecule('FAN', 'FCC#N', charge=0)
atom_charges = assign_charges_to_mol(rd_mol, charges)
generate_ff_xml('FAN', symbols, coords, atom_charges, rd_mol,
                 'data/molecules/FAN/monomer_FAN/ff.xml')
```

Or batch:

```bash
python native/resp/AutoFF_Builder.py --limit 10
```

---

## Active Learning Workflow

```bash
# 1. Prepare candidates
echo "CCCC\nCCCCCO\nc1ccccc1" > candidates.txt

# 2. (Optional) Training data: smiles + property to predict
echo "smiles,property" > training_data.csv
echo "CCO,0.5" >> training_data.csv

# 3. Run one iteration
python ht_workflow/run_ht_screening.py -c candidates.txt -t training_data.csv -n 3
# → batches/batch0/suggested_smiles.txt
# → data/molecules/<NAME>/monomer_<NAME>/ff.xml
```

Active learning uses **GPR + Expected Improvement** over Morgan fingerprints (2048-bit, radius=2) with PCA → 50D reduction.

⚠️ **Single-shot by default.** Multi-round iteration requires implementing `ht_utils.py` to read simulation results back into `training_data.csv`.

---

## Molecular Database

```python
from data.mol_db import mol_db, salts, solvents
print(f"Total: {len(mol_db)} molecules ({len(salts)} salts + {len(solvents)} solvents)")
```

---

## Force Field XML Format

Both pipelines produce OpenMM-compatible XML:

```xml
<ForceField>
  <AtomTypes>
    <Type name="FAN_1" class="FAN_1" element="F" mass="19.00"/>
  </AtomTypes>
  <Residues>
    <Residue name="FAN">
      <Atom name="F00" type="FAN_1" charge="-0.147000"/>
    </Residue>
  </Residues>
  <HarmonicBondForce>...</HarmonicBondForce>
  <HarmonicAngleForce>...</HarmonicAngleForce>
  <PeriodicTorsionForce>...</PeriodicTorsionForce>
  <NonbondedForce coulomb14scale="0.5" lj14scale="0.5">...</NonbondedForce>
</ForceField>
```

---

## Comparing Web vs Local LigParGen

| Aspect | Web LigParGen | Local LigParGen |
|--------|--------------|-----------------|
| Charge | CM1A (server) | CM1A / CM1A-LBCC / CM5 |
| CM5 charges | ❌ | ✅ (`-q QORCA`) |
| Network | Required | Not after setup |
| Batch | Rate-limited | Unlimited |
| Speed | ~5-10s/mol | ~5s/mol |

---

## Known Issues

### traken.chem.yale.edu parsing
The server uses separate `<input type="hidden">` (fileout) and `<input type="submit">` (go) tags — not grouped inside `<form>` elements. Regex `zip(go_vals, fileout_vals)` is the reliable parsing method.

### BOSS on macOS
BOSS 5.1 binaries are **32-bit Linux ELF**. Cannot run natively on macOS (Apple Silicon or Intel). Use Docker or a Linux VM/HPC.

---

## Citation

- **LigParGen**: Dodda et al., *Nucleic Acids Research* 2017, W331-W336
- **CM1A-LBCC**: Dodda et al., *J. Phys. Chem. B* 2017, 121, 3864-3870
- **OPLS-AA**: Jorgensen et al., *J. Am. Chem. Soc.* 1996, 118, 11225-11236
