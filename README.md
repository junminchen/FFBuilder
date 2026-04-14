# FFBuilder

> Automated Force Field Builder for OPLS-AA MD simulations — with **BOSS-dependent** (LigParGen local) and **BOSS-free** (Web LigParGen, RESP) pipelines.

---

## Project Structure

```
FFBuilder/
├── README.md              ← This file
├── .gitignore
├── handoff.md             ← Legacy notes
│
├── data/                  ← Molecule & force field data
│   ├── mol_db.py          ← 107 electrolyte molecules (salts + solvents)
│   ├── molecules/         ← .pdb / .mol structure files
│   ├── forcefields/        ← Output XML force field files
│   └── compare_FAN/        ← FAN comparison results (Web LigParGen)
│
├── utils/                  ← Shared utilities
│   ├── FFutils.py         ← Charge detection, XML generation, PDB skeleton
│   ├── forcefield.py      ← ForceFieldManager class
│   ├── build_electrolyte.py
│   ├── run_MD_bulk.py
│   ├── compare_charges.py
│   ├── generate_parity_svg.py
│   ├── update_opls_resp.py
│   └── mol_list.py
│
├── native/                 ← BOSS-dependent (local LigParGen + RESP)
│   ├── ligpargen_local/
│   │   ├── LigParGen_2.3/  ← LigParGen v2.3 Python package
│   │   └── AutoFF_LPG_Builder.py
│   └── resp/
│       ├── RESP_Workflow.py
│       └── AutoFF_Builder.py
│
├── web/                    ← BOSS-free (traken.chem.yale.edu)
│   └── MD_OPLS_Workflow.py
│
└── ht_workflow/            ← Active learning screening (GPR + EI)
    ├── screening.py         ← Morgan FP + PCA + GPR + Expected Improvement
    ├── ht_workflow.py      ← Batch orchestration
    ├── ht_utils.py         ← Simulation result checker (template)
    ├── ligpargen_local.py   ← Local LigParGen wrapper
    └── run_ht_screening.py ← CLI entry point
```

---

## Quick Start

```bash
git clone https://github.com/junminchen/FFBuilder.git
cd FFBuilder

# Create environment
conda create -n ffbuilder -c rdkit -c conda-forge rdkit openbabel pandas numpy
conda activate ffbuilder

# Install
pip install -e ./native/ligpargen_local/LigParGen_2.3
```

---

## Three Pipelines

| Pipeline | Charge Method | BOSS Required | Output |
|----------|--------------|---------------|--------|
| **Web LigParGen** | CM1A (server-side) | ❌ No | Per-molecule XML |
| **LigParGen Local** | CM1A / CM1A-LBCC | ✅ Yes | Per-molecule XML |
| **RESP** | RESP / AM1-BCC | ❌ No | Per-molecule XML |

---

### Pipeline 1: Web LigParGen (BOSS-free)

Uses Yale's web server — no BOSS needed, just network access.

```bash
# Single molecule
python web/MD_OPLS_Workflow.py
# → data/compare_FAN/web/molecules/FAN/FAN.xml

# Batch from mol_db
python -c "
from web.MD_OPLS_Workflow import OPLSWorkflow
from data.mol_db import mol_db
wf = OPLSWorkflow(mol_db)
wf.run()
"
```

**Requirements:** `requests`, `beautifulsoup4`, `MDAnalysis`, `rdkit`

> ⚠️ Known issue: traken.chem.yale.edu HTML uses double-quote attributes (not escaped). Use raw regex extraction:
> ```python
> import re
> fileouts = re.findall(r'name="fileout" value="([^"]+)"', res.text)
> ```

---

### Pipeline 2: LigParGen Local (BOSS required)

Runs LigParGen 2.3 locally. **BOSS software required.**

#### BOSS Installation

```bash
# BOSS is a Fortran/C quantum chemistry program (Yale Jorgensen lab)
# Binaries are 32-bit Linux ELF → requires Linux or Docker

# Option A: Docker
docker run -it -v $(pwd):/work i386/ubuntu:18.04
apt-get install -y csh openjdk python3
# Install BOSS inside container

# Option B: Remote SSH
# Point to your HPC with BOSS installed
```

#### Usage

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

### Pipeline 3: RESP (BOSS-free, custom QM charges)

Uses RESP/AM1-BCC charges from your own QM calculations. No BOSS needed.

```bash
# With external QM charge files
python -c "
from native.resp.RESP_Workflow import load_qm_charges, process_molecule, assign_charges_to_mol
from utils.FFutils import generate_ff_xml

charges = load_qm_charges('FAN_orca.chg', method='RESP')
symbols, coords, _, rd_mol = process_molecule('FAN', 'FCC#N', charge=0)
atom_charges = assign_charges_to_mol(rd_mol, charges)
generate_ff_xml('FAN', symbols, coords, atom_charges, rd_mol, 'FAN.xml')
"
```

---

## Active Learning Workflow

```bash
# 1. Prepare candidates
echo "CCCC\nCCCCCO\nc1ccccc1" > candidates.txt

# 2. (Optional) Provide training data
echo "smiles,property" > training_data.csv
echo "CCO,0.5" >> training_data.csv

# 3. Run one iteration
python ht_workflow/run_ht_screening.py -c candidates.txt -t training_data.csv -n 3
# → batches/batch0/suggested_smiles.txt
# → data/forcefields/opls_ht_update.xml
```

Active learning uses **GPR + Expected Improvement** over Morgan fingerprints (2048-bit, radius=2) with PCA → 50D reduction. Currently a **single-shot** loop — multi-round iteration requires implementing `ht_utils.py` to read simulation results back into `training_data.csv`.

---

## Molecular Database

```python
from data.mol_db import mol_db, salts, solvents
print(f"Total: {len(mol_db)} molecules ({len(salts)} salts + {len(solvents)} solvents)")
```

SMILES use RDKit atom indices (`[C:1]`, `[O:2]`) for charge assignment.

---

## Force Field Output

Both pipelines produce **OpenMM-compatible OPLS-AA XML**:

```xml
<ForceField>
  <AtomTypes>
    <Type name="FAN_1" class="FAN_1" element="F" mass="19.00"/>
  </AtomTypes>
  <Residues>
    <Residue name="FAN">
      <Atom name="F1" type="FAN_1" charge="-0.1325"/>
    </Residue>
  </Residues>
  <HarmonicBondForce>...</HarmonicBondForce>
  <HarmonicAngleForce>...</HarmonicAngleForce>
  <PeriodicTorsionForce>...</PeriodicTorsionForce>
  <NonbondedForce>...</NonbondedForce>
</ForceField>
```

---

## Comparing Web vs Local LigParGen

| Aspect | Web LigParGen | Local LigParGen |
|--------|--------------|-----------------|
| Charge method | CM1A (server-side) | CM1A / CM1A-LBCC (BOSS) |
| CM5 charges | ❌ | ✅ Via ORCA (`-q QORCA`) |
| Network | Required | Not after setup |
| Batch size | Rate-limited | Unlimited |

---

## Citation

- **LigParGen**: Dodda et al., *Nucleic Acids Research* 2017, W331-W336
- **CM1A-LBCC**: Dodda et al., *J. Phys. Chem. B* 2017, 121, 3864-3870
- **OPLS-AA**: Jorgensen et al., *J. Am. Chem. Soc.* 1996, 118, 11225-11236
