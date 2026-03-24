# FFBuilder Handoff - High-Throughput & LigParGen Integration

This document summarizes the integration of High-Throughput (HT) screening and local LigParGen workflows into the FFBuilder repository.

## 1. New Components Added

### Active Learning Screening (`ht_screening/`)
- **`screening.py`**: A `torch`-based active learning module.
  - Implements **Gaussian Process Regressor (GPR)**.
  - Uses **RDKit Morgan Fingerprints** for feature extraction.
  - Employs **Expected Improvement (EI)** to suggest the next batch of molecules.
- **`ligpargen_local.py`**: Local automation for LigParGen.
  - Supports generating polymer chains from monomers (`[Cu]`/`[Au]` logic).
  - Handles local CLI calls to `LigParGen` and patches/extracts OpenMM XML parameters.
- **`ht_workflow.py`**: The main orchestrator for batch-based high-throughput runs.
  - Automatically manages batch folders (`batches/batchX`).
  - Merges new parameters into a global `forcefields/opls_ht_update.xml`.
- **`ht_utils.py`**: Utility functions for tracking simulation success across batches.

### Automation Scripts (Root)
- **`run_ht_screening.py`**: CLI entry point for the active learning workflow.
- **`AutoFF_LPG_Builder.py`**: A batch builder that uses LigParGen to generate full OPLS-AA force fields for all molecules in `mol_db.py`.

## 2. Integration with Existing Architecture
- **`mol_db.py` Integration**: `HighThroughputWorkflow` and `AutoFF_LPG_Builder.py` now directly pull SMILES from the unified `mol_db.py`.
- **`FFutils.py` Integration**: New scripts use the standard `get_molecule_charge` and `ensure_mol_dir` from your latest architecture.
- **Force Field Merging**: The system is designed to be additive. New molecules generated via LigParGen are merged into `forcefields/opls_lpg_update.xml` without overwriting base force fields.

## 3. Setup & Requirements
- **LigParGen 2.3**: A patched version is included in `LigParGen_2.3/` (fixed for `networkx` 3.x compatibility).
- **Environment Variables**:
  - `BOSSdir`: Must point to your local BOSS installation for LigParGen to calculate CM1A charges.
- **Dependencies**:
  - `torch`, `rdkit`, `numpy`, `pandas`, `openbabel-wheel`.

## 4. Usage Examples

### Run Active Learning Iteration
```bash
python run_ht_screening.py --nsuggest 5 --repeats 1
```

### Build Full Forcefields for mol_db
```bash
python AutoFF_LPG_Builder.py --limit 10
```

### Compare Charges (Existing Tool)
```bash
python generate_parity_svg.py
```

## 5. Summary of Patches
- **LigParGen**: Patched `CreatZmat.py` to use `G.nodes[]` instead of the deprecated `G.node[]` to support Python 3.13 and modern NetworkX.
