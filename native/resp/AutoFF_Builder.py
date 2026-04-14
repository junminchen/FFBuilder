import os
import sys

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

try:
    from data.mol_db import mol_db
except ImportError:
    from mol_db import mol_db

from native.resp.RESP_Workflow import process_molecule
from utils.FFutils import get_molecule_charge, generate_ff_xml, ensure_mol_dir

def run_workflow(limit=10):
    count = 0
    for name, smiles in mol_db.items():
        if count >= limit:
            break
            
        print(f"\n--- Batch processing: {name} ---")
        
        # Determine total charge from FFutils
        charge = get_molecule_charge(name, smiles)
            
        try:
            # Step 1: Process molecule (Structure + Charges)
            symbols, coords, charges, rd_mol = process_molecule(name, smiles, charge)
            
            # Step 2: Ensure directory exists from FFutils
            out_dir = ensure_mol_dir(name)
            
            # Step 3: Generate XML from FFutils
            xml_path = f"{out_dir}/{name}.xml"
            generate_ff_xml(name, symbols, coords, charges, rd_mol, xml_path)
            
            count += 1
        except Exception as e:
            print(f"Error processing {name}: {e}")

if __name__ == "__main__":
    run_workflow(limit=10)
