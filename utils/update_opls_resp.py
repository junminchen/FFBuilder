import os
import sys
import numpy as np
import xml.etree.ElementTree as ET

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from data.mol_db import mol_db
except ImportError:
    from mol_db import mol_db

from native.resp.RESP_Workflow import process_molecule
from utils.FFutils import get_molecule_charge

def update_xml_charges(in_path, out_path):
    print(f"\n--- Updating {in_path} -> {out_path} ---")
    if not os.path.exists(in_path):
        print(f"Error: {in_path} not found.")
        return
        
    tree = ET.parse(in_path)
    root = tree.getroot()
    
    # 映射 Residue 名字到 mol_db 中的名字 (处理命名差异)
    name_map = {
        'LiA': 'Li',
        'NaA': 'Na',
        'DFO': 'DFOB',
        'TFS': 'TFSI',
        'ECA': 'EC',
        'PCA': 'PC',
        'PPA': 'PP',
        'PSA': 'PS',
        'EPA': 'EP',
        'DFE': 'DFEA',
        'FEM': 'FEMC',
        'SLA': 'SL',
        'HFD': 'HFDEC',
        'TER': 'TTE', # 可能的对应关系
    }
    
    residues = root.findall(".//Residue")
    
    for res in residues:
        res_name = res.get("name")
        target_name = name_map.get(res_name, res_name)
        
        if target_name in mol_db:
            smiles = mol_db[target_name]
            total_charge = get_molecule_charge(target_name, smiles)
            
            try:
                # 重新计算 RESP 电荷 (SCALE=1.0)
                symbols, coords, charges, _ = process_molecule(target_name, smiles, total_charge)
                
                # 对应修改 XML 中的 Atom charge
                atoms = res.findall("Atom")
                if len(atoms) != len(charges):
                    print(f"Warning: Atom count mismatch for {res_name} ({len(atoms)} in XML vs {len(charges)} in calculation). Skipping.")
                    continue
                
                for i, atom in enumerate(atoms):
                    atom.set("charge", str(round(charges[i], 4)))
                print(f"Successfully updated charges for {res_name} using RESP.")
                
            except Exception as e:
                print(f"Failed to update {res_name}: {e}")
        else:
            print(f"No SMILES found for {target_name}, skipping.")

    # 格式化并保存
    ET.indent(tree, space="    ")
    tree.write(out_path, encoding="utf-8", xml_declaration=True)
    print(f"New forcefield file saved to {out_path}")

if __name__ == "__main__":
    update_xml_charges("forcefields/opls_salt.xml", "forcefields/opls_salt_resp.xml")
    update_xml_charges("forcefields/opls_solvent.xml", "forcefields/opls_solvent_resp.xml")
