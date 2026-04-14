#!/usr/bin/env python
import os
import sys
import argparse
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

try:
    from data.mol_db import mol_db
except ImportError:
    from mol_db import mol_db

from utils.FFutils import get_molecule_charge, ensure_mol_dir
from ht_workflow.ligpargen_local import run_ligpargen_local, extract_parameters_from_lpg_xml

def build_all_lpg(limit=10, base_ff="opls_solvent.xml"):
    update_ff_path = os.path.join("data", "forcefields", "opls_lpg_update.xml")
    base_ff_path = os.path.join("data", "forcefields", base_ff)
    
    # Load base or existing update
    source_path = update_ff_path if os.path.exists(update_ff_path) else base_ff_path
    if os.path.exists(source_path):
        tree = ET.parse(source_path)
        root = tree.getroot()
    else:
        root = ET.Element("ForceField")
        ET.SubElement(root, "AtomTypes")
        ET.SubElement(root, "Residues")
        
    atom_types_node = root.find("AtomTypes")
    residues_node = root.find("Residues")
    
    existing_residues = [res.get("name") for res in residues_node.findall("Residue")]
    
    count = 0
    for name, smiles in mol_db.items():
        if count >= limit:
            break
            
        res_name = name.upper()
        if res_name in existing_residues:
            print(f"Skipping {name}: already in forcefield.")
            continue
            
        print(f"\n--- LigParGen processing: {name} ---")
        
        # Determine total charge
        charge = get_molecule_charge(name, smiles)
        out_dir = ensure_mol_dir(name)
        
        # We use a short name for residue in LPG to keep it 3 chars if possible
        lpg_res_name = name[:3].upper()
        
        try:
            if run_ligpargen_local(smiles, name, work_dir=out_dir, charge=charge, res_name=lpg_res_name):
                xml_path = os.path.join(out_dir, f"{lpg_res_name}.xml")
                if os.path.exists(xml_path):
                    params = extract_parameters_from_lpg_xml(xml_path, name)
                    if params:
                        # Append to global
                        for at in params["AtomTypes"]:
                            atom_types_node.append(at)
                        for res in params["Residues"]:
                            residues_node.append(res)
                        for tag, f_node in params["Forces"].items():
                            target_f = root.find(tag)
                            if target_f is None:
                                target_f = ET.SubElement(root, tag)
                                if tag == "NonbondedForce":
                                    target_f.attrib = f_node.attrib
                            for item in f_node:
                                target_f.append(item)
                        
                        print(f"Successfully integrated {name} into global forcefield.")
                        count += 1
        except Exception as e:
            print(f"Error processing {name}: {e}")

    # Save merged XML
    xml_str = minidom.parseString(ET.tostring(root)).toprettyxml(indent="  ")
    clean_xml = "\n".join([line for line in xml_str.split('\n') if line.strip()])
    
    with open(update_ff_path, "w") as f:
        f.write(clean_xml)
    print(f"\nDone! Integrated forcefield saved to {update_ff_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build complete OPLS-AA forcefields using LigParGen for mol_db")
    parser.add_argument("-l", "--limit", type=int, default=5, help="Limit number of molecules to process")
    args = parser.parse_args()
    
    build_all_lpg(limit=args.limit)
