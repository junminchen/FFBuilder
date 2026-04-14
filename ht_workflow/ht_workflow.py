import os
import time
import json
import argparse
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom
from typing import List, Dict, Any, Optional
import pandas as pd

from .ligpargen_local import create_long_smiles, run_ligpargen_local, extract_parameters_from_lpg_xml
from .screening import suggest_next
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from utils.forcefield import ForceFieldManager

try:
    from data.mol_db import mol_db
    MOL_DB_AVAILABLE = True
except ImportError:
    try:
        from mol_db import mol_db
        MOL_DB_AVAILABLE = True
    except ImportError:
        MOL_DB_AVAILABLE = False

class BatchManager:
    def __init__(self, root_dir: str = "batches"):
        self.root_dir = root_dir
        if not os.path.exists(root_dir):
            os.makedirs(root_dir)
            
    def get_latest_batch_id(self) -> int:
        batches = [d for d in os.listdir(self.root_dir) if d.startswith("batch")]
        if not batches:
            return -1
        ids = [int(d.replace("batch", "")) for d in batches if d.replace("batch", "").isdigit()]
        return max(ids) if ids else -1
        
    def create_next_batch(self) -> str:
        next_id = self.get_latest_batch_id() + 1
        batch_dir = os.path.join(self.root_dir, f"batch{next_id}")
        os.makedirs(batch_dir)
        return batch_dir

class HighThroughputWorkflow:
    def __init__(self, candidate_file: str = "candidates.txt", training_data_file: str = "training_data.csv", base_ff="opls_solvent.xml"):
        self.candidate_file = candidate_file
        self.training_data_file = training_data_file
        self.manager = ForceFieldManager()
        self.base_ff_path = os.path.join("forcefields", base_ff)
        self.update_ff_path = os.path.join("forcefields", "opls_ht_update.xml")
        self.batch_manager = BatchManager()
        
        # Load candidates from file or mol_db
        self.candidates = []
        if os.path.exists(candidate_file):
            with open(candidate_file, 'r') as f:
                self.candidates = [line.strip() for line in f if line.strip()]
        
        if not self.candidates and MOL_DB_AVAILABLE:
            print("No candidates file found, loading from mol_db.py")
            self.candidates = list(mol_db.values())
            
        # Load training data if available
        self.train_smiles = []
        self.train_props = []
        self._load_training_data()

    def _load_training_data(self):
        if self.training_data_file and os.path.exists(self.training_data_file):
            try:
                df = pd.read_csv(self.training_data_file)
                if 'smiles' in df.columns and 'property' in df.columns:
                    self.train_smiles = df['smiles'].tolist()
                    self.train_props = df['property'].tolist()
            except Exception as e:
                print(f"Error loading training data: {e}")

    def run_iteration(self, n_suggest=5, repeats=1, batch_dir: Optional[str] = None):
        print(f"--- Starting HT Iteration ---")
        
        if not batch_dir:
            batch_dir = self.batch_manager.create_next_batch()
        
        # 1. Suggest next molecules
        available_candidates = [s for s in self.candidates if s not in self.train_smiles]
        
        if not available_candidates:
            print("No more candidates to screen.")
            return
            
        suggested = suggest_next(available_candidates, self.train_smiles, self.train_props, n_suggest=n_suggest)
        print(f"Suggested molecules: {suggested}")
        
        # Save suggested to batch dir
        with open(os.path.join(batch_dir, "suggested_smiles.txt"), "w") as f:
            for s in suggested:
                f.write(f"{s}\n")
        
        # 2. Process each suggested molecule
        new_params = []
        batch_id = self.batch_manager.get_latest_batch_id()
        for i, smiles in enumerate(suggested):
            # Create a 3-letter residue name based on batch and index
            # e.g., B00, B01... for batch 0
            res_name = f"B{batch_id:01d}{i:01d}"[:3].upper()
            name = f"HT_Batch{batch_id}_{res_name}_{i}"
            print(f"Processing {smiles} as {name} (Residue: {res_name})...")
            
            # Generate polymer if requested
            if repeats > 1:
                proc_smiles, _ = create_long_smiles(smiles, repeats=repeats)
            else:
                proc_smiles = smiles
                
            mol_dir = os.path.join(batch_dir, name)
            if run_ligpargen_local(proc_smiles, name, work_dir=mol_dir, res_name=res_name):
                # The actual residue name in the XML will be res_name
                xml_path = os.path.join(mol_dir, f"{res_name}.xml")
                if os.path.exists(xml_path):
                    params = extract_parameters_from_lpg_xml(xml_path, res_name)
                    if params:
                        new_params.append(params)
                        print(f"Successfully extracted parameters for {res_name}")
            else:
                print(f"Failed to generate parameters for {name}")
                
        # 3. Merge new parameters into global force field
        if new_params:
            self.merge_to_global_ff(new_params)
            
    def merge_to_global_ff(self, new_params_list: List[Dict[str, Any]]):
        print(f"Merging {len(new_params_list)} new molecules into {self.update_ff_path}...")
        
        source_path = self.update_ff_path if os.path.exists(self.update_ff_path) else self.base_ff_path
        
        if os.path.exists(source_path):
            tree = ET.parse(source_path)
            root = tree.getroot()
        else:
            root = ET.Element("ForceField")
            ET.SubElement(root, "AtomTypes")
            ET.SubElement(root, "Residues")
            
        atom_types_node = root.find("AtomTypes")
        residues_node = root.find("Residues")
        
        for params in new_params_list:
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
                    
        xml_str = minidom.parseString(ET.tostring(root)).toprettyxml(indent="  ")
        clean_xml = "\n".join([line for line in xml_str.split('\n') if line.strip()])
        
        with open(self.update_ff_path, "w") as f:
            f.write(clean_xml)
        print(f"Global force field updated at {self.update_ff_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="High Throughput Workflow for FFBuilder")
    parser.add_argument("-c", "--candidates", default="candidates.txt", help="Path to candidates SMILES file")
    parser.add_argument("-t", "--training", default="training_data.csv", help="Path to training data CSV (smiles, property)")
    parser.add_argument("-n", "--nsuggest", type=int, default=5, help="Number of molecules to suggest")
    parser.add_argument("-r", "--repeats", type=int, default=1, help="Number of repeat units for polymers")
    
    args = parser.parse_args()
    
    workflow = HighThroughputWorkflow(args.candidates, args.training)
    workflow.run_iteration(n_suggest=args.nsuggest, repeats=args.repeats)
