import os
import subprocess
import math
import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def create_long_smiles(
    smile,
    repeats=None,
    add_end_Cs=True,
    reaction="[Cu][*:1].[*:2][Au]>>[*:1]-[*:2]",
    product_index=0,
):
    """
    Creates a polymer SMILES string from a monomer SMILES with [Cu] and [Au] as attachment points.
    Mimics HiTPoly's create_long_smiles.
    """
    if "Cu" in smile:
        if repeats is None:
            repeats = 1
        
        if repeats > 1:
            try:
                mol = Chem.MolFromSmiles(smile)
                new_mol = mol
                rxn = AllChem.ReactionFromSmarts(reaction)
                
                for _ in range(repeats - 1):
                    results = rxn.RunReactants((mol, new_mol))
                    if not results or not results[0]:
                        break
                    new_mol = results[product_index][0]
                
                new_smile = Chem.MolToSmiles(new_mol)
            except Exception as e:
                print(f"Error generating polymer: {e}")
                return "None", 0
        else:
            new_smile = smile
            
        if add_end_Cs:
            new_smile = new_smile.replace("[Cu]", "C").replace("[Au]", "C").replace("[Ca]", "")
        else:
            new_smile = new_smile.replace("[Cu]", "").replace("[Au]", "").replace("[Ca]", "")
    else:
        new_smile = smile
        repeats = 1

    try:
        long_smile = Chem.MolToSmiles(Chem.MolFromSmiles(new_smile))
        return long_smile, repeats
    except:
        return "None", 0

def run_ligpargen_local(smiles, name, work_dir=".", charge=0, optimization=0, res_name=None):
    """
    Runs LigParGen locally. 
    Assumes LigParGen is installed and 'LigParGen' command is in path.
    Requires BOSSdir environment variable to be set for CM1A/OPLS-AA generation.
    """
    if 'BOSSdir' not in os.environ:
        print("Warning: $BOSSdir is not set. LigParGen will likely fail for local OPLS-AA generation.")
        # return False # Uncomment if you want to enforce BOSS presence
    
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES: {smiles}")
        return False
    mol = Chem.AddHs(mol)
    try:
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    except:
        print(f"Warning: Molecule embedding failed for {smiles}")
    
    mol_file = os.path.join(work_dir, f"{name}.mol")
    with open(mol_file, "w") as f:
        f.write(Chem.MolToMolBlock(mol))
    
    if not res_name:
        res_name = (name.split('_')[-1] if '_' in name else name)[:3].upper()
    if not res_name: res_name = "MOL"
    
    # LigParGen 2.3 CLI usage: LigParGen -m phenol.mol -r PHN -c 0 -o 0 -l
    command = f"LigParGen -m {name}.mol -o {optimization} -c {charge} -r {res_name} -l"
    print(f"Running command: {command} in {work_dir}")
    
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True, cwd=work_dir)
        if result.returncode != 0:
            print(f"LigParGen failed for {name}: {result.stderr}")
            return False
        
        # Check if output XML exists
        xml_path = os.path.join(work_dir, f"{res_name}.xml")
        if not os.path.exists(xml_path):
            print(f"Expected output XML not found at {xml_path}")
            return False
        return True
    except Exception as e:
        print(f"Error running LigParGen: {e}")
        return False

def extract_parameters_from_lpg_xml(xml_path, name):
    """
    Extracts parameters from LigParGen generated XML and formats them for OpenMM ForceField XML.
    """
    if not os.path.exists(xml_path):
        return None
    
    tree = ET.parse(xml_path)
    root = tree.getroot()
    
    # LigParGen XML usually has:
    # <AtomTypes>, <Residues>, <HarmonicBondForce>, <HarmonicAngleForce>, <PeriodicTorsionForce>, <NonbondedForce>
    
    # We need to rename types to avoid conflicts when merging
    type_map = {}
    
    # 1. Process AtomTypes
    atom_types = []
    at_node = root.find("AtomTypes")
    if at_node is not None:
        for t in at_node.findall("Type"):
            old_name = t.get("name")
            new_name = f"{name}_{old_name}"
            type_map[old_name] = new_name
            t.set("name", new_name)
            t.set("class", new_name)
            atom_types.append(t)
            
    # 2. Process Residues
    residues = []
    res_node = root.find("Residues")
    if res_node is not None:
        for res in res_node.findall("Residue"):
            res.set("name", name.upper())
            for atom in res.findall("Atom"):
                atom.set("type", type_map[atom.get("type")])
            residues.append(res)
            
    # 3. Process Forces
    forces = {}
    for tag in ["HarmonicBondForce", "HarmonicAngleForce", "PeriodicTorsionForce"]:
        f_node = root.find(tag)
        if f_node is not None:
            for item in f_node:
                for attr in ["type1", "type2", "type3", "type4", "class1", "class2", "class3", "class4"]:
                    if item.get(attr) in type_map:
                        item.set(attr, type_map[item.get(attr)])
            forces[tag] = f_node
            
    # 4. NonbondedForce
    nb_node = root.find("NonbondedForce")
    if nb_node is not None:
        for atom in nb_node.findall("Atom"):
            if atom.get("type") in type_map:
                atom.set("type", type_map[atom.get("type")])
        forces["NonbondedForce"] = nb_node
        
    return {
        "AtomTypes": atom_types,
        "Residues": residues,
        "Forces": forces
    }

if __name__ == "__main__":
    # Test polymer generation
    monomer = "[Cu]COC[Au]" # PEO monomer
    poly_smi, reps = create_long_smiles(monomer, repeats=3)
    print(f"Polymer SMILES (3 units): {poly_smi}")
