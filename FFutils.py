import os
import numpy as np
import xml.etree.ElementTree as ET

def get_molecule_charge(name, smiles):
    """Automatically determine the total charge of a molecule."""
    charge = 0
    if "+" in smiles: charge = 1
    elif "-" in smiles: charge = -1
    
    # Special cases for salts
    name_upper = name.upper()
    if "LI" in name_upper or "NA" in name_upper:
        charge = 1
    elif any(s in name_upper for s in ["PF6", "BF4", "FSI", "TFSI", "DFP", "BOB", "DFOB"]):
        charge = -1
    return charge

def read_pdb_skeleton(pdb_path):
    """Simple PDB parser for symbols and coords."""
    symbols = []
    coords = []
    if not os.path.exists(pdb_path):
        return None, None
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                sym = line[76:78].strip()
                if not sym:
                    # Fallback to atom name
                    name_part = line[12:16].strip()
                    sym = ''.join([c for c in name_part if c.isalpha()])[:1]
                symbols.append(sym)
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return symbols, np.array(coords)

def generate_ff_xml(name, symbols, coords, charges, rd_mol, out_path):
    """Generate a standard OpenMM XML forcefield file."""
    root = ET.Element("ForceField")
    atom_types = ET.SubElement(root, "AtomTypes")
    residues = ET.SubElement(root, "Residues")
    res = ET.SubElement(residues, "Residue", name=name)
    
    for i, (s, q) in enumerate(zip(symbols, charges)):
        t_name = f"opls_{s}_{i}" 
        ET.SubElement(atom_types, "Type", name=t_name, class_="XX", element=s, mass="1.0")
        ET.SubElement(res, "Atom", name=f"{s}{i+1}", type=t_name, charge=str(round(q, 4)))
    
    # Guess connectivity if rd_mol is not provided
    if rd_mol:
        for bond in rd_mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            ET.SubElement(res, "Bond", atom1=f"{symbols[a1]}{a1+1}", atom2=f"{symbols[a2]}{a2+1}")
    else:
        natm = len(symbols)
        for i in range(natm):
            for j in range(i + 1, natm):
                dist = np.linalg.norm(coords[i] - coords[j])
                if dist < 1.65:
                    if symbols[i] == 'H' and symbols[j] == 'H': continue
                    ET.SubElement(res, "Bond", atom1=f"{symbols[i]}{i+1}", atom2=f"{symbols[j]}{j+1}")
        
    tree = ET.ElementTree(root)
    ET.indent(tree, space="  ")
    tree.write(out_path, encoding="utf-8", xml_declaration=True)
    print(f"Forcefield XML saved: {out_path}")

def ensure_mol_dir(name):
    """Ensure a directory for the molecule exists in molecules/."""
    path = f"molecules/{name}"
    os.makedirs(path, exist_ok=True)
    return path
