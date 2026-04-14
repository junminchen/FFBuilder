#!/usr/bin/env python
import os
import sys
import time
import re
import shutil
import requests
import urllib3
import numpy as np
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom
from rdkit import Chem
from rdkit.Chem import AllChem
import MDAnalysis as mda

# Try loading user-defined molecule list
try:
    import sys as _sys
    _sys.path.insert(0, os.path.dirname(__file__))
    import utils.mol_list as mol_list
    if hasattr(mol_list, 'smiles_batch'):
        molecules = mol_list.smiles_batch
    elif hasattr(mol_list, 'molecules'):
        molecules = mol_list.molecules
    else:
        molecules = {"FAN": "FCC#N"}
except ImportError:
    molecules = {"FAN": "FCC#N"}

# RESP support (optional)
try:
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from native.resp.RESP_Workflow import process_molecule as calc_resp_charges
    from utils.FFutils import get_molecule_charge, read_pdb_skeleton
    RESP_AVAILABLE = True
except ImportError:
    RESP_AVAILABLE = False

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

def print_info(stage, message):
    print(f"[{time.strftime('%H:%M:%S')}] [{stage}] {message}")

def padding(n):
    return str(n).zfill(2)

class OPLSWorkflow:
    def __init__(self, mol_dict):
        self.mol_dict = mol_dict
        # Paths are relative to project root (parent of web/)
        project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.mol_root = os.path.join(project_root, "data", "molecules")
        self.base_url = "https://traken.chem.yale.edu"
        self.session = requests.Session()
        
        if not os.path.exists(self.mol_root): os.makedirs(self.mol_root)

    def _get_symmetry_dict(self, pdb_path):
        mol = Chem.MolFromPDBFile(pdb_path, sanitize=False)
        if mol is None: return None
        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
        unique_ranks = sorted(list(set(ranks)))
        result_dict = {}
        for atom_idx, rank in enumerate(ranks):
            new_id = unique_ranks.index(rank) + 1
            key = str(new_id)
            if key not in result_dict: result_dict[key] = []
            result_dict[key].append(atom_idx)
        return result_dict

    def stage_1_download(self, name, smiles):
        print_info(name, "Stage 1: Downloading...")
        payload = {'smiData': (None, smiles), 'molpdbfile': ('', b''), 'checkopt': (None, ' 0 '), 'chargetype': (None, 'cm1abcc'), 'dropcharge': (None, ' 0 ')}
        try:
            res = self.session.post(f"{self.base_url}/cgi-bin/results_lpg.py", files=payload, timeout=300, verify=False)
            # traken.yale.edu: go (submit) and fileout (hidden) are SEPARATE <input> tags
            # → pair them by appearance order (zip), not by form grouping
            go_vals = re.findall(r'<input type="submit"[^>]+name="go"[^>]+value="([^"]+)"', res.text)
            fileout_vals = re.findall(r'<input type="hidden"[^>]+name="fileout"[^>]+value="([^"]+)"', res.text)
            if not fileout_vals or "XML" not in go_vals:
                print_info(name, f"Download failed: no XML in go_vals {go_vals}")
                return False
            file_map = dict(zip(go_vals, fileout_vals))
            mol_dir = os.path.join(self.mol_root, name)
            if not os.path.exists(mol_dir): os.makedirs(mol_dir)
            # Save raw XML to temp, write PDB directly (will be fixed in stage 2)
            for label in ["XML", "PDB"]:
                r = self.session.post(f"{self.base_url}/cgi-bin/download_lpg.py", data={'go': label, 'fileout': file_map[label]}, verify=False)
                if label == "XML":
                    with open(os.path.join(mol_dir, f"{name}.xml"), 'wb') as f: f.write(r.content)
                else:
                    with open(os.path.join(mol_dir, f"{name}.pdb"), 'wb') as f: f.write(r.content)
            return True
        except Exception as e:
            print_info(name, f"Error: {e}"); return False

    def stage_2_fix_pdb(self, name):
        print_info(name, "Stage 2: Fixing PDB...")
        mol_dir = os.path.join(self.mol_root, name)
        pdb_file = os.path.join(mol_dir, f"{name}.pdb")
        try:
            u = mda.Universe(pdb_file)
            if not hasattr(u, 'elements') or u.atoms.elements[0] == '':
                u.add_TopologyAttr('elements', [mda.topology.guessers.guess_atom_element(a.name) for a in u.atoms])
            resid = (name + 'A')[:3].upper()
            for i, atom in enumerate(u.atoms): atom.name = atom.element + padding(i)
            u.atoms.residues.resnames = resid
            u.atoms.write(pdb_file)
            return True
        except Exception: return False

    def stage_3_consolidate_ff(self, name):
        print_info(name, "Stage 3: Consolidating Symmetry...")
        mol_dir = os.path.join(self.mol_root, name)
        xml_file = os.path.join(mol_dir, f"{name}.xml")
        target_xml = os.path.join(mol_dir, f"{name}.xml")

        dic_atypes = self._get_symmetry_dict(os.path.join(mol_dir, f"{name}.pdb"))
        if not dic_atypes: return False
        
        # 使用分子名作为前缀，确保全局唯一性
        new_dic = {f"{name}_{int(k)}": v for k, v in dic_atypes.items()}
        
        try:
            root = ET.parse(xml_file).getroot()
            map_types = {}; map_classes = {}
            new_root = ET.Element("ForceField")
            new_types_node = ET.SubElement(new_root, "AtomTypes")
            
            for node in next(c for c in list(root) if c.tag == 'AtomTypes'):
                old_idx = int(node.attrib['name'][5:]) - 800
                new_type_name = next(k for k, v in new_dic.items() if old_idx in v)
                map_types[node.attrib['name']] = new_type_name
                map_classes[node.attrib['class']] = new_type_name
                if not any(t.attrib['name'] == new_type_name for t in new_types_node):
                    ET.SubElement(new_types_node, "Type", attrib={'name': new_type_name, 'class': new_type_name, 'element': node.attrib['element'], 'mass': node.attrib['mass']})

            new_res_node = ET.SubElement(new_root, "Residues")
            res_val = ET.SubElement(new_res_node, "Residue", attrib={'name': name.upper()})
            for i, atom_node in enumerate(next(c for c in list(root) if c.tag == 'Residues')[0].findall("Atom")):
                ET.SubElement(res_val, "Atom", attrib={'name': atom_node.attrib['name'][0]+padding(i), 'type': map_types[atom_node.attrib['type']]})
            for b in next(c for c in list(root) if c.tag == 'Residues')[0].findall("Bond"): ET.SubElement(res_val, "Bond", attrib=b.attrib)

            for tag in ["HarmonicBondForce", "HarmonicAngleForce", "PeriodicTorsionForce"]:
                old_f = next(c for c in list(root) if c.tag == tag)
                new_f_node = ET.SubElement(new_root, tag)
                seen = set()
                for item in old_f:
                    attr = item.attrib.copy()
                    class_keys = [k for k in attr if k.startswith('class')]
                    classes = [map_classes[attr[ck]] for ck in class_keys]
                    for ck, cv in zip(class_keys, classes): attr[ck] = cv
                    if tuple(classes) not in seen and tuple(classes)[::-1] not in seen:
                        ET.SubElement(new_f_node, item.tag, attrib=attr); seen.add(tuple(classes))

            nb_node = next(c for c in list(root) if c.tag == 'NonbondedForce')
            new_nb = ET.SubElement(new_root, "NonbondedForce", attrib=nb_node.attrib)
            nb_params = {}
            for atom in nb_node:
                at = map_types[atom.attrib['type']]
                if at not in nb_params: nb_params[at] = []
                nb_params[at].append([float(atom.attrib['charge']), float(atom.attrib['sigma']), float(atom.attrib['epsilon'])])
            for at, vals in nb_params.items():
                avg = np.average(np.array(vals), axis=0)
                ET.SubElement(new_nb, "Atom", attrib={'type': at, 'charge': f"{avg[0]:.6f}", 'sigma': f"{avg[1]:.6f}", 'epsilon': f"{avg[2]:.6f}"})

            xml_str = minidom.parseString(ET.tostring(new_root)).toprettyxml(indent="  ")
            with open(target_xml, 'w') as f: f.write("\n".join(l for l in xml_str.split("\n") if l.strip()))
            return True
        except Exception as e:
            print_info(name, f"FF Error: {e}"); return False

    def stage_3b_resp_charges(self, name, smiles):
        """Replace CM1A charges in <name>.xml with RESP charges."""
        if not RESP_AVAILABLE:
            print_info(name, "RESP not available, skipping.")
            return False
        print_info(name, "Stage 3b: Computing RESP charges...")
        mol_dir = os.path.join(self.mol_root, name)
        xml_file = os.path.join(mol_dir, f"{name}.xml")
        if not os.path.exists(xml_file):
            print_info(name, f"XML not found at {xml_file}, skipping RESP.")
            return False
        try:
            total_charge = get_molecule_charge(name, smiles)
            symbols, coords, resp_charges, _ = calc_resp_charges(name, smiles, total_charge)
            tree = ET.parse(xml_file)
            root = tree.getroot()
            res = root.find(".//Residue")
            if res is None:
                print_info(name, "No Residue found in XML, skipping.")
                return False
            atoms = res.findall("Atom")
            if len(atoms) != len(resp_charges):
                print_info(name, f"Atom count mismatch ({len(atoms)} XML vs {len(resp_charges)} RESP), skipping.")
                return False
            for i, atom in enumerate(atoms):
                atom.set("charge", f"{resp_charges[i]:.6f}")
            xml_str = minidom.parseString(ET.tostring(root)).toprettyxml(indent="  ")
            with open(xml_file, 'w') as f:
                f.write("\n".join(l for l in xml_str.split("\n") if l.strip()))
            print_info(name, f"RESP charges applied (total={sum(resp_charges):.4f}).")
            return True
        except Exception as e:
            print_info(name, f"RESP Error: {e}"); return False

    def stage_4_incremental_merge(self):
        print("\n" + "="*50)
        print_info("GLOBAL", f"Stage 4: Incrementally merging into {self.update_xml_path}...")
        
        # 1. 加载现有的底座 XML (如果不存在则初始化空字典)
        merged_dict = {}
        existing_residues = set()
        
        if os.path.exists(self.base_xml_path):
            print_info("GLOBAL", f"Loading base: {self.base_xml_path}")
            base_root = ET.parse(self.base_xml_path).getroot()
            for elem in base_root:
                tag = elem.tag
                if tag not in merged_dict: merged_dict[tag] = []
                # 记录已存在的 Residue
                if tag == "Residues":
                    for res in elem.findall("Residue"): existing_residues.add(res.get("name"))
                # 将旧元素存入字典
                merged_dict[tag].append(elem)
        
        # 2. 将 mol_dict 中的新分子追加进去
        added_count = 0
        for name in self.mol_dict.keys():
            residue_name = name.upper()
            if residue_name in existing_residues:
                print_info(name, "Already exists in base XML, skipping merge.")
                continue
            
            xml_file = os.path.join(self.mol_root, name, f"{name}.xml")
            if os.path.exists(xml_file):
                print_info(name, "Appending new parameters to categories...")
                mol_root = ET.parse(xml_file).getroot()
                added_count += 1
                for elem in mol_root:
                    if elem.tag not in merged_dict: merged_dict[elem.tag] = []
                    merged_dict[elem.tag].append(elem)

        if added_count == 0 and not existing_residues:
            print_info("GLOBAL", "No parameters found to merge.")
            return

        # 3. 构造最终的力场 DOM
        final_root = ET.Element('ForceField')
        for tag, elems in merged_dict.items():
            # 为每一个类别创建一个统一的父标签（如 <AtomTypes>, <Residues>）
            container = ET.SubElement(final_root, tag)
            # 继承第一个找到的属性（尤其是 NonbondedForce 的比例因子）
            if tag == "NonbondedForce": container.attrib = elems[0].attrib
            for parent_node in elems:
                for child in parent_node: container.append(child)

        # 4. 保存为 _update 版本
        xml_str = minidom.parseString(ET.tostring(final_root)).toprettyxml(indent="    ")
        with open(self.update_xml_path, 'w') as f:
            f.write("\n".join(l for l in xml_str.split("\n") if l.strip()))
        
        print_info("GLOBAL", f"Incremental Merge Complete! Updated file: {self.update_xml_path}")
        print("="*50 + "\n")

    def run(self, use_resp=False):
        for name, smiles in self.mol_dict.items():
            xml_cache = os.path.join(self.mol_root, name, f"{name}.xml")
            pdb_cache = os.path.join(self.mol_root, name, f"{name}.pdb")
            if os.path.exists(xml_cache) and os.path.exists(pdb_cache):
                print_info(name, "Cache hit, skipping.")
                continue
            if self.stage_1_download(name, smiles):
                self.stage_2_fix_pdb(name)
                self.stage_3_consolidate_ff(name)
                if use_resp:
                    self.stage_3b_resp_charges(name, smiles)
                time.sleep(5)
        print_info("DONE", f"XML files → {self.mol_root}/<name>/<name>.xml")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="OPLS ForceField Workflow via traken.yale.edu")
    parser.add_argument("-r", "--resp", action="store_true", help="Compute RESP charges and replace CM1A charges in XML")
    args = parser.parse_args()
    workflow = OPLSWorkflow(molecules)
    workflow.run(use_resp=args.resp)
