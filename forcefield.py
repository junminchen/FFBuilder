#!/usr/bin/env python
import os
import json
import xml.etree.ElementTree as ET
from typing import Optional, Dict, Any, List

class ForceFieldManager:
    def __init__(self, ff_dir: str = "forcefields"):
        self.ff_dir = ff_dir
        if not os.path.exists(ff_dir):
            os.makedirs(ff_dir)
    
    def get(self, name: str, format: str = "xml") -> Optional[Any]:
        if format == "xml":
            return self._get_xml(name)
        elif format == "dict":
            return self._get_as_dict(name)
        elif format == "json":
            d = self._get_as_dict(name)
            return json.dumps(d, indent=2) if d else None
        return None
    
    def _get_xml(self, name: str) -> Optional[str]:
        for filename in [f"{name}.xml", f"opls_{name}.xml", f"opls_solvent_{name}.xml"]:
            path = os.path.join(self.ff_dir, filename)
            if os.path.exists(path):
                with open(path, 'r') as f:
                    return f.read()
        update_path = os.path.join(self.ff_dir, "opls_solvent_update.xml")
        if os.path.exists(update_path):
            return self._extract_molecule_from_xml(update_path, name.upper())
        return None
    
    def _get_as_dict(self, name: str) -> Optional[Dict[str, Any]]:
        xml_str = self._get_xml(name)
        if not xml_str:
            return None
        return self._xml_to_dict(xml_str)
    
    def _get_as_json(self, name: str) -> Optional[str]:
        d = self._get_as_dict(name)
        return json.dumps(d, indent=2) if d else None
    
    def _extract_molecule_from_xml(self, xml_path: str, mol_name: str) -> Optional[str]:
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
            residues = root.find(".//Residues")
            if residues is None:
                return None
            target_res = None
            for res in residues.findall("Residue"):
                if res.get("name") == mol_name:
                    target_res = res
                    break
            if target_res is None:
                return None
            new_root = ET.Element("ForceField")
            atom_types = root.find("AtomTypes")
            if atom_types is not None:
                new_root.append(ET.fromstring(ET.tostring(atom_types)))
            res_node = ET.SubElement(new_root, "Residues")
            res_node.append(ET.fromstring(ET.tostring(target_res)))
            for tag in ["HarmonicBondForce", "HarmonicAngleForce", "PeriodicTorsionForce", "NonbondedForce"]:
                parent = root.find(f".//{tag}")
                if parent is not None:
                    new_root.append(ET.fromstring(ET.tostring(parent)))
            return ET.tostring(new_root, encoding='unicode')
        except Exception as e:
            print(f"Error extracting {mol_name}: {e}")
            pass
        return None
    
    def _xml_to_dict(self, xml_str: str) -> Dict[str, Any]:
        def parse_element(elem):
            result = dict(elem.attrib)
            children = list(elem)
            if children:
                result[elem.tag] = [parse_element(c) for c in children]
            return result
        root = ET.fromstring(xml_str)
        return parse_element(root)
    
    def list_available(self) -> List[str]:
        available = []
        for f in os.listdir(self.ff_dir):
            if f.endswith(".xml"):
                try:
                    tree = ET.parse(os.path.join(self.ff_dir, f))
                    root = tree.getroot()
                    for res in root.findall(".//Residue"):
                        available.append(res.get("name"))
                except Exception:
                    pass
        return sorted(set(available))
    
    def get_by_molecule(self, mol_name: str) -> Optional[str]:
        return self._get_xml(mol_name)

def load_forcefield(name: str = "opls_solvent_update", ff_dir: str = "forcefields") -> Optional[str]:
    manager = ForceFieldManager(ff_dir)
    return manager.get(name)

if __name__ == "__main__":
    mgr = ForceFieldManager()
    print("Available molecules:", mgr.list_available())
    print("\nExample - get FAN as XML (first 200 chars):")
    xml = mgr.get("FAN", format="xml")
    print(xml[:200] if xml else "Not found")
