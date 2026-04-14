import xml.etree.ElementTree as ET

def extract_charges(xml_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    charges = {}
    for residue in root.findall(".//Residue"):
        res_name = residue.get("name")
        res_charges = []
        for atom in residue.findall("Atom"):
            c = atom.get("charge")
            if c is not None:
                res_charges.append({
                    "atom": atom.get("name"),
                    "charge": float(c)
                })
        charges[res_name] = res_charges
    return charges

def create_parity():
    old_salt = extract_charges("forcefields/opls_salt.xml")
    new_salt = extract_charges("forcefields/opls_salt_resp.xml")
    old_solv = extract_charges("forcefields/opls_solvent.xml")
    new_solv = extract_charges("forcefields/opls_solvent_resp.xml")
    
    print(f"{'Molecule':<10} | {'Atom':<6} | {'Original':>10} | {'New(CHelpG)':>12} | {'Diff':>10}")
    print("-" * 65)
    
    # Compare Salt Molecules
    for res in ["BF4", "LiA", "NaA", "PF6", "TFS"]:
        o_list = old_salt.get(res, [])
        n_list = new_salt.get(res, [])
        if o_list and n_list:
            for i in range(min(len(o_list), 2)):
                o, n = o_list[i], n_list[i]
                print(f"{res:<10} | {o['atom']:<6} | {o['charge']:10.4f} | {n['charge']:12.4f} | {n['charge']-o['charge']:10.4f}")
    
    print("-" * 65)
    
    # Compare Solvents
    for res in ["ECA", "DMC", "PCA", "DEC"]:
        o_list = old_solv.get(res, [])
        n_list = new_solv.get(res, [])
        if o_list and n_list:
            for i in range(min(len(o_list), 3)):
                o, n = o_list[i], n_list[i]
                print(f"{res:<10} | {o['atom']:<6} | {o['charge']:10.4f} | {n['charge']:12.4f} | {n['charge']-o['charge']:10.4f}")

if __name__ == "__main__":
    create_parity()
