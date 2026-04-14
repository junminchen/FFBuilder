import xml.etree.ElementTree as ET
import math

def extract_all_charges(old_xml, new_xml):
    def get_charges(path):
        tree = ET.parse(path)
        data = {}
        for res in tree.findall(".//Residue"):
            res_name = res.get("name")
            for atom in res.findall("Atom"):
                c = atom.get("charge")
                if c is not None:
                    data[(res_name, atom.get("name"))] = float(c)
        return data

    old_data = get_charges(old_xml)
    new_data = get_charges(new_xml)
    
    pairs = []
    for key in old_data:
        if key in new_data:
            pairs.append((old_data[key], new_data[key], key[0] + ":" + key[1]))
    return pairs

def generate_svg(points, filename="charge_parity.svg"):
    # SVG Settings
    width, height = 600, 600
    margin = 60
    
    # Scale calculation
    all_vals = [p[0] for p in points] + [p[1] for p in points]
    min_v, max_v = min(all_vals) - 0.1, max(all_vals) + 0.1
    
    def transform(v, size, margin, v_min, v_max):
        return margin + (v - v_min) / (v_max - v_min) * (size - 2 * margin)

    svg = [f'<svg width="{width}" height="{height}" xmlns="http://www.w3.org/2000/svg" style="background:white">']
    
    # Axes and Grid
    svg.append(f'<line x1="{margin}" y1="{height-margin}" x2="{width-margin}" y2="{height-margin}" stroke="black" />') # X
    svg.append(f'<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height-margin}" stroke="black" />') # Y
    
    # Diagonal Line (y=x)
    x1 = transform(min_v, width, margin, min_v, max_v)
    y1 = height - transform(min_v, height, margin, min_v, max_v)
    x2 = transform(max_v, width, margin, min_v, max_v)
    y2 = height - transform(max_v, height, margin, min_v, max_v)
    svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#ccc" stroke-dasharray="5,5" />')
    
    # Labels
    svg.append(f'<text x="{width/2}" y="{height-10}" text-anchor="middle" font-family="Arial" font-size="14">Original Charges (OPLS)</text>')
    svg.append(f'<text x="15" y="{height/2}" transform="rotate(-90,15,{height/2})" text-anchor="middle" font-family="Arial" font-size="14">New Charges (CHelpG)</text>')
    svg.append(f'<text x="{width/2}" y="30" text-anchor="middle" font-family="Arial" font-size="18" font-weight="bold">Charge Parity Plot</text>')

    # Points
    for old, new, label in points:
        px = transform(old, width, margin, min_v, max_v)
        py = height - transform(new, height, margin, min_v, max_v)
        
        # Color based on value (positive/negative)
        color = "red" if old > 0 else "blue"
        svg.append(f'<circle cx="{px}" cy="{py}" r="4" fill="{color}" fill-opacity="0.6">')
        svg.append(f'<title>{label}: {old:.3f} -> {new:.3f}</title>')
        svg.append('</circle>')
        
        # Label large changes
        if abs(new - old) > 0.15:
            svg.append(f'<text x="{px+5}" y="{py-5}" font-family="Arial" font-size="10" fill="black">{label}</text>')

    svg.append('</svg>')
    
    with open(filename, "w") as f:
        f.write("\n".join(svg))
    print(f"Parity plot saved to {filename}")

if __name__ == "__main__":
    salt_points = extract_all_charges("forcefields/opls_salt.xml", "forcefields/opls_salt_resp.xml")
    solv_points = extract_all_charges("forcefields/opls_solvent.xml", "forcefields/opls_solvent_resp.xml")
    generate_svg(salt_points + solv_points)
