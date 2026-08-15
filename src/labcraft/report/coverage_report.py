import json
import os
from typing import List, Dict
from jinja2 import Environment, FileSystemLoader

from labcraft.lamp.coverage import StrainVerdict, SiteVerdict

def render_coverage_report(
    verdicts: List[StrainVerdict],
    output_path: str,
    panel_name: str,
    dg_threshold: float,
    max_mismatches_count: int
):
    # Prepare data
    total_strains = len(verdicts)
    if total_strains == 0:
        with open(output_path, 'w') as f:
            f.write("<html><body>No strains analyzed</body></html>")
        return
        
    covered_thermo = sum(1 for v in verdicts if v.is_amplifiable_thermo)
    covered_count = sum(1 for v in verdicts if v.is_amplifiable_count)
    
    divergences = [v for v in verdicts if v.is_amplifiable_thermo != v.is_amplifiable_count]
    
    # Limiting primers thermo
    limiting_counts = {}
    for v in verdicts:
        if not v.is_amplifiable_thermo and v.limiting_primer_thermo:
            limiting_counts[v.limiting_primer_thermo] = limiting_counts.get(v.limiting_primer_thermo, 0) + 1
            
    sorted_limiting = sorted(limiting_counts.items(), key=lambda x: x[1], reverse=True)
    worst_primer = sorted_limiting[0][0] if sorted_limiting else None
    
    # Collect all unique primers across all evaluations for the matrix header
    all_primers = set()
    for v in verdicts:
        all_primers.update(v.evaluations.keys())
    # Sort them generically (F3, B3, FIP, BIP, LF, LB)
    role_order = {"F3":0, "B3":1, "FIP":2, "BIP":3, "LF":4, "LB":5}
    def sort_key(p_name):
        base = p_name.split("_")[0]
        return role_order.get(base, 99)
    sorted_primers = sorted(list(all_primers), key=sort_key)
    
    # Prepare JSON/CSV output
    json_path = output_path.replace(".html", ".json")
    json_data = {
        "panel_name": panel_name,
        "total_strains": total_strains,
        "covered_thermo": covered_thermo,
        "covered_count": covered_count,
        "dg_threshold": dg_threshold,
        "worst_primer": worst_primer,
        "strains": []
    }
    for v in verdicts:
        s_data = {
            "strain_id": v.strain_id,
            "thermo": v.is_amplifiable_thermo,
            "count": v.is_amplifiable_count,
            "evals": {}
        }
        for pname, ev in v.evaluations.items():
            s_data["evals"][pname] = {
                "verdict": ev.verdict.value,
                "dg": ev.dg_hyb,
                "mismatches": ev.n_mismatches_count
            }
        json_data["strains"].append(s_data)
        
    with open(json_path, 'w') as f:
        json.dump(json_data, f, indent=2)
        
    # Render HTML
    env = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), 'templates')))
    template = env.get_template('coverage.html.j2')
    
    html = template.render(
        panel_name=panel_name,
        total_strains=total_strains,
        covered_thermo=covered_thermo,
        covered_count=covered_count,
        thermo_rate=round(covered_thermo / total_strains * 100, 1),
        count_rate=round(covered_count / total_strains * 100, 1),
        divergences=divergences,
        worst_primer=worst_primer,
        sorted_primers=sorted_primers,
        verdicts=verdicts,
        dg_threshold=dg_threshold,
        max_mismatches_count=max_mismatches_count,
        SiteVerdict=SiteVerdict
    )
    
    with open(output_path, 'w') as f:
        f.write(html)
