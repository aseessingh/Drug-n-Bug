import os
import csv
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from kegg_pull.rest import KEGGrest

# List of input files
file_paths = [
    "drug_non_drug_base_KEGG_KO_stats.csv",
    "drug_non_drug_0_base_KEGG_KO_stats.csv",
    "drug_non_drug_base_KEGG_pathway_stats.csv",
    "drug_non_drug_0_base_KEGG_pathway_stats.csv",
    "drug_non_drug_base_KEGG_module_stats.csv",
    "drug_non_drug_0_base_KEGG_module_stats.csv"
]

# Collect all DisplayNames
all_display_names = []
for file in file_paths:
    if os.path.exists(file):
        with open(file, newline='') as f:
            reader = csv.DictReader(f)
            if 'DisplayName' in reader.fieldnames:
                all_display_names.extend(row['DisplayName'] for row in reader)
            else:
                print(f"Warning: 'DisplayName' column not found in {file}")
    else:
        print(f"Warning: File not found: {file}")

# Find duplicates
duplicates = [item for item, count in Counter(all_display_names).items() if count > 1]

# KEGG functions
kegg = KEGGrest()

def get_kegg_name(entry_id, prefix):
    try:
        response = kegg.get(entry_ids=f"{prefix}:{entry_id}")
        if response.status.successful:
            for line in response.text_body.splitlines():
                if line.startswith("NAME"):
                    return line.replace("NAME", "").strip()
        return "not found"
    except Exception:
        return "not found"

def get_kegg_links(target, source_id, prefix):
    try:
        response = kegg.link(target_database=target, source_database=f"{prefix}:{source_id}")
        if response.status.successful:
            ids = [line.split("\t")[1].split(":")[1] for line in response.text_body.strip().splitlines()]
            return sorted(set(ids))
        return []
    except Exception:
        return []

def annotate_entry(entry_id):
    if entry_id.startswith("K"):
        entry_type = "KO"
        name = get_kegg_name(entry_id, "ko")
        modules = get_kegg_links("module", entry_id, "ko")
        pathways = get_kegg_links("pathway", entry_id, "ko")
        module_names = [get_kegg_name(mid, "md") for mid in modules]
        pathway_names = [get_kegg_name(pid, "path") for pid in pathways]
    elif entry_id.startswith("M"):
        entry_type = "Module"
        name = get_kegg_name(entry_id, "md")
        modules = [entry_id]
        module_names = [name]
        pathways = get_kegg_links("pathway", entry_id, "md")
        pathway_names = [get_kegg_name(pid, "path") for pid in pathways]
    elif entry_id.startswith("map") or entry_id.startswith("path"):
        entry_type = "Pathway"
        name = get_kegg_name(entry_id, "path")
        modules = []
        module_names = []
        pathways = [entry_id]
        pathway_names = [name]
    else:
        entry_type = "Unknown"
        name = "unknown"
        modules = []
        module_names = []
        pathways = []
        pathway_names = []

    return {
        "Input_ID": entry_id,
        "Type": entry_type,
        "Name": name,
        "Modules": "; ".join(modules),
        "Module_Names": "; ".join(module_names),
        "Pathways": "; ".join(pathways),
        "Pathway_Names": "; ".join(pathway_names)
    }

# Parallel annotation
with ThreadPoolExecutor() as executor:
    results = list(executor.map(annotate_entry, duplicates))

# Write output
with open("kegg_hierarchy_output.csv", "w", newline='') as f:
    fieldnames = ["Input_ID", "Type", "Name", "Modules", "Module_Names", "Pathways", "Pathway_Names"]
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

print("KEGG annotation completed and saved to kegg_hierarchy_output.csv.")
