import json
import pandas as pd
with open('m5.go.bp.v2024.1.Mm.json', 'r') as f:
    json_data = json.load(f)
GOTERM = []
GOID = []
term = list(json_data.keys())
for key, val in json_data.items():
    GOTERM.append(key)
    GOID.append(val["exactSource"])
    
term_id_mapping = pd.DataFrame({"GOTERM" : GOTERM, "GOID" : GOID})

term_id_mapping.to_csv("20241213_m5.go.bp.v2024.1.Mm.term.id_mapping_hj.txt", index = False, mode = 'w', header = True, sep = "\t")
