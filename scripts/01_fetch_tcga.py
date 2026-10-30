import requests, json, os, csv
API="https://api.gdc.cancer.gov"
filters={"op":"and","content":[
 {"op":"in","content":{"field":"cases.project.project_id","value":["TCGA-LIHC"]}},
 {"op":"in","content":{"field":"data_category","value":["Transcriptome Profiling"]}},
 {"op":"in","content":{"field":"data_type","value":["Gene Expression Quantification"]}},
 {"op":"in","content":{"field":"analysis.workflow_type","value":["STAR - Counts"]}},
 {"op":"in","content":{"field":"access","value":["open"]}}]}
params={"filters":json.dumps(filters),"fields":"file_id,file_name,file_size,cases.submitter_id,cases.samples.sample_type,cases.samples.submitter_id",
        "format":"JSON","size":2000}
r=requests.post(f"{API}/files",json=params,timeout=120).json()
hits=r["data"]["hits"]
print("files:",len(hits))
with open("data/raw/tcga/file_manifest.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["file_id","file_name","size","case","sample","sample_type"])
    for h in hits:
        c=h["cases"][0]; s=c["samples"][0]
        w.writerow([h["file_id"],h["file_name"],h["file_size"],c["submitter_id"],s["submitter_id"],s["sample_type"]])
r2=requests.get(f"{API}/projects/TCGA-LIHC",params={"expand":"summary,diagnoses,demographic","format":"JSON"},timeout=120)
open("data/raw/tcga/clinical.json","w").write(r2.text)
print("clinical saved", len(r2.text))
