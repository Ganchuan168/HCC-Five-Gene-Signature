import requests, json
API="https://api.gdc.cancer.gov/cases"
fields=["submitter_id","demographic.sex_at_birth","demographic.age_at_index","demographic.vital_status",
        "demographic.days_to_death","diagnoses.age_at_diagnosis","diagnoses.ajcc_pathologic_stage",
        "diagnoses.ajcc_pathologic_t","diagnoses.ajcc_pathologic_n","diagnoses.ajcc_pathologic_m",
        "diagnoses.tumor_grade","diagnoses.primary_diagnosis","diagnoses.days_to_last_follow_up",
        "diagnoses.vascular_invasion","diagnoses.child_pugh_classification","diagnoses.ishak_fibrosis_score"]
p={"filters":json.dumps({"op":"in","content":{"field":"project.project_id","value":["TCGA-LIHC"]}}),
   "fields":",".join(fields),"format":"JSON","size":1000}
r=requests.get(API,params=p,timeout=120).json()
hits=r["data"]["hits"]; print("cases",len(hits))
json.dump(hits,open("data/raw/tcga/cases.json","w"))
