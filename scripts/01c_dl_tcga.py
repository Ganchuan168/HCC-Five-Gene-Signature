import os, time, requests, sys
ROOT = os.environ.get("HCC_ROOT", os.getcwd())
base = os.path.join(ROOT, "data/raw/tcga")
man=[l.rstrip("\n").split("\t") for l in open(f"{base}/file_manifest.tsv")][1:]
os.makedirs(f"{base}/counts",exist_ok=True)
ok=0; skipped=0
for fid,fname,size,*_ in man:
    p=f"{base}/counts/{fname}"; size=int(size)
    if os.path.exists(p) and os.path.getsize(p)==size: skipped+=1; continue
    for attempt in range(6):
        try:
            have=os.path.getsize(p) if os.path.exists(p) else 0
            h={"Range":f"bytes={have}-"} if have else {}
            r=requests.get(f"https://api.gdc.cancer.gov/data/{fid}",headers=h,stream=True,timeout=120)
            mode="ab" if (have and r.status_code==206) else "wb"
            with open(p,mode) as f:
                for chunk in r.iter_content(65536): f.write(chunk)
            if os.path.getsize(p)==size: ok+=1; break
        except Exception as e:
            time.sleep(3)
    sys.stdout.flush()
print("downloaded",ok,"skipped",skipped,"total",len(man))
