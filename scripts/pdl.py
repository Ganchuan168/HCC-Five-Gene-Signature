#!/usr/bin/env python3
"""分段并行下载器：每段独立校验长度，合并后校验 gzip。"""
import os, sys, gzip, time, requests
from concurrent.futures import ThreadPoolExecutor

def head_size(url):
    r = requests.head(url, timeout=60, allow_redirects=True)
    n = r.headers.get("Content-Length")
    if not n:
        r = requests.get(url, headers={"Range": "bytes=0-0"}, timeout=60)
        n = r.headers["Content-Range"].split("/")[-1]
    return int(n)

def fetch_part(url, start, end, path, tries=8):
    if os.path.exists(path) and os.path.getsize(path) == end - start + 1:
        return True
    for _ in range(tries):
        try:
            have = os.path.getsize(path) if os.path.exists(path) else 0
            r = requests.get(url, headers={"Range": f"bytes={start+have}-{end}"}, stream=True, timeout=180)
            mode = "ab" if (have and r.status_code == 206) else "wb"
            if mode == "wb":
                have = 0
            with open(path, mode) as f:
                for c in r.iter_content(262144):
                    f.write(c)
            if os.path.getsize(path) == end - start + 1:
                return True
        except Exception:
            pass
        time.sleep(2)
    return False

def download(url, out, nparts=12):
    total = head_size(url)
    if os.path.exists(out) and os.path.getsize(out) == total and gzip_ok(out):
        print("skip", out); return True
    step = total // nparts + 1
    ranges = [(i*step, min((i+1)*step-1, total-1)) for i in range(nparts)]
    tmp = out + ".parts"; os.makedirs(tmp, exist_ok=True)
    with ThreadPoolExecutor(max_workers=8) as ex:
        res = list(ex.map(lambda r: fetch_part(url, r[0], r[1], f"{tmp}/p{r[0]}"), ranges))
    if not all(res):
        print("part fail", out); return False
    with open(out, "wb") as w:
        for r in ranges:
            with open(f"{tmp}/p{r[0]}", "rb") as f:
                w.write(f.read())
    for r in ranges: os.remove(f"{tmp}/p{r[0]}")
    os.rmdir(tmp)
    ok = gzip_ok(out)
    print(("OK " if ok else "CORRUPT ") + out, os.path.getsize(out))
    return ok

def gzip_ok(p):
    try:
        with gzip.open(p, "rb") as f:
            while f.read(1 << 24):
                pass
        return True
    except Exception:
        return False

if __name__ == "__main__":
    out, url = sys.argv[1], sys.argv[2]
    for attempt in range(5):
        if download(url, out, int(sys.argv[3]) if len(sys.argv) > 3 else 12):
            break
