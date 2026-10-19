"""Mirror FlowRepository FR-FCM-ZZPH file by file over HTTP, with retries and MD5 checks.

FlowRepository's TLS certificate had expired on 2026-10-18, so this uses plain HTTP.
The per-file MD5 prefixes shown on the experiment page are compared after download.
Usage: python3 mirror_zzph.py <data_dir>
"""
import hashlib, json, sys, time, urllib.request
from pathlib import Path

BASE = "http://flowrepository.org"
ATTACHMENTS = {  # page order matched the attachment id order
    "population_names_Levine_32dim.txt": "/experiments/817/attachments/1514/download",
    "population_names_Levine_13dim.txt": "/experiments/817/attachments/1515/download",
    "population_names_Samusik.txt": "/experiments/817/attachments/1516/download",
    "GatingML_Nilsson_rare.xml": "/experiments/817/attachments/1517/download",
    "GatingML_Mosmann_rare.xml": "/experiments/817/attachments/1518/download",
}

def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def fetch(url, dest, tries=4):
    for t in range(1, tries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "mirror_zzph/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
                ctype = r.headers.get("Content-Type", "")
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            if "text/html" in ctype and dest.suffix == ".fcs":
                raise RuntimeError(f"got HTML instead of FCS ({ctype})")
            return ctype
        except Exception as e:
            print(f"  try {t} failed: {e}", flush=True)
            time.sleep(10 * t)
    raise SystemExit(f"gave up on {url}")

def main(data_dir):
    data_dir = Path(data_dir); data_dir.mkdir(parents=True, exist_ok=True)
    files = json.load(open(Path(__file__).with_name("zzph_files.json")))
    for f in files:
        if f["link"] is None:
            f["link"] = ATTACHMENTS[f["name"]]
    report = []
    for f in files:
        dest = data_dir / f["name"]
        if dest.exists() and md5(dest).startswith(f["md5_prefix"] or ""):
            print(f"have {f['name']}", flush=True)
        else:
            print(f"get  {f['name']} ({f['size']})", flush=True)
            fetch(BASE + f["link"], dest)
        sum_ = md5(dest)
        ok = sum_.startswith(f["md5_prefix"]) if f["md5_prefix"] else None
        report.append((f["name"], dest.stat().st_size, sum_, f["md5_prefix"], ok))
        print(f"  md5 {sum_} expected prefix {f['md5_prefix']} match={ok}", flush=True)
    with open(data_dir / "CHECKSUMS.md5", "w") as out:
        for name, size, sum_, prefix, ok in report:
            out.write(f"{sum_}  {name}\n")
    bad = [r for r in report if r[4] is False]
    print("DONE", "all checksums match" if not bad else f"MISMATCH: {[b[0] for b in bad]}")

if __name__ == "__main__":
    main(sys.argv[1])
