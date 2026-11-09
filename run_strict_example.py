import csv
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parent
with (root/"data/replay-edges.csv").open(newline="") as stream:
    edges = [(int(row["u"]),int(row["v"])) for row in csv.DictReader(stream)]
adjacency = [[] for _ in range(max(max(e) for e in edges)+1)]
for u,v in edges:
    adjacency[u].append(v)
    adjacency[v].append(u)
offsets = [0]
for row in adjacency:
    offsets.append(offsets[-1]+len(row))
with tempfile.TemporaryDirectory() as tmp:
    work = Path(tmp)
    graph = work/"replay.bin"
    graph.write_bytes(struct.pack("<QQ",len(adjacency),len(edges))
                      + b"".join(struct.pack("<Q",o) for o in offsets)
                      + b"".join(struct.pack("<I",v) for row in adjacency for v in sorted(row)))
    for arm in ("maintained","warm","cold"):
        target = work/arm
        target.mkdir()
        env = {k:v for k,v in os.environ.items() if not k.startswith(("KCDS_","OMP_","TRACE_"))}
        env.update(OMP_NUM_THREADS="1",TRACE_ARM=arm,TRACE_DIR=str(target),TRACE_PYTHON=sys.executable,
                   TRACE_CHECKER=str(root/"verify_checkpoint.py"),TRACE_NATIVE=str(root/"build/verify_k4.exe"),
                   KCDS_NOCHECK="1")
        args = [root/"build/trace.exe",graph,"4",".01","2","4","replay","1","10",".15","hybrid","32"]
        result = subprocess.run(args,check=True,capture_output=True,text=True,env=env,timeout=180)
        reports = [json.loads(s) for s in result.stdout.splitlines() if s.startswith("{")]
        if len(reports)!=5 or not all(r["rigorous"] and r["epsilon"]=="1/100" for r in reports):
            raise ValueError(f"Missing verified checkpoint for {arm}")
        print(f"{arm}: {len(reports)} verified checkpoints")
