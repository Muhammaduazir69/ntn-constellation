#!/usr/bin/env python3
"""ns-3 campaign: route recomputation policies over Walker shells (resumable, one jsonl line per run)."""
import itertools, json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.expanduser("~/6g_ntn_ns3/ns-3-dev/build/contrib/ntn-constellation/examples/"
                         "ns3.43-ntn-constellation-route-lifetime-optimized")
OUT = os.path.join(HERE, "out", "ns3_routes.jsonl")
SHELLS = {  # filings: Starlink Gen1 shell 1, Kuiper shell 1, OneWeb phase 1 (phasing F assumed 1 where unpublished)
    "starlink": dict(planes=72, perPlane=22, phasing=39, inc=53.0, altKm=550, star=0),
    "kuiper": dict(planes=34, perPlane=34, phasing=1, inc=51.9, altKm=630, star=0),
    "oneweb": dict(planes=12, perPlane=49, phasing=1, inc=87.9, altKm=1200, star=1),
}
PAIRS = {"lon-nyc": (51.51, -0.13, 40.71, -74.01), "fra-sin": (50.11, 8.68, 1.35, 103.82),
         "syd-lax": (-33.87, 151.21, 34.05, -118.24)}
POLICIES = [dict(policy="snap", period=0.1), dict(policy="snap", period=1.0),
            dict(policy="ondemand", detect=0.05, setup=0.0), dict(policy="ondemand", detect=0.05, setup=1.0),
            dict(policy="asr", window=30.0),
            dict(policy="life", eps=0.01), dict(policy="life", eps=0.05), dict(policy="life", eps=0.10)]
STARTS = [0, 1500, 3000, 4500]
DUR = 1205.0


def key(sh, pr, po, st):
    return json.dumps(dict(shell=sh, pair=pr, start=st, **po), sort_keys=True)


def run(job):
    sh, pr, po, st = job
    a = dict(SHELLS[sh]); la, loa, lb, lob = PAIRS[pr]
    a.update(latA=la, lonA=loa, latB=lb, lonB=lob, duration=DUR, tick=0.1, guard=0.5, startOffset=st, **po)
    cmd = ["nice", "-n", "5", BIN] + [f"--{k}={v}" for k, v in a.items()]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=7200,
                           cwd=os.path.expanduser("~/6g_ntn_ns3/ns-3-dev"))
        res = json.loads(r.stdout.strip().splitlines()[-1])
    except Exception as e:  # keep the campaign going; the failure is recorded
        res = dict(error=str(e)[:300])
    for k in ("eps", "period", "window"):
        res.pop(k, None)
    res.update(json.loads(key(sh, pr, po, st)))
    with open(OUT, "a") as f:
        f.write(json.dumps(res) + "\n")
    return res


def main(par=4):
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT):
            r = json.loads(l)
            if "error" not in r:
                done.add(key(r["shell"], r["pair"], {k: r[k] for k in ("policy",) } | {k: r[k] for k in ("period", "detect", "setup", "window", "eps") if k in r and k in _po_keys(r["policy"])}, r["start"]))
    jobs = [j for j in itertools.product(SHELLS, PAIRS, POLICIES, STARTS) if key(*j) not in done]
    # interleave policies so partial results cover every policy
    jobs.sort(key=lambda j: (STARTS.index(j[3]), list(SHELLS).index(j[0])))
    if REVERSE:
        jobs.reverse()
    print(len(jobs), "runs to do", flush=True)
    with ThreadPoolExecutor(par) as ex:
        for r in ex.map(run, jobs):
            print(r.get("shell"), r.get("pair"), r.get("policy"), r.get("start"), "ERR" if "error" in r else "ok", flush=True)


def _po_keys(policy):
    return {"snap": ("period",), "ondemand": ("detect", "setup"), "asr": ("window",), "life": ("eps",)}[policy]


REVERSE = "--reverse" in sys.argv

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    main(int(args[0]) if args else 4)
