#!/usr/bin/env python3
"""ns-3 campaign v4 of the TAES study (after the first hostile review). Resumable, one jsonl line per run.

Every run uses a common ground-link acquisition time (tacq) for all policies, a 0.5 s switch guard, a random
phase for snapshot clocks, counters from the start of traffic, and four start instants a quarter day apart (three for the stress parts), 15 min of traffic per run, so
that the Earth's rotation varies the geometry. Parts:
  main     all shells, pairs and starts; snap 1 s (0.1 s on Starlink London-New York only), on-demand, merged 30/120 s, life eps 0/0.01/0.05/0.10,
           life eps 0.05 with the hold constraint, maximum remaining visibility (mrv); tacq = 1 s
  tacq     tacq 0 and 3 s; Starlink and OneWeb, London-New York; snap 1 s, on-demand, merged 30 s, life 0.05
  pred     guard 0.1/0.5/2 s and along-track ephemeris error 0/5/20 km for life 0.05; error 5/20 km for merged 30 s
  cross    OneWeb with cross-plane links off above 70 deg latitude; snap 1 s, on-demand, merged 30 s, life 0.05
  reach    life 0.05 with the per-second dump of the longest admissible lifetime (offline optimum of (P1))
"""
import json, os, subprocess, sys, threading
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
NS3 = os.path.expanduser("~/6g_ntn_ns3/ns-3-dev")
BIN = os.path.join(NS3, "build/contrib/ntn-constellation/examples/ns3.43-ntn-constellation-route-lifetime-optimized")
OUT = os.path.join(HERE, "out", "ns3_routes4.jsonl")
REACH = os.path.join(HERE, "out", "reach")
SHELLS = {"starlink": dict(planes=72, perPlane=22, phasing=39, inc=53.0, altKm=550, star=0),
          "kuiper": dict(planes=34, perPlane=34, phasing=1, inc=51.9, altKm=630, star=0),
          "oneweb": dict(planes=12, perPlane=49, phasing=1, inc=87.9, altKm=1200, star=1)}
PAIRS = {"lon-nyc": (51.51, -0.13, 40.71, -74.01), "fra-sin": (50.11, 8.68, 1.35, 103.82),
         "syd-lax": (-33.87, 151.21, 34.05, -118.24)}
STARTS = [0, 21600, 43200, 64800]          # four geometries a quarter day apart
MAIN = [dict(policy="snap", period=1.0), dict(policy="ondemand"),
        dict(policy="asr", window=30.0), dict(policy="asr", window=120.0),
        dict(policy="life", eps=0.0), dict(policy="life", eps=0.01), dict(policy="life", eps=0.05),
        dict(policy="life", eps=0.10), dict(policy="life", eps=0.05, hold=1), dict(policy="mrv")]
CORE = [dict(policy="snap", period=1.0), dict(policy="ondemand"), dict(policy="asr", window=30.0),
        dict(policy="life", eps=0.05)]


def jobs():
    J = []
    for st in STARTS:
        for sh in SHELLS:
            for pr in PAIRS:
                for po in MAIN:
                    J.append(("main", sh, pr, st, dict(po, tAcq=1.0)))
        # the 100 ms clock as a reference point only (it costs 36,000 route computations an hour)
        J.append(("main", "starlink", "lon-nyc", st, dict(policy="snap", period=0.1, tAcq=1.0)))
    for st in STARTS[:3]:
        for sh in ("starlink", "oneweb"):
            for ta in (0.0, 3.0):
                for po in CORE:
                    J.append(("tacq", sh, "lon-nyc", st, dict(po, tAcq=ta)))
            for gd in (0.1, 0.5, 2.0):
                for sg in (0.0, 5.0, 20.0):
                    if gd == 0.5 and sg == 0.0:
                        continue  # identical to the main run
                    J.append(("pred", sh, "lon-nyc", st, dict(policy="life", eps=0.05, tAcq=1.0, guard=gd, predSigmaKm=sg)))
            for sg in (5.0, 20.0):
                J.append(("pred", sh, "lon-nyc", st, dict(policy="asr", window=30.0, tAcq=1.0, predSigmaKm=sg)))
        for pr in PAIRS:
            for po in CORE:
                J.append(("cross", "oneweb", pr, st, dict(po, tAcq=1.0, crossLatDeg=70.0)))
    for st in STARTS[:2]:
        for sh in SHELLS:
            for pr in PAIRS:
                J.append(("reach", sh, pr, st, dict(policy="life", eps=0.05, tAcq=1.0)))
    return J


def key(j):
    return json.dumps([j[0], j[1], j[2], j[3], sorted(j[4].items())])


def run(j):
    part, sh, pr, st, po = j
    a = dict(SHELLS[sh]); la, loa, lb, lob = PAIRS[pr]
    a.update(latA=la, lonA=loa, latB=lb, lonB=lob, duration=905.0, tick=0.1, guard=0.5, startOffset=st,
             RngRun=1 + STARTS.index(st))
    a.update(po)
    if part == "reach":
        os.makedirs(REACH, exist_ok=True)
        a["reachOut"] = os.path.join(REACH, f"{sh}_{pr}_{st}.csv")
    cmd = ["nice", "-n", "5", BIN] + [f"--{k}={v}" for k, v in a.items()]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=10800, cwd=NS3)
        res = json.loads(r.stdout.strip().splitlines()[-1])
    except Exception as e:
        res = dict(error=f"{type(e).__name__}: {str(e)[:200]}")
    res.update(part=part, shell=sh, pair=pr, start=st, cfg=sorted(po.items()))
    return res


def main(par):
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT):
            r = json.loads(l)
            if "error" not in r:
                done.add(json.dumps([r["part"], r["shell"], r["pair"], r["start"], r["cfg"]]))
    todo = [j for j in jobs() if key(j) not in done]
    print(len(todo), "runs to do", flush=True)
    lock = threading.Lock()

    def work(j):
        r = run(j)
        with lock:
            open(OUT, "a").write(json.dumps(r) + "\n")
        print(j[0], j[1], j[2], j[3], "ERR" if "error" in r else "ok", flush=True)

    with ThreadPoolExecutor(par) as ex:
        list(ex.map(work, todo))


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
