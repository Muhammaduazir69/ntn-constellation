"""ns-3 campaign v5 of the TAES study (after the second review). Resumable, one jsonl line per run.

Changes from v4: six start instants drawn at random over ten days (seeded), not commensurate with the shells'
plane spacing; 45 min of traffic per run (4 to 9 ground passes), so change rates are not quantized by the
horizon; a measurement-only elevation-trigger baseline (elev: switch when an end satellite sinks below the
25 deg up threshold, while the old link stays up to 23 deg); merged snapshots swept over 30, 120 and 300 s;
an acquisition-time sweep. Parts:
  main   all shells, London-New York and Sydney-Los Angeles; snap 1 s, on-demand, elev, merged 30/120/300 s,
         predictive shortest (life eps 0), life eps 0.05 and 0.10, mrv; tacq = 1 s
  tacq   tacq 0, 0.2, 1 (main) and 3 s; Starlink London-New York; snap 1 s, on-demand, elev, merged 120 s, life 0.05
  pred   guard 0.5 and 2 s x along-track ephemeris error 5 and 20 km for life 0.05; merged 120 s with 5 and 20 km;
         Starlink and OneWeb, London-New York
  cross  OneWeb with cross-plane links off above 70 deg latitude, London-New York; on-demand, elev, merged 120 s, life 0.05
"""
import json, os, subprocess, sys, threading
import numpy as np
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
NS3 = os.path.expanduser("~/6g_ntn_ns3/ns-3-dev")
BIN = os.path.join(NS3, "build/contrib/ntn-constellation/examples/ns3.43-ntn-constellation-route-lifetime-optimized")
OUT = os.path.join(HERE, "out", "ns3_routes6.jsonl")
DUR = 2705.0
SHELLS = {"starlink": dict(planes=72, perPlane=22, phasing=39, inc=53.0, altKm=550, star=0),
          "kuiper": dict(planes=34, perPlane=34, phasing=1, inc=51.9, altKm=630, star=0),
          "oneweb": dict(planes=12, perPlane=49, phasing=1, inc=87.9, altKm=1200, star=1)}
PAIRS = {"lon-nyc": (51.51, -0.13, 40.71, -74.01), "syd-lax": (-33.87, 151.21, 34.05, -118.24)}
STARTS = sorted(int(x) for x in np.random.default_rng(20261030).uniform(0, 10 * 86400, 6))
MAIN = [dict(policy="snap", period=1.0), dict(policy="ondemand"), dict(policy="elev"),
        dict(policy="asr", window=30.0), dict(policy="asr", window=120.0), dict(policy="asr", window=300.0),
        dict(policy="life", eps=0.0), dict(policy="life", eps=0.05), dict(policy="life", eps=0.10), dict(policy="mrv")]
TACQ = [dict(policy="snap", period=1.0), dict(policy="ondemand"), dict(policy="elev"), dict(policy="asr", window=120.0),
        dict(policy="life", eps=0.05)]


def jobs():
    """v6: rerun the policies whose code changed after the third review (merged snapshots repair a break and fall
    back to the longest-lived route; lifetime-aware search no longer polls when no route outlives the guard) and
    add the elevation trigger with a delay budget (elevb). Unchanged policies are taken from v5."""
    MAIN6 = [dict(policy="asr", window=30.0), dict(policy="asr", window=120.0), dict(policy="asr", window=300.0),
             dict(policy="elevb", eps=0.05), dict(policy="elevb", eps=0.10),
             dict(policy="life", eps=0.0), dict(policy="life", eps=0.05), dict(policy="life", eps=0.10)]
    J = []
    for st in STARTS:
        for sh in SHELLS:
            for pr in PAIRS:
                for po in MAIN6:
                    J.append(("main", sh, pr, st, dict(po, tAcq=1.0)))
        for ta in (0.0, 0.2, 3.0):
            for po in (dict(policy="asr", window=120.0), dict(policy="elevb", eps=0.05), dict(policy="life", eps=0.05)):
                J.append(("tacq", "starlink", "lon-nyc", st, dict(po, tAcq=ta)))
        for sh in ("starlink", "oneweb"):
            for gd in (0.5, 2.0):
                for sg in (5.0, 20.0):
                    J.append(("pred", sh, "lon-nyc", st, dict(policy="life", eps=0.05, tAcq=1.0, guard=gd, predSigmaKm=sg)))
            for sg in (5.0, 20.0):
                J.append(("pred", sh, "lon-nyc", st, dict(policy="asr", window=120.0, tAcq=1.0, predSigmaKm=sg)))
        for po in (dict(policy="elevb", eps=0.05), dict(policy="asr", window=120.0), dict(policy="life", eps=0.05)):
            J.append(("cross", "oneweb", "lon-nyc", st, dict(po, tAcq=1.0, crossLatDeg=70.0)))
    return J


def key(j):
    return json.dumps([j[0], j[1], j[2], j[3], sorted(j[4].items())])


def run(j):
    part, sh, pr, st, po = j
    a = dict(SHELLS[sh]); la, loa, lb, lob = PAIRS[pr]
    a.update(latA=la, lonA=loa, latB=lb, lonB=lob, duration=DUR, tick=0.1, guard=0.5, startOffset=st,
             RngRun=1 + STARTS.index(st))
    a.update(po)
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
