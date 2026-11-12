"""Offline optimum of (P1) from the per-second reach dumps of the lifetime-aware policy (eps = 0.05).

At second t the longest admissible route lives L(t); a route chosen at t can be kept until
R(t) = t + L(t) - g - t_a, when its replacement must start. The greedy policy of Algorithm 1 switches at
R(t) every time. The optimum may switch earlier to a route that reaches further:
    f(t) = 0 if R(t) >= t_end,   f(t) = 1 + min_{t < t' <= R(t)} f(t')   otherwise,
the fewest further switches from a route chosen at t (dynamic programming backward on the 1 s grid).
Writes out/offline_opt.json with changes per hour of greedy and optimum per run, and the floor."""
import csv, glob, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
G, TA = 0.5, 1.0


def solve(path):
    rows = list(csv.DictReader(open(path)))
    t = np.array([float(r["t"]) for r in rows])
    L = np.array([float(r["lifetime_s"]) for r in rows])
    n = len(t)
    if n < 10:
        return None
    t_end = t[-1]
    R = t + L - G - TA
    # index of the last grid point <= R(t)
    reach_idx = np.searchsorted(t, R, side="right") - 1
    # a route that cannot be sustained for even one grid step moves on at the next step
    reach_idx = np.maximum(reach_idx, np.minimum(np.arange(n) + 1, n - 1))
    INF = 10 ** 9
    f = np.full(n, INF, dtype=np.int64)
    for k in range(n - 1, -1, -1):
        if R[k] >= t_end:
            f[k] = 0
            continue
        j = reach_idx[k]
        if j > k:
            f[k] = (1 if R[k] > t[k] else 0) + f[k + 1:j + 1].min()  # a route that cannot be sustained is not a switch yet
    # greedy: from the first point, switch at R each time
    k, greedy = 0, 0
    while R[k] < t_end:
        j = reach_idx[k]
        if j <= k:
            greedy = None
            break
        greedy += 1 if R[k] > t[k] else 0
        k = j
    hours = (t_end - t[0]) / 3600.0
    return dict(opt_per_h=float(f[0]) / hours if f[0] < INF else None,
                greedy_per_h=(greedy / hours) if greedy is not None else None)


if __name__ == "__main__":
    res = {}
    for p in sorted(glob.glob(os.path.join(HERE, "out", "reach8", "*.csv"))):
        r = solve(p)
        if r:
            res[os.path.basename(p)[:-4]] = r
            print(os.path.basename(p), r)
    json.dump(res, open(os.path.join(HERE, "out", "offline_opt.json"), "w"), indent=1)
