"""Lemma 1 against direct vector geometry on random pairs of the Starlink shell; writes out/pair_check.json."""
import json, numpy as np, isl_core as C
sh = C.walker(72, 22, 39, 53, 550e3); rng = np.random.default_rng(1); err = 0.0; n = 2000
for _ in range(n):
    a, b = rng.integers(0, 72 * 22, 2); u = rng.uniform(0, 2 * np.pi)
    X = C.pos(sh, u); d = np.linalg.norm(X[a] - X[b])
    M, N = C.MN(sh["Om"][b] - sh["Om"][a], sh["u0"][b] - sh["u0"][a], sh["i"])
    err = max(err, abs(d - sh["r"] * np.sqrt(max(2 - 2 * (M - N * np.cos(2 * (sh["u0"][a] + u) + sh["u0"][b] - sh["u0"][a])), 0))))
json.dump(dict(pairs=n, max_err_m=err), open("out/pair_check.json", "w")); print(n, err)
