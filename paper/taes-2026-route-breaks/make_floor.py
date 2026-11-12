"""Recomputation floor per shell, the same number as the macro (changes per hour)."""
import json, os
import isl_core as _C
from law_tables import pass_max
HERE = os.path.dirname(os.path.abspath(__file__))
L = json.load(open(os.path.join(HERE, "out", "law.json")))
FLOOR = {}
for s in ("starlink", "kuiper", "oneweb"):
    d = L[s]
    sh = _C.walker(d["P"], d["S"], d["F"], d["inc"], d["alt_km"] * 1e3, s == "oneweb")
    FLOOR[s] = 3600 / max(pass_max(sh, lat, 23.0) for lat in (51.51, 40.71, -33.87, 34.05))
