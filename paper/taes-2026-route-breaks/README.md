# Code and data: route breaks and lifetime-aware routing over +Grid laser links (IEEE TAES, 2026)

Companion material for M. Uzair, "Route Breaks and Lifetime-Aware Routing in Low Earth Orbit Constellations With Range-Limited +Grid Laser Links," correspondence submitted to IEEE Transactions on Aerospace and Electronic Systems, 2026.

The routing policies (longest-lived-route search, elevation ranking, ISASR, merged snapshots) and the real-TLE shells are in this module (`model/contact-graph-router.*`, `model/contact-graph-scheduler.*`). This folder holds the paper pipeline.

| Path | Contents |
|---|---|
| `run_all.sh` | The whole pipeline, in order |
| `ns3_patch/` | Source of record of the campaign example and law check, with `build.sh` and `build_law.sh` |
| `law_tables.py`, `grid_allF.py`, `pair_check.py`, `check_track*.py`, `real_shell_check.py` | Closed-form range law, +Grid range cap for every phasing factor, and their checks |
| `run_p7.py`, `run_p7_r3.py`, `run_ns3_routes*.py`, `offline_opt.py` | Packet-level ns-3 campaigns and the offline optimum |
| `data/` | Public Starlink and OneWeb element sets (CelesTrak) and the matched +Grid edge lists |
| `out/` | Every result file the paper uses (`p7_v8.jsonl`, `p7_v9.jsonl`, `law*.json`, `real_shells.json`, `reach*/` dumps, statistics) |
| `make_macros_corr.py`, `make_figs.py` | Numbers and figures of the paper, generated from `out/` |

## Reproduce
`./run_all.sh` (ns-3.43 with this module in `contrib/`). To regenerate only the numbers and figures from the stored results, run `python3 make_macros_corr.py` and `python3 make_figs.py`.

License: as the repository. The element sets are public CelesTrak data.
Contact: Muhammad Uzair, muhammaduzairr69@gmail.com, ORCID 0009-0002-4104-2680.
