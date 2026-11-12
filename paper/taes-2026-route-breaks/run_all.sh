#!/bin/bash
# Paper 7 pipeline, in order. ns3_patch/build.sh builds the one binary against the existing ns-3 libraries.
set -e
cd "$(dirname "$0")"
python3 law_tables.py          # closed-form law vs sampling -> out/law.json
python3 grid_allF.py           # +Grid range over every phasing factor -> out/grid_allF.json
python3 pair_check.py          # Lemma 1 vs vector geometry -> out/pair_check.json
python3 real_shell_check.py 1700   # law on the real 53 deg Starlink shell -> out/real_shell.json
python3 check_track.py         # Lemma 2 by exact-model Monte Carlo -> out/track_check.json
python3 gen_real_shell.py      # real element sets and +Grid matching -> data/real_*, out/real_shells.json
ns3_patch/build.sh             # bin/p7-route-lifetime-r3 (r2: same sources with the round-3 block removed)
ns3_patch/build_law.sh         # bin/p7-law-check (library propagator only)
(cd ~/6g_ntn_ns3/ns-3-dev && "$OLDPWD"/bin/p7-law-check > "$OLDPWD"/out/law_check.jsonl)
python3 run_p7.py 3            # campaign on bin/p7-route-lifetime-r2 -> out/p7_v8.jsonl
python3 run_p7_r3.py 3         # r3: reproduction, ISASR (Algorithm 3), real-shell elevation sweep, held-out tests -> out/p7_v9.jsonl
python3 offline_opt.py         # dynamic program of (P1) from the reach dumps -> out/offline_opt.json
python3 make_macros_corr.py    # ../macros_p7.tex (correspondence)
python3 make_figs.py           # ../figures/*.pdf
