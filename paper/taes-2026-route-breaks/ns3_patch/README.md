# Paper 7 patch to the ntn-constellation module (Nov 2026)

`ns-3-dev` is shared with running campaigns, so this revision did not touch it or its build. All ns-3 data of the
article come from one standalone binary, `../bin/p7-route-lifetime`, built by `build.sh` from the files here against
the existing optimized libraries in `ns-3-dev/build/lib`. Its sha256 is recorded in every row of `../out/p7_v7.jsonl`.

- `orig/`: the five source files as copied from `contrib/ntn-constellation` on 4 Nov 2026.
- `contact-graph-{scheduler,router}.{h,cc}`, `ntn-constellation-route-lifetime.cc`: patched copies.
  `p7_route_lifetime.patch` is the diff against `orig/`. `apply_edits.py` and `edit_example.py` produce them from `orig/`.

The patched classes compile under other names (`-DContactGraphScheduler=ContactGraphSchedulerP7`,
`-DContactGraphRouter=ContactGraphRouterP7`, and `P7_STANDALONE` renames their log components and TypeIds), because a
second registration under the same name is fatal and `-fno-semantic-interposition` makes an override unreliable. There is
no symbol override: the library keeps its classes and the executable uses its own.

What changed.
1. Speed, no change of meaning: ground links sampled only for candidate satellites, skipping a satellite until it could
   reach the mask at a 2 deg/s elevation-rate bound, on the same 0.1 s grid; ISL ranges refreshed every `--islRefresh`
   seconds (5 s in the campaign; every tick with a latitude cutoff, where they do drop); flat-array Dijkstra with the same
   relaxation order. About 200 s of CPU per run falls to about 8 s.
   Checked against the shared-build rows: snapshot Starlink London-New York start 169873 has identical changes (130.667/h)
   and loss (0.0267407); lifetime eps 0.05 has identical changes, Dijkstra count, received count and mean delay, and
   mean stretch 2.296 % against 2.299 %.
2. Admissibility (the M9 bug): a route must outlive g + 2 t_a, not g + t_a. A route with lifetime in between needed its
   successor's acquisition to start while its own was still pending, and `Switch` drops a decision during an acquisition.
   That caused the 21 lost packets on Starlink Sydney-Los Angeles start 218123; with the fix that run loses none.
3. New policies `asrk`, `track`, `alpr`, `isasr`; bounded fallback `--fallbackEps`; one-terminal mode `--bbm`;
   control-plane mode `--noData`; real shells `--tle --edges` with Vallado SGP4; counters for ingress, egress, ISL-only
   changes, acquisitions, breaks by end, route age, residual at install, search CPU, stretch quantiles and CPU time.

**To merge at the shared build:** apply `p7_route_lifetime.patch` under `contrib/ntn-constellation` (model and example
paths), drop the `P7_STANDALONE` blocks, change the two includes in the example back to `"ns3/contact-graph-*.h"`,
then rebuild. Nothing else depends on the rename; other users of the scheduler keep the full scan unless they call
`SetGslCandidates` and `SetIslRefreshEvery`. The router's weight map is no longer refreshed by `UpdateEdgeWeight`, so
`ShortestPathWeightedRegenOnly` must be moved to the flat adjacency when merging.

## Round 2 (Nov 2026): binary `../bin/p7-route-lifetime-r2`, rows in `../out/p7_v8.jsonl`

`build.sh` now builds `p7-route-lifetime-r2` at nice 19; `edit_example.py` holds the round-2 edits after the marker "round 2 of review".
1. Residual cache (C1): a link that does not end within the 1800 s horizon is remembered as saturated with the time of that
   prediction and re-predicted after 900 s, instead of being cached as `now + 1800` and aging to a false end. Checked on the reach
   dump of Kuiper Sydney-Los Angeles start 218123: the lifetime no longer falls to 1.2 s at t = 1800 s.
2. ISASR (C2): activeness cost is state per edge as in Algorithm 3 of Bhattacharjee et al.; links that do not end within the horizon
   carry no cost, as the article states.
3. Elevation noise (`--elevNoiseDeg`), least-squares track over `--trackSamples` readings, elevation-rank floor grid (`--floorGrid`),
   merged-snapshot fallback share, track-estimate failures.
The merge note above still applies; the patch file is regenerated from `orig/`.
