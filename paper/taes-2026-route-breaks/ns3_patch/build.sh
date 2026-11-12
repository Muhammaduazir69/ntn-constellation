#!/bin/bash
# Standalone build of the paper's binary; touches nothing under ns-3-dev. The patched scheduler and router are
# compiled under other class names, so the library keeps its own classes and this executable uses its own.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
B=/home/uzair/6g_ntn_ns3/ns-3-dev/build
OUT="$HERE/../bin"
FLAGS="-DEIGEN_MPL2_ONLY -DHAVE_BOOST -DHAVE_BOOST_UNITS -DHAVE_EIGEN3 -DHAVE_GSL -DHAVE_LIBXML2 -DHAVE_PACKET_H -DHAVE_SQLITE3 -DNS3_BUILD_PROFILE_OPTIMIZED -D__LINUX__ -DP7_STANDALONE -DContactGraphScheduler=ContactGraphSchedulerP7 -DContactGraphRouter=ContactGraphRouterP7 -I$HERE -I$B/include -iquote $B/include/ns3 -I/usr/include/eigen3 -I/usr/include/libxml2 -O3 -DNDEBUG -std=c++20 -fPIE -fno-semantic-interposition -march=native -mtune=native -Wno-misleading-indentation"
LIBS=$(ls $B/lib/libns3.43-*-optimized.so | grep -v -- "-test-" | grep -E "libns3.43-(ntn-constellation|ntn-cho|applications|internet|point-to-point|bridge|traffic-control|network|mobility|stats|core)-optimized" | sort)
T=$(mktemp -d "$OUT/obj.XXXX")
for f in contact-graph-scheduler contact-graph-router ntn-constellation-route-lifetime; do
  nice -n 19 g++ $FLAGS -c "$HERE/$f.cc" -o "$T/$f.o"
done
nice -n 19 g++ -o "$OUT/p7-route-lifetime-r3" "$T"/*.o $LIBS -Wl,-rpath,$B/lib -Wl,-rpath-link,$B/lib
rm -rf "$T"
sha256sum "$OUT/p7-route-lifetime-r3"
