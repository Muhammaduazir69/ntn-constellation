#!/bin/bash
# Standalone check of Lemma 1 against the library propagator; touches nothing under ns-3-dev.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
B=/home/uzair/6g_ntn_ns3/ns-3-dev/build
FLAGS="-DNS3_BUILD_PROFILE_OPTIMIZED -D__LINUX__ -I$B/include -iquote $B/include/ns3 -O2 -DNDEBUG -std=c++20 -fPIE -march=native -mtune=native"
LIBS=$(ls $B/lib/libns3.43-*-optimized.so | grep -E "libns3.43-(ntn-constellation|mobility|network|core|stats|bridge|internet|traffic-control)-optimized" | sort)
nice -n 10 g++ $FLAGS "$HERE/law_check.cc" -o "$HERE/../bin/p7-law-check" $LIBS -Wl,-rpath,$B/lib -Wl,-rpath-link,$B/lib
sha256sum "$HERE/../bin/p7-law-check"
