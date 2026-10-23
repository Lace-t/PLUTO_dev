#!/usr/bin/env bash
set -euo pipefail
test_src=$(cd "$(dirname "$0")" && pwd)
pluto_src=$(cd "$test_src/../../../Src" && pwd)
test_out=$(mktemp -d /tmp/pluto-flux-cache.XXXXXX)
sed 's/#define DIVB_CONTROL CONSTRAINED_TRANSPORT/#define DIVB_CONTROL DIV_CLEANING/' \
  "$test_src/definitions.h" > "$test_out/definitions.h"
cc -std=c17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
  -ffunction-sections -fdata-sections -Wl,--gc-sections \
  -DPCP_FOFC=YES -DUPDATE_VECTOR_POTENTIAL=NO \
  -I"$test_out" -I"$pluto_src" -I"$pluto_src/RMHD" -I"$pluto_src/EOS/Ideal" \
  "$test_src/test_flux_cache.c" "$pluto_src/RMHD/pcp_stage.c" \
  "$pluto_src/arrays.c" -lm -o "$test_out/test"
ASAN_OPTIONS=detect_leaks=0 "$test_out/test"
printf 'Artifacts: %s\n' "$test_out"
