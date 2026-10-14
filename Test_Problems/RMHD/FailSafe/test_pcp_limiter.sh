#!/usr/bin/env bash
set -euo pipefail
test_src=$(cd "$(dirname "$0")" && pwd)
pluto_src=$(cd "$test_src/../../../Src" && pwd)
test_out=$(mktemp -d /tmp/pluto-rmhd-pcp-limiter.XXXXXX)
for energy in 0 1; do
  cc -std=c17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
    -ffunction-sections -fdata-sections -Wl,--gc-sections \
    -I"$test_src" -I"$pluto_src" -I"$pluto_src/RMHD" -I"$pluto_src/EOS/Ideal" \
    -DRMHD_REDUCED_ENERGY="$energy" \
    "$test_src/test_pcp_limiter.c" "$pluto_src/RMHD/pcp_nr.c" \
    -lm -o "$test_out/test_$energy"
  ASAN_OPTIONS=detect_leaks=0 "$test_out/test_$energy"
done
printf 'Artifacts: %s\n' "$test_out"
