#!/usr/bin/env bash
set -euo pipefail
test_src=$(cd "$(dirname "$0")" && pwd)
pluto_src=$(cd "$test_src/../../../Src" && pwd)
test_out=$(mktemp -d /tmp/pluto-local-ct.XXXXXX)
cc -std=c17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
  -ffunction-sections -fdata-sections -Wl,--gc-sections \
  -I"$test_src" -I"$pluto_src" -I"$pluto_src/RMHD" -I"$pluto_src/EOS/Ideal" \
  "$test_src/test_local_ct.c" "$pluto_src/arrays.c" -lm -o "$test_out/test"
ASAN_OPTIONS=detect_leaks=0 "$test_out/test"
printf 'Artifacts: %s\n' "$test_out"
