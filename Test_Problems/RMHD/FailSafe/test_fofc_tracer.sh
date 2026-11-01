#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
src=$(cd "$here/../../../Src" && pwd)
out=$(mktemp -d /tmp/pluto-fofc-tracer.XXXXXX)
sed -e 's/#define NTRACER 0/#define NTRACER 1/' \
    -e 's/#define DIVB_CONTROL CONSTRAINED_TRANSPORT/#define DIVB_CONTROL DIV_CLEANING/' \
    "$here/definitions.h" > "$out/definitions.h"
cc -std=c17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
  -ffunction-sections -fdata-sections -Wl,--gc-sections \
  -DPCP_FOFC=YES -DUPDATE_VECTOR_POTENTIAL=NO \
  -I"$out" -I"$src" -I"$src/RMHD" -I"$src/EOS/Ideal" \
  "$here/test_fofc_tracer.c" "$src/RMHD/pcp_stage.c" \
  "$src/adv_flux.c" "$src/arrays.c" -lm -o "$out/test"
ASAN_OPTIONS=detect_leaks=0 "$out/test"

