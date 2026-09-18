#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
src=$(cd "$here/../../../Src" && pwd)
out=$(mktemp -d /tmp/pluto-ct-correction.XXXXXX)
sed 's/#define NTRACER 0/#define NTRACER 1/' "$here/definitions.h" > "$out/definitions.h"
printf '\n#define PCP_NR YES\n' >> "$out/definitions.h"
for energy in 0 1; do
  cc -std=c17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
    -ffunction-sections -fdata-sections -Wl,--gc-sections \
    -I"$out" -I"$src" -I"$src/RMHD" -I"$src/EOS/Ideal" \
    -DRMHD_REDUCED_ENERGY="$energy" \
    "$here/test_ct_correction.c" "$src/MHD/CT/ct_rmhd_correct.c" \
    "$src/RMHD/mappers_loc.c" "$src/RMHD/rmhd_energy_solve.c" \
    "$src/RMHD/rmhd_pressure_fix.c" "$src/RMHD/pcp_nr.c" -lm \
    -o "$out/test_$energy"
  ASAN_OPTIONS=detect_leaks=0 "$out/test_$energy"
done
printf 'Artifacts: %s\n' "$out"
