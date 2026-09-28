#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
src=$(cd "$here/../../../Src" && pwd)
out=$(mktemp -d /tmp/pluto-split-recovery.XXXXXX)
for divb in CONSTRAINED_TRANSPORT DIV_CLEANING; do
 sed -e 's/#define NTRACER 0/#define NTRACER 1/' \
     -e "s/#define DIVB_CONTROL CONSTRAINED_TRANSPORT/#define DIVB_CONTROL $divb/" \
     "$here/definitions.h" > "$out/definitions.h"
 for energy in 0 1; do for strict in 0 1; do for nr in 0 1; do
  cc -std=c17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
   -ffunction-sections -fdata-sections -Wl,--gc-sections \
   -Wl,--wrap=RMHD_EnergySolve -Wl,--wrap=RMHD_PCP_NR -Wl,--wrap=RMHD_PressureFix \
   -DPCP_FOFC=1 -DUPDATE_VECTOR_POTENTIAL=0 -DPCP_NR="$nr" \
   -DRMHD_REDUCED_ENERGY="$energy" -DRMHD_FAILSAFE_PCP_STRICT="$strict" \
   -I"$out" -I"$src" -I"$src/RMHD" -I"$src/EOS/Ideal" \
   "$here/test_split_recovery.c" "$src/RMHD/mappers.c" "$src/RMHD/mappers_loc.c" \
   "$src/RMHD/rmhd_energy_solve.c" "$src/RMHD/rmhd_pressure_fix.c" \
   "$src/RMHD/pcp_nr.c" "$src/RMHD/pcp_stage.c" "$src/mappers3D.c" \
   "$src/arrays.c" "$src/rbox.c" -lm -o "$out/test"
  ASAN_OPTIONS=detect_leaks=0 "$out/test"
 done; done; done
done
printf 'PASS 16 configurations. Artifacts: %s\n' "$out"
