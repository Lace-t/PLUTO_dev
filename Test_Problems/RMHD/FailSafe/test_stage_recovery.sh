#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
src=$(cd "$here/../../../Src" && pwd)
out=$(mktemp -d /tmp/pluto-stage-recovery.XXXXXX)
for failsafe in 0 1; do for energy in 0 1; do
 cc -std=c17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
  -ffunction-sections -fdata-sections -Wl,--gc-sections \
  -DPCP_FOFC=1 -DUPDATE_VECTOR_POTENTIAL=0 -DPCP_NR=1 \
  -DFAILSAFE="$failsafe" -DRMHD_REDUCED_ENERGY="$energy" \
  -I"$here" -I"$src" -I"$src/RMHD" -I"$src/EOS/Ideal" \
  "$here/test_stage_recovery.c" "$src/RMHD/pcp_stage.c" "$src/mappers3D.c" \
  "$src/RMHD/mappers.c" "$src/RMHD/mappers_loc.c" \
  "$src/RMHD/rmhd_energy_solve.c" "$src/RMHD/rmhd_pressure_fix.c" \
  "$src/RMHD/pcp_nr.c" "$src/failsafe.c" "$src/arrays.c" "$src/rbox.c" \
  -lm -o "$out/test_${failsafe}_${energy}"
 ASAN_OPTIONS=detect_leaks=0 "$out/test_${failsafe}_${energy}"
done; done
printf 'Artifacts: %s\n' "$out"
