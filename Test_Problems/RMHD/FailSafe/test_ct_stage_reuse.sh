#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
src=$(cd "$here/../../../Src" && pwd)
out=$(mktemp -d /tmp/pluto-ct-stage-reuse.XXXXXX)
sed 's/#define NTRACER 0/#define NTRACER 1/' "$here/definitions.h" > "$out/definitions.h"
for rk in RK2 RK3; do for energy in 0 1; do
 cc -std=c17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
  -ffunction-sections -fdata-sections -Wl,--gc-sections \
  -Wl,--wrap=ConsToPrim3D -Wl,--wrap=RMHD_PCP_ConsToPrim3D_CheckModified \
  -DPCP_FOFC=1 -DUPDATE_VECTOR_POTENTIAL=0 -DCT_RMHD_CONSISTENT=1 \
  -DCT_RMHD_DIAGNOSTICS=0 -DPCP_NR=1 -DTIME_STEPPING="$rk" -DRMHD_REDUCED_ENERGY="$energy" \
  -I"$out" -I"$src" -I"$src/RMHD" -I"$src/EOS/Ideal" \
  "$here/test_ct_stage_reuse.c" "$src/Time_Stepping/rk_step.c" \
  "$src/MHD/CT/ct_field_average.c" "$src/MHD/CT/ct_rmhd_correct.c" \
  "$src/RMHD/pcp_stage.c" "$src/RMHD/pcp_nr.c" "$src/mappers3D.c" \
  "$src/RMHD/mappers.c" "$src/RMHD/mappers_loc.c" \
  "$src/RMHD/rmhd_energy_solve.c" "$src/RMHD/rmhd_pressure_fix.c" \
  "$src/failsafe.c" "$src/arrays.c" "$src/rbox.c" "$src/EOS/Ideal/eos.c" -lm -o "$out/test"
 stages=2; if [[ $rk == RK3 ]]; then stages=3; fi
 for stage in $(seq 0 "$stages"); do
  ASAN_OPTIONS=detect_leaks=0 "$out/test" "$stage"
 done
done; done
printf 'Artifacts: %s\n' "$out"
