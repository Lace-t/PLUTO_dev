#!/usr/bin/env bash
set -euo pipefail
test_src=$(cd "$(dirname "$0")" && pwd)
pluto_src=$(cd "$test_src/../../../Src" && pwd)
test_out=$(mktemp -d /tmp/pluto-stage-mapper.XXXXXX)
flags=(-std=c17 -O1 -g -ffunction-sections -fdata-sections -Wl,--gc-sections
 -DPCP_STAGE_LIMITING=1 -DUPDATE_VECTOR_POTENTIAL=0
 -I"$test_src" -I"$pluto_src" -I"$pluto_src/RMHD" -I"$pluto_src/EOS/Ideal")
sources=("$test_src/test_stage_mapper.c" "$pluto_src/mappers3D.c"
 "$pluto_src/RMHD/pcp_stage.c" "$pluto_src/RMHD/pcp_nr.c"
 "$pluto_src/arrays.c" "$pluto_src/rbox.c")
for energy in 0 1; do
 for strict in 0 1; do
  cc "${flags[@]}" -DRMHD_REDUCED_ENERGY="$energy" -DRMHD_FAILSAFE_PCP_STRICT="$strict" \
   -fsanitize=address,undefined -fno-omit-frame-pointer "${sources[@]}" -lm -o "$test_out/serial"
  ASAN_OPTIONS=detect_leaks=0 "$test_out/serial"
 done
done
if [[ "${RUN_MPI:-0}" == 1 ]]; then
 mpicc "${flags[@]}" -DPARALLEL -I"$pluto_src/Parallel" "${sources[@]}" -lm -o "$test_out/mpi"
 mpirun -np 2 "$test_out/mpi"
fi
printf 'Artifacts: %s\n' "$test_out"
