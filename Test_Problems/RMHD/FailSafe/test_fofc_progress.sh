#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
src=$(cd "$here/../../../Src" && pwd)
out=$(mktemp -d /tmp/pluto-fofc-progress.XXXXXX)
sed 's/#define DIVB_CONTROL CONSTRAINED_TRANSPORT/#define DIVB_CONTROL DIV_CLEANING/' \
  "$here/definitions.h" > "$out/definitions.h"
flags=(-std=c17 -O1 -g -ffunction-sections -fdata-sections -Wl,--gc-sections
  -DPCP_FOFC=YES -DUPDATE_VECTOR_POTENTIAL=NO -DPCP_PRESSURE_FIX=YES
  -I"$out" -I"$src" -I"$src/RMHD" -I"$src/EOS/Ideal")
sources=("$here/test_fofc_progress.c" "$src/RMHD/pcp_stage.c"
  "$src/RMHD/pcp_nr.c" "$src/arrays.c" "$src/rbox.c")
for cap in default 4; do
cap_flags=()
if [[ $cap != default ]]; then cap_flags=(-DPCP_STAGE_LOCAL_PASSES="$cap"); fi
cc "${cap_flags[@]}" "${flags[@]}" -fsanitize=address,undefined -fno-omit-frame-pointer \
  "${sources[@]}" -lm -o "$out/serial"
for mode in 0 1 2 3 4 5; do
  ASAN_OPTIONS=detect_leaks=0 "$out/serial" "$mode"
done
if [[ ${RUN_MPI:-0} == 1 ]]; then
  mpicc "${cap_flags[@]}" "${flags[@]}" -DPARALLEL -I"$src/Parallel" \
    "${sources[@]}" -lm -o "$out/mpi"
  for mode in 0 1 2 3 4 5; do mpirun -np 2 "$out/mpi" "$mode"; done
fi

done
