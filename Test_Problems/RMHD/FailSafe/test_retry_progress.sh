#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
src=$(cd "$here/../../../Src" && pwd)
out=$(mktemp -d /tmp/pluto-retry-progress.XXXXXX)
flags=(-std=c17 -O1 -g -ffunction-sections -fdata-sections -Wl,--gc-sections
 -DPCP_FOFC=1 -DUPDATE_VECTOR_POTENTIAL=0
 -I"$here" -I"$src" -I"$src/RMHD" -I"$src/EOS/Ideal")
sources=("$here/test_retry_progress.c" "$src/RMHD/pcp_stage.c" "$src/failsafe.c" "$src/arrays.c")
for passes in 1 4; do
 cc "${flags[@]}" -DRMHD_FAILSAFE_LOCAL_PASSES="$passes" \
  -fsanitize=address,undefined -fno-omit-frame-pointer "${sources[@]}" -lm -o "$out/test"
 for mode in 0 1 2 3 4; do
  mkdir "$out/p${passes}m${mode}"
  status=0
  (cd "$out/p${passes}m${mode}"; ASAN_OPTIONS=detect_leaks=0 "$out/test" "$mode") > "$out/p${passes}m${mode}.log" 2>&1 || status=$?
  expected=0; if [[ $mode == 4 ]]; then expected=1; fi
  if [[ $status != "$expected" ]]; then cat "$out/p${passes}m${mode}.log"; exit 1; fi
 done
done
if [[ ${RUN_MPI:-0} == 1 ]]; then
 mpicc "${flags[@]}" -DPARALLEL -I"$src/Parallel" "${sources[@]}" -lm -o "$out/mpi"
 for mode in 0 1 2 3; do
  (cd "$out"; mpirun -np 2 ./mpi "$mode")
 done
fi
printf 'PASS retry progress (10 serial cases). Artifacts: %s\n' "$out"
