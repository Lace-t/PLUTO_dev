#!/usr/bin/env bash
set -euo pipefail
test_src=$(cd "$(dirname "$0")" && pwd)
pluto_src=$(cd "$test_src/../../../Src" && pwd)
test_out=$(mktemp -d /tmp/pluto-local-stage.XXXXXX)
flags=(-std=c17 -O1 -g -ffunction-sections -fdata-sections -Wl,--gc-sections
 -DPCP_STAGE_LIMITING=1 -DUPDATE_VECTOR_POTENTIAL=0
 -I"$test_src" -I"$pluto_src" -I"$pluto_src/RMHD" -I"$pluto_src/EOS/Ideal")
sources=("$test_src/test_local_stage.c" "$pluto_src/Time_Stepping/rk_step.c"
 "$pluto_src/RMHD/pcp_stage.c" "$pluto_src/RMHD/pcp_nr.c"
 "$pluto_src/failsafe.c" "$pluto_src/arrays.c" "$pluto_src/rbox.c")
for rk in RK2 RK3; do
cc -DTIME_STEPPING="$rk" "${flags[@]}" -fsanitize=address,undefined -fno-omit-frame-pointer "${sources[@]}" -lm -o "$test_out/serial"
targets=(1 2)
if [[ "$rk" == RK3 ]]; then targets+=(3); fi
for target in "${targets[@]}"; do for mode in 0 1 2; do
 ASAN_OPTIONS=detect_leaks=0 "$test_out/serial" "$mode" "$target"
done; done
done
if [[ "${RUN_MPI:-0}" == 1 ]] && command -v mpicc >/dev/null; then
 mpicc "${flags[@]}" -DPARALLEL -I"$pluto_src/Parallel" "${sources[@]}" -lm -o "$test_out/mpi"
 for mode in 0 1 2; do mpirun -np 2 "$test_out/mpi" "$mode" 2; done
fi
printf 'Artifacts: %s\n' "$test_out"
