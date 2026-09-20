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
run_case() {
 local executable=$1 mode=$2 target=$3 failsafe=$4 status=0
 shift 4
 "$@" "$executable" "$mode" "$target" > "$test_out/case.log" 2>&1 || status=$?
 if [[ "$failsafe" == 0 && ( "$mode" == 2 || "$mode" -ge 4 ) ]]; then
   test "$status" -eq 1
   rg -q 'FOFC and primitive recovery failed with FAILSAFE=NO' "$test_out/case.log"
 else
   test "$status" -eq 0
   rg -q 'PASS local stage' "$test_out/case.log"
 fi
 printf 'PASS failsafe=%s mode=%s target=%s (%s)\n' "$failsafe" "$mode" "$target" "$executable"
}
for failsafe in 0 1; do
for rk in RK2 RK3; do
cc -DFAILSAFE="$failsafe" -DTIME_STEPPING="$rk" "${flags[@]}" -fsanitize=address,undefined -fno-omit-frame-pointer "${sources[@]}" -lm -o "$test_out/serial_${failsafe}_${rk}"
targets=(1 2)
if [[ "$rk" == RK3 ]]; then targets+=(3); fi
for target in "${targets[@]}"; do for mode in 0 1 2 3 4 5 6 7; do
 ASAN_OPTIONS=detect_leaks=0 run_case "$test_out/serial_${failsafe}_${rk}" "$mode" "$target" "$failsafe"
done; done
done
if [[ "${RUN_MPI:-0}" == 1 ]] && command -v mpicc >/dev/null; then
 mpicc -DFAILSAFE="$failsafe" "${flags[@]}" -DPARALLEL -I"$pluto_src/Parallel" "${sources[@]}" -lm -o "$test_out/mpi_${failsafe}"
 for mode in 0 1 2 3 4 5 6 7; do run_case "$test_out/mpi_${failsafe}" "$mode" 2 "$failsafe" mpirun -np 2; done
fi
done
printf 'Artifacts: %s\n' "$test_out"
