#!/usr/bin/env bash
set -euo pipefail
test_src=$(cd "$(dirname "$0")" && pwd)
pluto_src=$(cd "$test_src/../../../Src" && pwd)
test_out=$(mktemp -d /tmp/pluto-rmhd-failsafe-test.XXXXXX)
cd "$test_out"
for passes in 1 4; do
for pcp in 0 1; do
 for retries in 0 2; do
  cc -std=c17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
    -I"$test_src" -I"$pluto_src" -I"$pluto_src/RMHD" -I"$pluto_src/EOS/Ideal" \
    -DRMHD_FAILSAFE_LOCAL_PASSES="$passes" -DPCP_STAGE_LIMITING="$pcp" -DRMHD_FAILSAFE_MAX_RETRIES="$retries" \
    "$test_src/test_failsafe.c" "$pluto_src/failsafe.c" "$pluto_src/arrays.c" \
    -lm -o "test_${pcp}_$retries"
  modes=(0 1 -1 -2 -3)
  if (( retries > 0 )); then modes+=(2 3 "$((passes+1))"); fi
  for mode in "${modes[@]}"; do
    case_out="$test_out/passes${passes}_pcp${pcp}_retries${retries}_mode${mode}"
    mkdir -p "$case_out"
    cd "$case_out"
    status=0
    log="passes${passes}_pcp${pcp}_retries${retries}_mode${mode}.log"
    ASAN_OPTIONS=detect_leaks=0 "$test_out/test_${pcp}_$retries" "$mode" > "$log" 2>&1 || status=$?
    expected=0
    if (( mode < 0 )); then expected=1; fi
    if (( status != expected )); then
      cat "$log"
      exit 1
    fi
    if (( mode < 0 )); then
      rg -q 'RMHD_FAILSAFE exhausted' "$log"
      # Diagnostic snapshots are disabled in this source configuration.
      # Validate their contents if present; exit/rollback checks always run.
      if (( mode == -1 )) && [[ -f bad_cells_failsafe_final_step12345678901_stage02_rank0000.dat ]]; then
        rg -q '# bad_cells=1 written_cells=1' bad_cells_failsafe_final_step12345678901_stage02_rank0000.dat
      else
        test ! -e bad_cells_failsafe_final_step12345678901_stage02_rank0000.dat
      fi
    else
      rg -q '^PASS' "$log"
    fi
    cd "$test_out"
  done
 done
done
done
printf 'PASS: PCP off/on, normal/local/reduced success, exhaustion, unflagged errors, invalid dt.\nArtifacts: %s\n' "$test_out"
