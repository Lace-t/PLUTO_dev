#!/usr/bin/env bash
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
src=$(cd "$here/../../../Src" && pwd)
out=$(mktemp -d /tmp/pluto-stage-exchange.XXXXXX)
mpicc -std=c17 -O1 -g -ffunction-sections -fdata-sections -Wl,--gc-sections \
 -Wl,--wrap=AL_Exchange -Wl,--wrap=MPI_Allreduce \
 -DPARALLEL -DPCP_STAGE_LIMITING=1 -DUPDATE_VECTOR_POTENTIAL=0 \
 -I"$here" -I"$src" -I"$src/RMHD" -I"$src/EOS/Ideal" -I"$src/Parallel" \
 "$here/test_stage_exchange.c" "$src/RMHD/pcp_stage.c" "$src/failsafe.c" \
 "$src/arrays.c" "$src/Parallel/al_exchange.c" -lm -o "$out/test"
printf 'Artifacts: %s\n' "$out"
mpirun -np 8 "$out/test"
