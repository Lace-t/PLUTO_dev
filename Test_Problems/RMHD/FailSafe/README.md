# RMHD FailSafe regression

With `PHYSICS=RMHD` and `FAILSAFE=YES`, the recovery order is:

1. Save the full-step fluid state and attempt the normal integration at `dt0`.
2. On failure, restore that state and retry at `dt0` with a first-order
   LF+CT fallback.
3. On further failures, restore the same state, retain the fallback regions,
   add any newly failed regions, and retry at `dt0/2`, `dt0/4`, etc.
4. On exhaustion, write per-rank bad-cell snapshots, flush logs, and abort.

`RMHD_FAILSAFE_MAX_RETRIES` controls the number of halvings, independently of
`PCP_STAGE_LIMITING`. It defaults to 4 and accepts 0 through 30. Zero still
allows the original-dt local fallback. An explicitly defined legacy
`PCP_STAGE_MAX_RETRIES` supplies the default if the new option is absent.
The former grid-scaled retry floor is no longer used; attempts stop early
if the proposed dt is nonfinite, nonpositive, or cannot advance floating-point
time. A successful smaller dt remains in `g_dt`; the existing `NextTimeStep`
CFL and growth rules govern later steps.

The failed cell receives `FLAG_HLL`. Its radius-one Chebyshev neighbourhood
receives `FLAG_FLAT`; the surrounding radius-two shell receives `FLAG_MINMOD`.
These remain useful diagnostics and mark the troubled region, but are not the
low-order operator: on every retry, `g_rmhd_failsafe_low_order` selects
piecewise-constant states with their shared staggered normal field,
`LF_Solver`, and piecewise-constant CT edge-EMF reconstruction. Thus the
fluid fluxes and the face magnetic fields are advanced from one compatible
first-order FV+CT construction. The mode is cleared before a successful
return. It is global for the retried full step because local CT limiting needs
owner-synchronized shared edge EMFs, which this failsafe does not cache.
Failure seeds are exchanged before dilation across MPI boundaries. Flags are
accumulated separately from the evolving data and restored before every retry.
No final cell state is spliced into an otherwise accepted update.

The saved state includes conservative and primitive arrays (including ghost
cells), initial flags, staggered magnetic fields, time/step/stage, the timeStep
structure and its pointed-to `cmax` array. Vector potential is saved when
enabled, and GLM speed is restored in GLM builds. Flux/EMF workspace is
recomputed by the normal integrator. Non-RMHD FailSafe keeps its existing path.
Local fallback requires `SHOCK_FLATTENING=MULTID`. Particle, forced-turbulence,
FARGO and HIGH_ORDER configurations are rejected at compilation because this
implementation does not back up their additional persistent state.

Stage-level failures must return collectively before any rank skips further
communication. Standard RK C2P and PCP checks already do this. A final error
reduction ensures all ranks also agree on the return from `Integrate`.
The final snapshot is named
`bad_cells_failsafe_final_step*_stage*_rank*.dat` and reports up to 1024 bad
active cells per rank; a failure confined to ghost cells may yield zero rows.

Run the serial regression with:

```sh
bash Test_Problems/RMHD/FailSafe/test_failsafe.sh
```

It compiles production `failsafe.c` and `arrays.c` against an integration stub
under AddressSanitizer and UndefinedBehaviorSanitizer. It exercises PCP off/on,
zero/two allowed halvings, immediate success, local fallback success, reduced
dt success, accumulated neighbourhood flags, full rollback, exhaustion with
and without C2P flags, diagnostics, and a dt unable to advance time. Temporary
binaries and logs are retained under `/tmp/pluto-rmhd-failsafe-test.*`.

For the two-rank test, from the PLUTO root:

```sh
mpicc -std=c17 -O1 -g -DPARALLEL \
  -ITest_Problems/RMHD/FailSafe -ISrc -ISrc/RMHD -ISrc/EOS/Ideal -ISrc/Parallel \
  Test_Problems/RMHD/FailSafe/test_failsafe.c Src/failsafe.c Src/arrays.c \
  -lm -o /tmp/test_rmhd_failsafe_mpi
cd /tmp
mpirun -np 2 ./test_rmhd_failsafe_mpi 3
mpirun -np 2 ./test_rmhd_failsafe_mpi -2
```

The first MPI run must succeed after two halvings. The second intentionally
exhausts retries without C2P flags and must exit nonzero after both ranks write
diagnostics. The test emulates ArrayLib's x-direction flag exchange using
real MPI messages; its integration stub reports errors only on rank zero.

Validation on 2026-09-09 also used fresh 2D CT RMHD shock-tube builds, with
PCP off/RK3 and PCP on/RK2. A test-only linker wrapper around `AdvanceStep`
injected failures into the first two attempts of a proposed final step.
Both runs tried `0.02`, `0.02` with local fallback, then accepted `0.01` and
continued to `tstop=0.02`, exercising the main-loop final-step correction.
These runs establish recovery control flow, not stability for arbitrary
production initial conditions.
# RK cell-average preflight (2026-09-16)

`bash test_stage_check.sh` exercises the real admissibility checker and the
RK pre-C2P gate with total/reduced energy and strict/compatibility modes,
under ASan/UBSan. The test checks valid input, negative density, negative q,
negative Psi, nonfinite momentum and a nonfinite CT face; U/V must remain
unchanged. `PCP_NR` and `PCP_STAGE_LIMITING` are not required.

`RMHD_FAILSAFE_PCP_CHECK=YES` enables the RK gate with failsafe.
`RMHD_FAILSAFE_PCP_STRICT=NO` is the default compatibility policy: failures
of the full check that pass finite/D/energy/q necessary conditions are
counted as `deferred_to_c2p` and left to native recovery. This is not full
PCP acceptance. Set STRICT=YES to reject every full-check failure.

The stage test also builds with `-DPARALLEL` using mpicc and `-ISrc/Parallel`;
link `pcp_stage.c`, `pcp_nr.c` and `arrays.c` with function/data sections and
`--gc-sections`. Run with two ranks: only rank zero injects invalid data,
but both ranks must return failure collectively.
