# RTX 5060 CUDA feasibility results

## Outcome

CUDA is worthwhile for this workload, but the data favors a thin native-real
FP64 cuDSS backend over a purpose-built iterative CUDA backend. On the exported
motor system, cuDSS reduced the complete 11-system cold nonlinear linear-solver
sequence from 7.634 s of six-thread CPU PCG time to 0.683 s, including a fresh
0.528 s symbolic analysis. The six-system continued sequence fell from 3.698 s
to 0.606 s, again including fresh symbolic analysis. Once a topology has been
analyzed, a typical value/RHS update costs about 0.8--1.0 ms to upload,
10.3--11.5 ms to refactor, and 1.65--1.68 ms to solve.

The recommended next step is an experimental `CudaLinearSystemBackend<double>`
implemented around cuDSS analysis/refactorization/solve. Retain the current CPU
assembly and full symmetric packed matrix, keep each nonlinear sequence on the
GPU, and cache cuDSS symbolic state by a topology identifier. Jacobi-PCG is a
useful fallback and validation path. Do not pursue GPU SpSV SSOR for this
matrix.

## Test system and methodology

- GPU: NVIDIA GeForce RTX 5060 Laptop GPU, compute capability 12.0, 8,080,064,512
  reported bytes of VRAM, direct PCIe attachment.
- Software: CUDA Toolkit 13.3.73, native `sm_120`, GCC/G++ 13 host compiler,
  and system cuDSS 0.8 CUDA-13 packages.
- Matrix: 150,405 real unknowns and 1,075,885 full symmetric CSR nonzeros.
- Tolerance: the exported XFEMM tolerance of `1e-8`; all iterative GPU results
  use the same preconditioned stopping metric.
- CPU oracle: the existing six-thread block-SSOR XFEMM backend. The export holds
  the incoming warm start and final CPU nodal solution for every PCG call.
- Cold point: 11 nonlinear systems, 7,046 CPU PCG iterations, 7.634 s measured
  CPU PCG wall time, and 9.115 s complete analysis time in this run.
- Continued 1-degree point: six nonlinear systems, 3,507 CPU PCG iterations,
  3.698 s CPU PCG wall time, and 4.558 s complete analysis time in this run.

The 11 cold matrices share topology hash `854ec6a1f2aa10bf`; the six continued
matrices share `6e7c172920fd6eeb`. The hash changes between the 0-degree and
1-degree points despite mesh reuse. Thus symbolic reuse works throughout each
Newton sequence today, but not across these two rotor angles under the current
ordering/coupling representation.

Timings below are GPU event timings except the explicitly labeled batch wall
time. Transfers, setup, analysis, factorization/preconditioner preparation, and
solve are separated. All GPU vectors stay resident during PCG iterations.

## Single-system comparison

The first cold matrix is the hardest individual system.

| Method | Preconditioner | Setup / analysis (ms) | Upload (ms) | Numeric prep / factor (ms) | Solve (ms) | Iterations | Estimated VRAM | Relative residual | Error vs CPU |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| CPU PCG | block-SSOR, 6 threads | matrix packing already done | -- | -- | 1,734.96 | 1,617 | host | 1.46e-8 | oracle |
| cuSPARSE PCG | Jacobi | 2.49 / 0 | 2.49 | 0.32 | 506.54 | 3,730 | 20.92 MiB | 9.95e-9 | 2.02e-9 |
| cuDSS | direct FP64 SPD | 0.28 / 527.87 | 1.24 | 16.08 | 3.14 | -- | 92.24 MiB | 2.90e-13 | 1.63e-9 |
| cuSPARSE PCG | SpSV SSOR-like | -- / 32.23 | 2.29 | 0.32 | 5,261.13 | 1,202 | 72.11 MiB | 1.68e-8 | 4.67e-9 |

The direct factor has 4,095,776 reported nonzeros. The input lower triangle has
613,145 stored entries, so factor fill is about 6.68x the stored triangle
(3.81x the full CSR NNZ). cuDSS's memory estimate remains small relative to the
available VRAM.

Jacobi needs 2.31x as many iterations as CPU block-SSOR on the first system but
is still 3.42x faster for the solve itself. The SSOR-like GPU preconditioner
reduces the iteration count by 26%, yet each iteration is dominated by two
serial sparse triangular solves and the result is 3.0x slower than CPU PCG and
10.4x slower than GPU Jacobi. Structural analysis reuse does not rescue it.

## Same-topology value/RHS changes

| Sequence / method | Systems | One-time analysis (ms) | Upload total (ms) | Numeric factor/prep total (ms) | Solve total (ms) | Measured linear-solver total (ms) | Speedup vs CPU PCG |
|---|---:|---:|---:|---:|---:|---:|---:|
| Cold / CPU block-SSOR PCG | 11 | -- | -- | -- | 7,634.48 | 7,634.48 | 1.00x |
| Cold / Jacobi PCG | 11 | 0 | 17.16 | 1.36 | 2,311.87 | 2,332.88 including setup | 3.27x |
| Cold / cuDSS | 11 | 527.87 | 9.86 | 125.18 | 19.81 | 682.71 including setup | 11.18x |
| Cold / SpSV SSOR PCG | 11 | 32.23 | about 17 | about 4 | 21,086.66 | about 21,140 | 0.36x |
| Continued / CPU block-SSOR PCG | 6 | -- | -- | -- | 3,697.63 | 3,697.63 | 1.00x |
| Continued / Jacobi PCG | 6 | 0 | 9.47 | 3.45 | 1,218.68 | 1,234.81 including setup | 2.99x |
| Continued / cuDSS | 6 | 519.33 | 5.79 | 69.92 | 11.08 | 606.12 including setup | 6.10x |

Without the one-time symbolic analysis, the cold cuDSS numeric/upload/solve work
is 154.84 ms (49.3x below CPU PCG), and the continued sequence is 86.79 ms
(42.6x below CPU PCG). This makes retaining node ordering and symbolic state a
high-value follow-up even though it was not required to demonstrate a win.

Replacing CPU PCG by the measured cuDSS sequence, without overlapping any CPU
work, projects the locally measured complete analyses from 9.115 s to about
2.164 s cold and from 4.558 s to about 1.466 s continued. Applied to the
original supplied six-core continued profile (2.589 s total, 1.768 s PCG), the
same 0.606 s GPU sequence projects about 1.427 s total, a 1.81x analysis-latency
improvement. If symbolic analysis can be retained across that operating-point
boundary, the projection is about 0.908 s, or 2.85x.

These complete-analysis numbers are projections, not integrated-backend wall
measurements. CPU assembly and GPU work can overlap across independent sweep
points, so they are deliberately not used as a throughput ceiling.

## Batch/concurrency throughput

Inputs were already resident and symbolic analysis was completed before the
timed direct-solver batch. Jacobi contexts share one CSR row/column structure;
each context has independent FP64 values and PCG vectors. The current cuDSS
test uses independent contexts and therefore conservatively duplicates each
context's structure and factor storage.

| Batch | Jacobi PCG systems/s | Jacobi avg solve latency (ms) | Jacobi estimated VRAM | cuDSS systems/s | cuDSS avg factor (ms) | cuDSS avg solve (ms) | cuDSS estimated VRAM |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.96 | 506.42 | 20.9 MiB | 53.83 | 10.32 | 1.51 | 92.2 MiB |
| 2 | 2.32 | 855.00 | 37.2 MiB | 72.46 | 18.64 | 2.27 | 184.5 MiB |
| 4 | 2.34 | 1,700.18 | 69.6 MiB | 80.44 | 37.85 | 4.65 | 369.0 MiB |
| 8 | 2.31 | 1,818.95 | 134.6 MiB | 88.92 | 72.01 | 10.53 | 738.0 MiB |
| 16 | 2.32 | 1,705.11 | 264.5 MiB | 91.21 | 96.77 | 15.70 | 1.44 GiB |
| 32 | 2.28 | 1,818.21 | 524.4 MiB | 91.63 | 108.51 | 26.55 | 2.88 GiB |

Jacobi throughput saturates by batch 2--4, consistent with sparse bandwidth and
kernel/reduction overhead being saturated by only a few systems. cuDSS scales
well through batch 8 and reaches practical saturation around 16 systems. Batch
8 is the best latency/throughput compromise; batch 16 is appropriate when
maximum sweep throughput matters. Batch 32 adds almost no throughput and uses
twice the memory.

The shared-index Jacobi implementation demonstrates that one resident CSR
structure can serve same-topology concurrent solves. A production cuDSS backend
should next compare independent contexts with cuDSS's batch CSR interface and
test whether shared row/column pointers reduce index storage without impairing
factorization concurrency.

## Correctness

- Five identical Jacobi-PCG solves produced the same 3,730 iterations, residual,
  nodal error, and bitwise FP64 solution fingerprint
  (`b01069056083f026`).
- Default-mode cuDSS is not bitwise deterministic: repeated parallel
  factorizations change a few low bits. Three repeats in one retained context
  differed from the first solution by at most 2.43e-12 in relative norm while
  retaining about 2.9e-13 residual. cuDSS 0.8 deterministic mode is evaluated
  separately in the AGE follow-up report.
- Every Jacobi-PCG and SpSV-PCG run converged at the exported `1e-8` tolerance.
- Jacobi conventional relative residuals are 3.29e-9 to 1.24e-8 across the 17
  matrices. Relative nodal errors versus XFEMM CPU are at most 1.99e-7.
- cuDSS conventional relative residuals are 7.83e-14 to 2.90e-13. Relative
  nodal errors versus XFEMM CPU are at most 2.05e-7. The larger solution
  difference than residual is consistent with the CPU oracle stopping at a
  much looser residual, not a direct-solver convergence failure.
- The source CPU 1-degree operating point retained its expected torque result
  (`-0.004371021721880083`). The AGE follow-up subsequently routes dumped GPU
  nodal fields through the unchanged XFEMM postprocessor and validates torque,
  flux linkage, and air-gap harmonics at three angles.
- The full serial repository suite passes all 51 enabled tests. Its analytical
  torque benchmark covers 0--90 degrees; the worst observed relative torque
  difference was 0.005553%, inside its 0.006% limit. The periodic air-gap torque
  variant produced the same reported values. Export is disabled throughout
  these tests, exercising the unchanged normal solver path.

No tolerance was loosened and all arithmetic is native real FP64.

The regression policy for an integrated backend should therefore require: the
existing CPU suite unchanged; exported-system residual and relative nodal-error
limits on every CUDA method; repeat-error limits rather than bitwise equality
for cuDSS; and end-to-end torque/force/flux comparisons on the motor sweep once
GPU solutions can flow through XFEMM postprocessing.

## Integration recommendation

1. Add an opt-in cuDSS implementation of `LinearSystemBackend<double>`; do not
   disturb Legacy/PETSc selection or normal non-CUDA builds.
2. Upload row/column structure once, perform `ANALYSIS` once per topology hash,
   and use `REFACTORIZATION` plus `SOLVE` for changing values/RHS.
3. Retain GPU allocations, factor storage, and previous solution across Newton
   iterations. Although a direct solve does not need a warm initial vector,
   keeping the result resident is useful for downstream/backend evolution.
4. Use the measured small angular union-bucket strategy from
   [AGE_TOPOLOGY_RESULTS.md](AGE_TOPOLOGY_RESULTS.md).  A complete-revolution
   union and an unbounded raw-topology cache are both inferior for this motor.
5. Schedule independent points in groups of roughly 8--16, with CPU assembly
   feeding resident GPU contexts asynchronously. Measure a true integrated
   motor sweep before productionizing the API.
6. Keep native-real FP64. Do not use the complex PETSc configuration as the
   performance baseline, and leave mixed precision for a separately validated
   experiment.

This evidence justifies integrating an existing library backend. It does not
justify writing a custom CUDA PCG/SSOR solver or a repository-wide CUDA
refactor.
