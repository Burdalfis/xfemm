# Persistent planar host assembly

This note records the Section-4 performance baseline and the lifetime analysis
used to optimize the zero-frequency planar motor path. It deliberately does
not apply to harmonic or axisymmetric assembly without a separate audit.

## ResidualOnly baseline

The baseline was measured on AC power at commit `1868a4d`, using the
Section-3 simulator at `b9c8326` and:

```sh
./bin/persistent_motor_session_benchmark \
  motor_compare/xfemm_benchmark_geprc_speedx2_1404_4600kv.fem \
  1.2 residual
```

The 12 same-bucket transient samples averaged 2.08333 magnetic Newton
iterations and 221.773 ms per evaluation. The first sample needed three
iterations; the remaining samples needed two. Mean timings in milliseconds
per evaluation were:

| Region | ms |
|---|---:|
| nonlinear B-H/material evaluation | 31.010 |
| volume-element matrix and distributed RHS | 110.984 |
| AGE matrix assembly | 2.267 |
| explicit nodal RHS | 1.168 |
| boundary/periodic constraints | 15.962 |
| complete matrix assembly excluding material | 133.140 |
| sparse pack/value scatter | 17.099 |
| H2D | 2.144 |
| cuDSS factorization | 20.254 |
| cuDSS solve | 3.093 |
| D2H | 0.516 |
| direct linkage projection | 0.278 |
| residual calculation | 5.674 |
| nonlinear bookkeeping | 1.438 |
| result packaging | 1.193 |
| other measured/unaccounted | 5.857 |

The initialized cold solve took 2761.40 ms for 11 Newton iterations, including
71.26 ms of bucket-definition/resource work and 512.93 ms of symbolic
analysis. The first hot two-iteration solve took 219.74 ms.

The benchmark's `+0.25 degree` transition does not leave a 1.2-degree bucket.
A separate 0.3-degree-bucket control was therefore used to measure bucket
lifetime behavior without changing the production hot baseline. A real new
bucket took 862.12 ms for three Newton iterations, including 26.16 ms of
definition work, 3.26 ms of device-resource work, and 493.27 ms of symbolic
analysis. Returning to the retained original bucket took 313.43 ms for three
iterations, with 0.150 ms switching cost and no definition, resource, or
symbolic work.

## Static/dynamic classification

`FSolver::Static2D` currently recomputes volume-element geometry and source
metadata inside every nonlinear Newton iteration. For a persistent sliding-AGE
motor the volume mesh never moves: rotor position changes the AGE interpolation
only. The element work is classified as follows.

### A. Immutable for the imported persistent mesh/model

- Element node indices, block/material index, block-label index, and edge
  boundary markers.
- Triangle `p`/`q` basis-gradient coefficients, signed area, centroid, and
  three edge lengths.
- The symmetric local geometric matrices `Mx`, `My`, and `Mxy`.
- Material lamination type/fill and whether the material has a B-H curve.
- Winding/circuit association and the element-area coefficient multiplying a
  circuit current-density source.
- Fixed block current density and conductivity association.
- Derivative-boundary local matrix/RHS coefficients.
- Permanent-magnet coercivity and orientation. A Lua magnetization-direction
  function depends on immutable element coordinates and is evaluated while
  preparing the cache; it must not be evaluated concurrently.
- Point-source node mapping, fixed-Dirichlet node/edge mapping, and ordinary
  periodic/antiperiodic node pairs.
- For the fixed sparse graph, every element's six upper-triangular matrix
  destinations and three RHS node destinations.

The cache must be invalidated after mesh import or node renumbering. Changes to
material, label, boundary, or circuit definitions require rebuilding it.

### B. AGE-bucket dependent

- AGE `InnerShift`/`OuterShift`, interpolating node weights, and their numeric
  matrix contributions.
- The bucket union sparsity and corresponding persistent cuDSS numeric-value
  destinations.

No ordinary volume-element geometry is bucket dependent.

### C. Nonlinear magnetic-state dependent

- Element flux density derived from the current nodal `A` iterate.
- B-H interpolation results, reluctivity/permeability, differential
  reluctivity, and the nonlinear tangent matrix `Mn`.
- The `Mn * A` Newton RHS contribution.
- Incremental/frozen-permeability data derived from a previous solution.

These values must be recomputed for every applicable Newton iteration. The
constitutive interpolation and arithmetic are unchanged by Section 4.

### D. Winding-current/source dependent

- Circuit `J` or voltage-gradient `dV` selected from the configured circuit
  constraint and the precomputed area/conductivity/source integrals.
- The resulting distributed winding-source RHS and `Jprev` bookkeeping.

The circuit-to-element association and geometric multiplier are immutable; the
current-dependent scalar is not.

### E. Otherwise dynamic

- Clearing matrix/RHS numerical values for a new Newton iteration.
- Applying fixed values and periodic/antiperiodic transformations to the newly
  assembled matrix/RHS.
- Numeric packing/upload, factorization, solve, residual verification, and
  nonlinear relaxation/bookkeeping.
- The solution/warm-start vector and requested convergence tolerance.

The optimized path must preserve the exact boundary transformation order and
must not assume that a numerically zero entry is absent from the persistent
sparsity graph.

## Serial precompute/direct-entry result

The optimized single-thread path caches immutable triangle geometry, local
geometric matrices, fixed boundary/source terms, circuit association, and the
evaluated permanent-magnet direction. Each element then retains the exact six
upper-triangular host sparse-entry addresses for its active bucket. Subsequent
Newton iterations update those entries directly, and packing refreshes numeric
values without walking the linked sparse topology again.

Using the same command and 12-sample workload as the baseline, the hot
ResidualOnly evaluation fell from 221.773 ms to 149.009 ms (1.488x). The first
transient sample required three magnetic Newton iterations and the remaining
samples required two, unchanged from the baseline.

| Region | baseline ms | optimized ms | speedup |
|---|---:|---:|---:|
| nonlinear B-H/material evaluation | 31.010 | 31.064 | 1.00x |
| volume-element matrix and distributed RHS | 110.984 | 42.077 | 2.64x |
| complete matrix assembly excluding material | 133.140 | 63.826 | 2.09x |
| sparse pack/value scatter | 17.099 | 13.516 | 1.27x |
| complete ResidualOnly evaluation | 221.773 | 149.009 | 1.49x |

cuDSS was deliberately held constant: factorization was 20.240 ms versus
20.254 ms at baseline, and solve was 3.092 ms versus 3.093 ms. The persistent
CPU/cuDSS regression retained a 1.34e-9 nodal relative difference between the
independent CPU and GPU solvers, maximum linkage error 3.69e-14 Wb-turn, and
torque error 1.27e-11 N m. Tangent validation retained 2.58e-11 H maximum
absolute centered-FD error and 2.87e-6 maximum relative error. The compact
sparse unit regression also compares direct six-entry element assembly against
the scalar reference before and after a numerical wipe.

## Threaded element/material assembly

Four accumulation designs were considered for the fixed graph:

- A private full CSR plus RHS per thread is simple but would replicate roughly
  6.1 MB per thread for this mesh and require a bandwidth-heavy whole-matrix
  reduction after every Newton iteration.
- Element coloring permits concurrent direct scatter, but changes the global
  accumulation order and adds a parallel launch for every color.
- Row ownership avoids matrix races but either duplicates cross-partition
  element arithmetic or requires a second contribution exchange.
- Retaining one compact local contribution per element permits both material
  interpolation and local matrix/RHS calculation to run independently, then
  performs a serial element-order scatter. This has bounded storage, no
  atomics, and preserves the legacy global summation order exactly.

The fourth strategy was selected. Material/B-H interpolation is now a distinct
parallel pass, so the material timing no longer includes two clock calls per
element. Local matrix/RHS values are computed in a second static-scheduled
parallel pass, and their retained values are accumulated in original element
order. `XFEMM_ASSEMBLY_THREADS` overrides only assembly; the older
`XFEMM_NUM_THREADS` is accepted as a fallback. With neither set, OpenMP uses up
to 16 threads, the measured optimum on the target 16-logical-CPU machine.
Concurrent multi-process workloads can explicitly divide the CPUs with
`XFEMM_ASSEMBLY_THREADS`.

The threaded algorithm has a one-thread cost of 162.484 ms because retaining
per-element values adds a second memory pass. Thread speedups below therefore
use that like-for-like algorithmic baseline; the separate 149.009 ms result
above remains the serial-algorithm baseline.

| threads | material ms | complete matrix ms | volume element ms | evaluation ms | program CPU | speedup vs threaded 1T |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 35.336 | 73.726 | 50.574 | 162.484 | 99% | 1.000x |
| 2 | 23.877 | 61.867 | 39.101 | 139.565 | 115% | 1.164x |
| 4 | 17.882 | 58.409 | 35.318 | 129.919 | 138% | 1.251x |
| 8 | 11.799 | 59.341 | 36.273 | 124.835 | 164% | 1.302x |
| 16 | 11.416 | 51.302 | 28.339 | 116.572 | 190% | 1.394x |

The CPU percentage is for the complete benchmark including cold setup,
serial postprocessing, GPU waits, and bucket work; it is not the utilization
inside the two parallel loops. At 16 threads, sparse packing was 13.200 ms,
cuDSS factorization 20.371 ms, and cuDSS solve 3.108 ms. Those remain consistent
with the single-thread values (13.175, 20.258, and 3.093 ms respectively).

Against the original 221.773 ms Section-3 baseline, the selected host path is
1.902x faster. The algorithmic precompute/direct-entry change accounts for
1.488x on its own; threading improves the new threaded representation by
1.394x, or the serial optimized path by 1.278x. These factors are reported
separately because their baselines differ.

## Final breakdown and CUDA decision gate

A final selected-default repeat measured 116.591 ms/evaluation. Its hot
ResidualOnly breakdown was:

| Region | ms | share |
|---|---:|---:|
| nonlinear material/B-H | 10.734 | 9.21% |
| volume element matrix/RHS | 29.092 | 24.95% |
| AGE matrix | 1.996 | 1.71% |
| explicit RHS | 1.171 | 1.00% |
| boundary constraints | 15.256 | 13.09% |
| other matrix-construction work | 4.673 | 4.01% |
| sparse pack/scatter | 13.217 | 11.34% |
| H2D / cuDSS factor / solve / D2H | 2.126 / 20.337 / 3.099 / 0.509 | 22.36% total |
| residual / nonlinear bookkeeping / packaging | 5.785 / 1.451 / 1.197 | 7.23% total |
| direct linkage | 0.282 | 0.24% |
| unaccounted | 5.588 | 4.79% |

The remaining material plus volume-element work is 39.826 ms, or 34.16% of
the evaluation. Eliminating only that work has an ideal Amdahl ceiling of
1.519x additional speedup. A persistent device assembly path could also avoid
most of the 13.217 ms host pack/scatter and 2.126 ms H2D stages; eliminating
all three has a less realistic upper ceiling of 1.90x. During a 30%-duty live
transient the GPU sampled 15-28% SM and 1-2% memory utilization, so device
headroom remains.

GPU assembly is therefore still worth a bounded Section-5 investigation, but
not an unmeasured rewrite. The simplified host path gives the implementation
shape: upload immutable element nodes, six CSR destinations, geometric
matrices, fixed sources, material/circuit IDs, and B-H tables once; keep the
nodal iterate and CSR values resident; clear values on device; evaluate
material and six local entries per element; and reduce directly into the
cuDSS numeric values. A deterministic segmented/color reduction should be
benchmarked against FP64 atomics. AGE and boundary transformations can remain
host-driven initially only if their contributions can be applied without a
full matrix round trip. cuDSS must consume the resulting resident values
directly. Section 4 does not implement this path.

The CPU/cuDSS persistent test, tangent centered-FD validation, compact sparse
tests, and live motor transient/cache/bridge tests pass with the selected
thread count. ASan plus UBSan unit/analytical runs also pass; LeakSanitizer is
not usable under the application ptrace environment and was disabled for the
ASan run.
