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
