# Sliding-air-gap cuDSS topology results

## Decision

Do not use one complete-revolution AGE superset as the default topology, and
do not cache all raw angle-specific topologies.  The cheapest measured design
is a persistent cuDSS context for a **small contiguous angular bucket** whose
CSR is the union of every AGE pattern possible in that bucket.  Move to the
next precomputed/cached bucket only when the rotor crosses its boundary.

For the representative motor, an exact continuous 0--1.2-degree bucket sampled
at every 0.4-degree graph transition added only 0.335% CSR NNZ.  It retained
one 97.1 MB cuDSS context instead of three contexts totaling 289.9 MB, and its
factor fill remained close to the angle-specific graphs.  A separate
two-position 0/1-degree union produced the same NNZ growth.  Across all 17
nonlinear matrices at those two operating
points it reduced cuDSS wall time from 1.789 s (two symbolic analyses) to
0.889 s (one analysis), while total numerical-factor time was slightly lower:
170.9 ms versus 179.5 ms.  This is the relevant result for a persistent
small-step transient evaluator.

The exact full-revolution union is valid and cuDSS does reuse its analysis when
retained entries change between zero and nonzero, but it is too expensive when
an angle needs many nonlinear factorizations.  It grows the full CSR by 149%,
factor fill by about 152%, and refactor time from about 10.7 ms to 77 ms.

## Experiment

The input was the same native-real FP64 motor system used for the CUDA
feasibility work:

- 150,405 unknowns and 1,075,885 full symmetric CSR entries;
- 900 nodes on each AGE ring;
- 8,100 undirected cross-ring entries at any one rotor position;
- RTX 5060 Laptop GPU (8 GB), CUDA 13.3, cuDSS 0.8;
- XFEMM tolerance `1e-8`, six-thread CPU result as the nodal oracle.

The diagnostic exporter retained the final nonlinear matrix by topology hash
for a mesh-reused, warm-started 0--125 degree sweep in 5-degree increments.
All completed angles had identical dimension and NNZ.  The benchmark inferred
the two AGE rings from the symmetric difference of two real CSR graphs.  The
complete-revolution topology adds the exact cross product of those rings;
inactive entries keep a numeric value of zero.  The observed/local topology is
the sorted union of the supplied real graphs.

The cuDSS objects, device CSR, values, RHS, solution, analysis state, and factor
allocations remain alive through each timed sequence.  On a topology hit, only
lower-triangle FP64 values and RHS are copied before factorization/solve.

## Topology orbit and recurrence

`xfemm-age-topology-enumerate` combines the actual CSR graphs with the actual
node coordinates from `.ans`.  It sorts both inferred rings geometrically,
extracts the real 0-degree cross-ring edge stencil, cyclically shifts that
stencil through every relative ring position, and hashes each resulting graph.

| Property | Result |
|---|---:|
| Nodes per ring | 900 |
| Cross-ring edges per angle | 8,100 |
| Distinct raw CSR patterns in `[0,360)` | 900 |
| Mechanical topology pitch | 0.4 degrees |
| 1-degree observed graph | exact cyclic shift 898 (equivalent to -2) |
| 5-degree observed graph | exact cyclic shift 888 (equivalent to -12) |
| Pattern recurrence | one mechanical revolution |

The 25 sampled positions from 0 through 120 degrees produced 25 distinct
hashes.  A separately solved 360-degree state reproduced the 0-degree physical
result, but floating-point angle positioning selected a different boundary
hash.  A production session must canonicalize angle modulo 360 degrees before
AGE positioning and topology lookup.

Caching effectiveness depends on the evaluation step:

- At 0.1-degree steps, each topology serves four evaluations.  A single active
  context avoids 75% of analyses even in the first revolution; a full cache
  approaches 100% reuse over repeated revolutions.
- At 1-degree or 5-degree steps, sampled positions do not repeat until the next
  revolution.  A small LRU cache gives essentially no first-revolution hits.
- Retaining 900 complete angle-specific contexts is not viable: each measured
  context needs about 96--97 MB, or roughly 87 GB for the full orbit.  Two
  retained contexts measured 193,590,504 bytes; the two-pattern union needed
  98,556,524 bytes.

For a monotonic fixed-step trajectory, the fraction of approximately 0.5 s
analyses avoided by an exact-topology cache is therefore:

| Mechanical step | First revolution | After `R` repeated revolutions |
|---|---:|---:|
| 0.1 degrees | 75% | `1 - 1/(4R)` |
| 1 degree | 0% | `1 - 1/R` |
| 5 degrees | 0% | `1 - 1/R` |

These fractions assume enough cache capacity to keep every topology needed for
later revolutions.  With a bounded small LRU, the 1- and 5-degree recurrence is
too distant to help.  A local union avoids analysis throughout its bucket with
one context and no such recurrence requirement.

This is why an angular union bucket is preferable to a large raw-topology
cache.  The cache remains useful for a small set of active motors, trial states,
or repeatedly visited buckets.

## Structural and cuDSS measurements

The following 25-angle table uses one final nonlinear system per sampled
angle.  Angle-specific contexts each pay analysis once.  Times are summed GPU
phase times; wall excludes the development-only host union builder.

| Metric | Angle-specific topology cache | Observed 0--120-degree sampled union | Complete-revolution union |
|---|---:|---:|---:|
| Full CSR NNZ | 1,075,885 each | 1,464,685 (+36.14%) | 2,679,685 (+149.07%) |
| Lower CSR NNZ | 613,145 each | 807,545 | 1,415,045 |
| Values + RHS update | 6,108,400 B | 7,663,600 B | 12,523,600 B |
| Symbolic analysis total | 12,799.3 ms | 627.5 ms | 772.5 ms |
| Numeric factor total | 268.7 ms | 1,659.8 ms | 1,933.9 ms |
| Solve total | 43.3 ms | 89.3 ms | 91.3 ms |
| Upload total | 31.2 ms | 28.9 ms | 43.9 ms |
| Factor NNZ | 4.067--4.116 million | 8,963,096 | 10,331,421 |
| cuDSS peak device estimate | about 96.2--96.9 MB each | 141,020,340 B | 161,285,276 B |
| Timed wall | 13,898.7 ms | 2,708.5 ms | 3,296.9 ms |

The observed union is exact only for the supplied 5-degree trajectory.  A
production bucket must include *all* intermediate 0.4-degree patterns within
its angular interval.  The topology enumerator provides the cyclic orbit needed
for a production bucket builder to generate those patterns without running a
field solve.

The narrow-bucket follow-up explicitly solved 0, 0.4, 0.8, and 1.2 degrees.
Those positions contained three distinct graphs (the boundary convention made
0 and 0.4 share one).  Their union results were:

| Continuous 0--1.2-degree bucket | Exact topology cache | Union bucket |
|---|---:|---:|
| Distinct cuDSS contexts | 3 | 1 |
| Full CSR NNZ | 1,075,885 each | 1,079,485 (+0.335%) |
| Factor NNZ | 4.070--4.096 million | 4,125,177 |
| Permanent/peak device estimate | 289,860,796 B total | 97,100,044 B |
| Symbolic analysis total | 3,206.9 ms | 1,066.3 ms |
| Numeric factors total | 40.7 ms | 29.1 ms |
| Solve total | 7.5 ms | 4.5 ms |
| Timed wall | 3,698.0 ms | 1,176.5 ms |

Absolute analysis times in this later run were about twice the earlier run due
to host load/power state; the within-run comparison and factor/fill behavior
are the relevant evidence.  Mesh reuse and warm starting reduced the small
angle updates to four Newton systems each, making this narrow union the closest
measured proxy for the intended transient hot path.

One final matrix per angle favors either union because it removes about 0.5 s
of analysis per new graph.  The nonlinear hot path changes the answer.  The
0--120-degree CPU sweep used 287 Newton linear systems (11 or 12 per angle).
Applying the measured per-factor costs projects roughly 16.7 s of cuDSS phase
work for angle-specific analyses, about 20.0 s for the sampled union, and about
24.5 s for the complete union.  These are projections, not integrated-backend
wall times, but the directly measured two-angle, 17-matrix sequence confirms
the crossover:

| 0/1-degree nonlinear sequence (17 matrices) | Angle-specific | Local union | Full-revolution union |
|---|---:|---:|---:|
| CSR growth | 0 | +0.335% | +149.07% |
| Analysis | 1,017.0 ms | 513.2 ms | 775.9 ms |
| Numeric factors | 179.5 ms | 170.9 ms | 1,330.0 ms |
| Solves | 27.6 ms | 25.6 ms | 69.4 ms |
| Uploads | 18.1 ms | 16.4 ms | 30.8 ms |
| Wall | 1,789.0 ms | 889.5 ms | 2,525.6 ms |

The development union/padding build took 336 ms for this 17-system test.  It
belongs in persistent session/bucket preparation, not in `evaluate()`; with it
included once, the local union still beat two independent analyses.

The real CPU 0--120-degree sweep (25 analyzed positions, 287 Newton systems)
took 242.320 s end to end, or 9.693 s per position, including CPU assembly,
PCG, and output/postprocessing.  For the standalone final-matrix sweep, adding
the one-time host union/padding step to the retained cuDSS wall gives 3.257 s
for the observed union and 3.990 s for the complete-revolution union, versus
13.899 s for one newly analyzed exact context per sampled angle.  File/process
startup is deliberately excluded because the future session keeps these
resources resident.

An integrated CUDA end-to-end motor sweep is intentionally not claimed here:
that would require the backend this benchmark is meant to justify.  The 287-
matrix phase projections above and the directly measured 17-matrix sequence
bound the linear-solver choice without folding simulated timings into the CPU
assembly/postprocessing total.

## Numerical and physical validation

All arithmetic was real FP64 and no tolerance was changed.

- Complete and observed union matrices repeatedly transitioned retained AGE
  entries between exact zero and nonzero values under one cuDSS analysis.
- Conventional relative residuals over the 25-angle run were approximately
  `1.8e-13` to `2.3e-13`.
- Maximum nodal relative error versus the six-thread XFEMM PCG oracle was
  `2.35e-7`; union and angle-specific direct solutions agreed at the scale
  expected from the oracle's looser `1e-8` stopping tolerance.
- A 0 -> 5 -> 0 -> 0 numeric refactor cycle returned the repeated field within
  `1.9e-14` relative norm for the union and `4.1e-14` for a topology-cache hit.
- cuDSS default mode is not bitwise reproducible.  cuDSS 0.8 deterministic mode
  produced identical hashes in the same cycle, at about 23% more factor time,
  3.2x solve time, and 29% more memory in this small test.  Use it as an
  optional regression/debug mode, not the performance default.

GPU solutions at 0, 60, and 120 degrees were scaled by XFEMM's unchanged
planar output conversion and inserted into copies of the CPU `.ans` files.
XFEMM's normal postprocessor then computed the physical quantities.  Across
the three positions, maximum CPU/GPU absolute differences were:

| Angle | CPU torque | GPU-field torque | Absolute difference |
|---:|---:|---:|---:|
| 0 degrees | -0.0022865878628 | -0.0022865881770 | 3.14e-10 |
| 60 degrees | -0.0022860776206 | -0.0022860787749 | 1.15e-9 |
| 120 degrees | -0.0022865759877 | -0.0022865771118 | 1.12e-9 |

| Quantity | Maximum absolute difference |
|---|---:|
| Torque | 1.16e-9 |
| Flux linkage A/B/C | 9.50e-12 / 8.17e-12 / 6.92e-12 |
| Air-gap vector potential sample | 2.49e-10 |
| First-harmonic gap components | 6.59e-8 |
| Sixth-harmonic gap components | 8.10e-8 |

At 0 and 360 degrees the CPU torque differed by only `1.64e-11`, and flux
linkages by about `1e-14`, independently confirming physical closure.

The complete repository regression suite passes all 51 enabled tests after the
export/benchmark additions, including the AGE torque and AGE snapshot parity
tests.  Export remains disabled unless its environment variable is present.

## Persistent `MotorSession` design

The benchmark points to the following ownership and state model:

1. `MotorSession` owns the loaded problem, fixed mesh, node ordering, AGE ring
   topology/orbit, material state, and eventually element-to-CSR scatter maps.
2. A small `AngleBucket` owns stable CSR row/column arrays containing every AGE
   pattern in its interval, a topology ID, cuDSS handle/config/data and matrix
   wrappers, device values, RHS/solution buffers, workspace, permutation,
   symbolic state, and factors.
3. `evaluate(theta, currents)` canonicalizes `theta`, selects a bucket, updates
   only numeric values/RHS, refactorizes/solves through all Newton trials, and
   commits the converged nodal vector only on acceptance.  Trial and committed
   vectors stay distinct and resident.
4. The final tangent values and cuDSS numerical factorization remain valid
   after convergence.  The existing factor can then solve additional RHS
   columns for future differential-inductance work; cuDSS natively supports
   single and multiple RHS.  No sensitivity formulation is added here.
5. Bucket construction, CSR unioning, analysis, allocation, and process startup
   must stay outside the hot path.  CPU assembly can feed one bucket while the
   GPU solves another motor/trajectory, but a single motor's accepted/trial
   state remains serialized explicitly.
6. Start with a narrow bucket (roughly the next few 0.4-degree patterns) and
   benchmark wider *continuous* buckets in the integrated backend.  Keep a
   bounded LRU of buckets subject to a VRAM budget.  Fall back to an exact
   topology context if an unforeseen pattern is absent rather than mutating an
   analyzed CSR graph.

This preserves all requested long-lived state and does not prevent later
multiple-RHS reuse of the final converged factorization.  It also avoids
optimizing around serialization, process startup, or rebuilding CUDA objects,
none of which belongs in the eventual transient hot path.

## cuDSS API evidence

NVIDIA documents that analysis comprises reordering plus symbolic
factorization, and explicitly permits changing matrix values followed only by
(re)factorization and solve.  The `cudssData_t` object owns internal state,
including factors.  The library also supports multiple RHS, reports permanent
and peak memory estimates after analysis, and offers an optional deterministic
mode.  See the official [cuDSS function documentation](https://docs.nvidia.com/cuda/cudss/functions.html),
[data types and phase definitions](https://docs.nvidia.com/cuda/cudss/types.html),
and [general description](https://docs.nvidia.com/cuda/cudss/general.html).
