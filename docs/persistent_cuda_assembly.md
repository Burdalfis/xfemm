# Persistent CUDA planar assembly

## Scope

This document records the data-lifetime decision for optional CUDA numerical
assembly in the existing zero-frequency, planar, persistent-motor cuDSS path.
The legacy and persistent host assembly paths remain the correctness oracle and
fallback. Harmonic, axisymmetric, and disposable-solve paths are out of scope.

## Phase-0 ownership audit

The current `CudssContext` already owns stable device buffers for CSR row
offsets, column indices, numerical values, right-hand side, and solution. Its
cuDSS CSR descriptor is created with the device value pointer, and the handle is
bound to the context's nonblocking CUDA stream. Consequently:

1. cuDSS numerical refactorization can consume numerical CSR values already on
   the GPU. No cuDSS API requires those values to originate in host memory.
2. The analyzed matrix descriptor retains the same device CSR pointers for the
   complete lifetime of a context. Repeated numeric refactorizations do not
   recreate the graph.
3. The host path currently reconstructs an exact upper CSR matrix, transposes
   it to the lower form used by cuDSS, scatters into the bucket-union host value
   array, and copies the full values and RHS to the device on every solve.
4. CUDA assembly can bypass that copy by writing the context's existing
   `d_values` and `d_rhs` buffers before cuDSS factorization on the same stream.
5. One resident AGE bucket must retain its CSR graph, device values/RHS/solution,
   cuDSS matrix descriptors, handle/config/data objects, analysis and factor
   storage, CUDA stream, and bucket-specific element-to-CSR mappings. Immutable
   host mesh/material metadata may be shared or retained cheaply outside the
   device context.
6. A single stream can strictly order value/RHS clear, volume assembly, small
   AGE updates, boundary transformations, factorization, and solve. Per-phase
   CUDA events may be recorded without synchronizing each phase; normal
   production operation should synchronize only where host control needs a
   result.
7. Switching AGE buckets selects a different context, stream, graph, mappings,
   and device buffers. A cached return reuses them. Eviction destroys those
   device/cuDSS resources. The committed nodal magnetic state remains in the
   session-owned host solution, independently of context eviction, so a rebuilt
   bucket can warm-start from it.

The gate therefore passes: device-assembled values can reach cuDSS without a
device-to-host-to-device matrix round trip.

## Intended device lifetime

Immutable data uploaded once per compatible persistent mesh comprises element
nodes, compact geometric coefficients, material/source identifiers and data,
and winding associations. Bucket-specific mappings translate each local
symmetric element entry and each AGE/constraint contribution to the lower CSR
owned by that bucket. Mutable context storage comprises nodal A, CSR values,
RHS, optional per-element constitutive state, and circuit-current scalars.

The first solve of a new bucket may use host assembly to establish and analyze
the bucket-union graph. Subsequent hot solves may use the device plan. This
bootstrap does not weaken the hot-path architecture and keeps exact-topology
fallback behavior unchanged.

## Constraints and nonlinear control

Late Dirichlet, periodic, and antiperiodic operations are numerical matrix/RHS
transformations; they do not require a host matrix. Their static node and CSR
position maps can be precomputed. Constraint operations must execute in the
same sequence as the host solver, with each independent row update parallelized
but constraint pairs ordered.

The current nonlinear controller computes convergence and relaxation on the
host and therefore still needs the solved nodal vector after each iteration.
That transfer measured about 0.5 ms and is not an initial optimization target.
The production CUDA path should compute the linear residual on the device and
return only its scalar norm. An explicit parity/debug mode may download final
matrix/RHS values, but ordinary solves must not do so.

The retained tangent API continues to use the context's last numerical
factorization and multi-RHS cuDSS solve. Device assembly must neither replace
that context nor refactorize solely for tangent extraction.

## Implemented bounded path

`CudssSessionOptions::assemblyBackend` now selects `Host`, `CudaAtomic`, or
`CudaDeterministic`; `Host` remains the default. The CUDA modes are rejected by
construction outside the persistent real-valued planar magnetostatic path.
Each AGE bucket owns a `CudaPlanarAssembly` beside its existing `CudssContext`.
The bucket's first exact host solve establishes and analyzes the validated
union graph; hot evaluations then write the context's existing device CSR
values and RHS directly.

The mesh-lifetime plan contains triangle nodes, geometric `p`/`q` terms,
`Mx`/`My`/`Mxy`, fixed matrix and permanent-magnet RHS terms, material and
circuit IDs, exact element-to-CSR destinations, flattened B-H tables, explicit
nodal RHS, and the ordered boundary-transform schedule. Mutable device storage
contains the nodal iterate, circuit sources/cases, element reluctivity and
differential-reluctivity values, local matrix/RHS contributions, AGE additions,
and residual reductions. The original FP64 Hermite B-H interpolation and
lamination formulas are used without approximation.

The default CUDA scatter uses FP64 atomics. The deterministic alternative
precomputes contributor lists for every CSR entry and RHS node and reduces each
destination in a stable order. Both preserve the lower-triangular cuDSS graph.
AGE arithmetic remains on the host because it costs about 2 ms; only its small
changed-entry list is uploaded. Dirichlet, periodic, and antiperiodic
transformations execute directly on the device in the same ordered sequence as
the host operations. Persistent CUDA events delimit the clear/upload,
material, element, scatter, AGE, and constraint phases. The volume batch has
one completion synchronization rather than one per kernel, and the constraint
completion orders cuDSS factorization on the same stream.

The nonlinear controller still downloads the solved nodal vector (about
0.50 ms) for relaxation and host control. Matrix/RHS downloads occur only when
`XFEMM_CUDA_ASSEMBLY_PARITY=1`. The normal path evaluates the residual on the
device and transfers two reduction scalars. The converged factorization is the
same object consumed by the retained multi-RHS magnetic-tangent API; tangent
extraction adds no matrix transfer or refactorization.

## Validation and measured result

On the 150,405-node GEPRC production mesh, host-versus-CUDA pre-solve parity
over hot, new-bucket, and cached-return evaluations gave maximum absolute CSR
and RHS differences of `1.88e-13` and `2.48e-12`; matrix symmetry difference
was zero. Relative entry errors are dominated by values close to zero (up to
about `5e-8`) and are not physically significant at these absolute errors.
Direct B-H sampling gave zero reluctivity difference and
`2.27e-13` maximum differential-reluctivity difference. A short coupled
transient differed from the host path by at most `5.81e-11 A`,
`6.55e-16 Wb-turn`, and `1.45e-13 N-m`, with identical accepted times, PWM,
active-set, and bucket sequences.

The 12-sample hot ResidualOnly benchmark measured:

| path | mean evaluation | material kernel | element kernel | scatter | constraints | factor + solve |
|---|---:|---:|---:|---:|---:|---:|
| Host, 16 assembly threads | 118.25 ms | 11.68 ms host | 29.09 ms host reference | 13.29 ms pack | 15.26 ms host reference | 23.43 ms |
| CUDA atomic | 51.95 ms | 1.02 ms | 0.98 ms | 0.32 ms | 5.91 ms | 23.35 ms |
| CUDA deterministic | 57.44 ms | 1.06 ms | 0.99 ms | 0.28 ms | 5.92 ms | 23.41 ms |

The atomic path is 2.28x faster than the fresh 118.25-ms host measurement
(2.24x against the 116.591-ms Section-4 reference). It reports 284,183,272
bytes of persistent bucket device allocations. A capacity-two transient
stayed around 844 MiB in `nvidia-smi` while crossing and evicting buckets;
committed linkage survived an evict/rebuild/rollback sequence to
`3.13e-17 Wb-turn` and torque to `4.07e-15 N-m`.

The successful toolchain was `/usr/local/cuda-13.3/bin/nvcc`, CUDA compilation
tools 13.3.73, with `/usr/bin/g++-13` 13.4.0 as host compiler. The target is an
RTX 5060 Laptop (`sm_120`) on Linux driver 580.173.02. The required configure
setting on this machine is:

```text
-DCMAKE_CUDA_COMPILER=/usr/local/cuda-13.3/bin/nvcc
-DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++-13
-DCMAKE_CUDA_ARCHITECTURES=120-real
```

`120-real` emits native SASS/cubin and deliberately omits PTX. NVIDIA's CUDA
13.x minor-version compatibility floor is an R580 driver, but that mode
requires an explicit native architecture and does not support PTX produced by
a newer toolkit. CUDA 13.3's corresponding full-feature driver is
610.43.02 or newer. Thus this R580 installation can run the tested native
cubin, but a build containing CUDA-13.3 PTX (`120` or `120-virtual`) fails at
runtime with an unsupported PTX/toolchain error unless the driver is upgraded
to the CUDA-13.3 level (or a supported compatibility package is installed).
Future builds on R580 must retain `120-real`; the command above is the supported
configuration and the PTX failure mode is intentional rather than silently
falling back.

The compatibility rules are documented in NVIDIA's CUDA Compatibility Guide
(`docs.nvidia.com/deploy/cuda-compatibility/minor-version-compatibility.html`)
and CUDA 13.3 release notes. Compute Sanitizer in the installed toolkit reports
the RTX 5060 Laptop as unsupported, so `memcheck` and `racecheck` could launch
but could not instrument this device. All kernel launches are nevertheless
checked immediately, and host-facing ASan/UBSan regressions pass.
