# XFEMM CUDA feasibility benchmark

This directory is intentionally independent of the XFEMM solver build. It
reads real FP64 magnetostatic systems exported by the legacy backend and does
not introduce CUDA into normal XFEMM binaries.

## Export format and use

Set `XFEMM_LINEAR_SYSTEM_EXPORT` to a filename prefix before running a real
magnetostatic analysis. Every PCG invocation writes
`PREFIX.NNNNNN.xfemm-system`. Export is disabled when the variable is absent or
empty.

For long sweeps, `XFEMM_LINEAR_SYSTEM_EXPORT_BY_TOPOLOGY=PREFIX` atomically
overwrites `PREFIX.HASH.xfemm-system` after each solve.  It therefore retains
the last nonlinear system seen for every full-CSR topology without producing
one large file per Newton iteration.  The hash covers row offsets and column
indices, not numeric values.

Version 1 is a little-endian binary stream with no structure padding:

1. magic `XFEMMLS\0` (8 bytes)
2. `uint32`: version, endian marker (`0x01020304`), scalar bytes, column bytes
3. `uint64`: dimension, full-CSR NNZ, process-wide solve index
4. `uint32`: flags (bit 0 warm start, bit 1 full symmetric CSR), CPU threads
5. `double`: XFEMM tolerance, SSOR relaxation
6. `int64`: CPU PCG iterations
7. `double`: CPU PCG preconditioned relative residual
8. arrays: `uint64 rowOffsets[n+1]`, `int32 columnIndices[nnz]`,
   `double values[nnz]`, `double rhs[n]`, `double initialSolution[n]`, and
   `double cpuSolution[n]`

The CPU solution is the validation oracle. The reported legacy PCG residual is
its preconditioned stopping metric; the tools also calculate the conventional
`||b-Ax||2/||b||2` residual.

Build the format inspector without CUDA:

```bash
cmake -S benchmarks/cuda -B build-cuda-inspect -DXFEMM_CUDA_BENCHMARK=OFF
cmake --build build-cuda-inspect
```

Build the GPU benchmark once a CUDA toolkit is installed:

```bash
cmake -S benchmarks/cuda -B build-cuda-benchmark \
  -DCMAKE_BUILD_TYPE=Release -DCUDSS_ROOT=/path/to/cudss
cmake --build build-cuda-benchmark -j
```

The RTX 5060 feasibility run used the system CUDA 13.3 and cuDSS packages. On
the tested Ubuntu installation, GCC 13 must be selected explicitly because the
system-default GCC is newer than the host compiler range supported by nvcc:

```bash
cmake -S benchmarks/cuda -B build-cuda-benchmark \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CUDA_COMPILER=/usr/local/cuda-13.3/bin/nvcc \
  -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++-13 \
  -DCUDAToolkit_ROOT=/usr/local/cuda-13.3 \
  -DCMAKE_CUDA_ARCHITECTURES=120
cmake --build build-cuda-benchmark -j
```

Run the native-real Jacobi-PCG and cuDSS sequence benchmark with one or more
same-topology files in Newton order:

```bash
build-cuda-benchmark/xfemm-cuda-benchmark --repetitions=3 \
  exported/system.000000.xfemm-system \
  exported/system.000001.xfemm-system
```

Select one method and measure resident independent solves at the requested
concurrency levels:

```bash
build-cuda-benchmark/xfemm-cuda-benchmark --method=pcg \
  --batch=1,2,4,8,16,32 exported/system.000000.xfemm-system
build-cuda-benchmark/xfemm-cuda-benchmark --method=cudss \
  --batch=1,2,4,8,16,32 exported/system.000000.xfemm-system
```

The separate SpSV executable measures the SSOR-like experiment without
complicating the main benchmark:

```bash
build-cuda-benchmark/xfemm-cuda-ssor-benchmark \
  exported/system.000000.xfemm-system
```

## Integrated persistent-session profile

`persistent_motor_session_benchmark.cpp` links the production
`PersistentMotorSession`; unlike the standalone feasibility tools above, it is
built from a cuDSS-enabled XFEMM tree. It reports session initialization,
same-bucket transient-like angle/current steps, a new-bucket symbolic event,
and a return to a cached bucket:

```bash
cmake --build build-cudss-session --target persistent_motor_session_benchmark
cfemm/bin/persistent_motor_session_benchmark \
  mfemm/testing/radial_machine/data/radial_machine_sliding.fem
```

Rows prefixed with `PROFILE` contain the complete per-evaluation breakdown;
`PROFILE_SUMMARY` reports the mean/min/max for the repeated hot path.

Sliding-air-gap graph experiments use the actual exported matrices:

```bash
# Complete cross-ring superset for an arbitrary mechanical revolution.
build-cuda-benchmark/xfemm-cuda-age-topology-benchmark \
  reference.xfemm-system shifted.xfemm-system system*.xfemm-system

# Union only the supplied trajectory/bucket patterns.
build-cuda-benchmark/xfemm-cuda-age-topology-benchmark --observed-only \
  reference.xfemm-system shifted.xfemm-system system*.xfemm-system

# Enumerate the exact cyclic topology orbit using real CSR and node geometry.
build-cuda-benchmark/xfemm-age-topology-enumerate \
  reference.xfemm-system shifted.xfemm-system reference.ans
```

The AGE benchmark keeps one cuDSS context per angle-specific topology for its
cache baseline and one context for the union topology.  Set
`XFEMM_CUDSS_DETERMINISTIC=1` to exercise cuDSS 0.8 deterministic mode.  The
optional `--dump-solutions PREFIX` output plus `xfemm-replace-ans-solution`
exists only to route GPU nodal fields through the unchanged XFEMM
postprocessor.  Planar magnetostatic internal unknowns require XFEMM's existing
`4*pi*1e-5` output scale when inserted into an `.ans` file.

The measurements and integration recommendation from the reference motor case
are in [RESULTS.md](RESULTS.md).
The sliding-air-gap topology/cache experiment and persistent-session design
recommendation are in [AGE_TOPOLOGY_RESULTS.md](AGE_TOPOLOGY_RESULTS.md).
