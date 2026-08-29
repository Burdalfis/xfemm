# Persistent planar motor tangent

For the supported zero-frequency planar stranded-series winding path, the
branch linkage at fixed geometry is a linear functional of the physical nodal
vector potential:

```
lambda_j = C_j A
```

`PersistentMotorSession` precomputes `C` from signed block turns, planar depth,
triangle area, and the linear nodal basis. This is the same expression as the
validated zero-current `FPProc::GetStrandedLinkage` calculation.

`FSolver::Static2D` solves for an internal variable `x` with
`A = (4 pi 1e-5) x`. For a private flat-current-density region created from an
original series circuit, differentiation of the assembled source term gives

```
d b / d i_j = (0.01 / depth) C_j^T.
```

The permanent-magnet, fixed boundary, and geometry/AGE terms do not vary with
circuit current. At a fixed rotor angle, the nonlinear magnetic Newton tangent
`K = dF/dx` therefore gives

```
K d x / d i = d b / d i
d lambda / d i = C (4 pi 1e-5) K^-1 (0.01 / depth) C^T.
```

The sign is positive because `Static2D` subtracts the element source into the
assembled RHS after its local `be` convention. Homogeneous Dirichlet and
periodic/antiperiodic transformations are applied to every sensitivity RHS in
the same order as the nonlinear solve.

The final cuDSS numeric factorization is reused. XFEMM's legacy nonlinear loop
assembles that tangent at the input to its final Newton update, then accepts the
updated field when the relative update is below tolerance. The extraction API
reports that final update norm so the small distinction is observable; it does
not silently assemble or refactorize another matrix.

Normal cuDSS mode solves all circuit columns as one dense multi-RHS solve.
cuDSS deterministic mode currently rejects changing a retained dense solve
from one RHS to multiple RHS, so that mode performs one retained-factor solve
per column. Neither path repeats magnetic assembly, symbolic analysis, or
numeric factorization.
