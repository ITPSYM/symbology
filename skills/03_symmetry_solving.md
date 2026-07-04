# Skill 03 — Symmetry invariant subspace solver

## Purpose

Given a target tensor and a discrete symmetry (cyclic, flip, parity),
compute the invariant subspace of the target's projection matrix. The
result is a square matrix `<target>_invariant.wxf` whose columns span
the symmetry-invariant subspace.

## CLI entry point

```bash
./bootstrap --solve-symmetry --symmetry <cyclic|flip|parity> --target <name>
```

## Flags

| Flag | Description |
|------|-------------|
| `--solve-symmetry` | Compute the invariant subspace. |
| `--symmetry <cyclic\|flip\|parity>` | Symmetry name. (`collinear` is excluded — it has its own solver, see [04_collinear_solving.md](04_collinear_solving.md).) |
| `--target <SEW_FpL \| FEC_W \| LEC_W>` | Target name. |

## How it works

1. **Locate the projection** for the target under
   `output/<symmetry>/<projection_filename(target)>` (e.g.
   `output/cyclic/SEW_5p1.wxf`).
2. If the projection file does not exist, **auto-invoke the projection
   pipeline** (`--project --symmetry <sym> --target <target>`) to generate
   it.
3. Compute the **invariant subspace** of the projection matrix:
   - The projection `P` is square and idempotent (`P^2 = P`).
   - The invariant subspace is the column space of `P` (the eigenspace
     with eigenvalue 1).
   - Computed via sparse RREF of `P - I`, then taking the kernel.
4. Write the invariant subspace as
   `output/<symmetry>/<target>_invariant.wxf`.

## Inputs

- The projection matrix at `output/<symmetry>/<projection_filename>`.
  Auto-generated if missing.
- No direct `data/` reads (the projection pipeline handles those).

## Outputs (to `output/<symmetry>/`)

| File | Description |
|------|-------------|
| `<target>_invariant.wxf` | The invariant subspace matrix. For `SEW_5p1`, dims `2 × 2`. |
| `summary.txt` | Updated by the projection pipeline if it was auto-invoked. |

## Key files

- `solve_symmetry.hpp` — `run_symmetry_solver(base_path, target, sym, ...)`.
- `run_solve.sh` — driver script that loops over `cyclic flip parity`.

## Conventions

- **The collinear symmetry is excluded** from this solver because its
  projection is dimension-shrinking (not square). Use the collinear
  solver instead.
- **Auto-invocation**: if the projection file is missing, the solver
  transparently runs `--project` to generate it. The user does not need
  to run `--project` separately.
- **Path resolution**: `base_path/"output"/sym.name`. No path override
  flags (yet).

## Smoke test

```bash
./bootstrap --solve-symmetry --symmetry cyclic --target SEW_5p1
```

Or via the driver script:

```bash
./run_solve.sh SEW_5p1   # runs cyclic, flip, parity
```

## Verified status

- `SEW_5p1` cyclic/flip/parity: all three projections are identity
  `diag(1, 1)` for `E6`, so the invariant subspace is the full 2-dim
  space. The `_invariant.wxf` files are byte-identical across the three
  symmetries.

This is consistent with the `E6` SEW construction enforcing the
symmetries by design; lower-weight symmetries (e.g. cyclic on
`first_w1`) are non-trivial.
