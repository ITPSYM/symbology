# Skill 04 — Collinear constraint solver

## Purpose

Solve the collinear constraint `c · A = boundary`, where:

- `A` is the expanded SEW collinear basis (rank `2L+1`:
  `{sew_dim, 11, ..., 11}`).
- `boundary` is the expected collinear limit at loop order `L` (rank
  `2L`: `{11, ..., 11}`).
- `c` is the unknown coefficient vector (length `sew_dim`).

The constraint is enforced in the projected letter subspace (each
11-dim letter slot projected via the `--letter-projection` matrix, by
default `colprojdiv_w1` which projects to the 2-dim divergent
subspace). Pass `--letter-projection identity` to solve in the full
letter space (no projection). The finite part is
`R* = c·A - boundary`, which must be divergent-free when a divergent
projection is used (verified by the indicator-vector method).

## CLI entry point

```bash
./bootstrap --solve-collinear --target <SEW_FpL> --rhs <rhs.wxf|0> \
    --projection <finite|divergent> --letter-projection <file|identity> \
    [--basis <basis.wxf> ...] [--data-dir <dir>] [--output-dir <dir>]
```

## Flags

| Flag | Description |
|------|-------------|
| `--solve-collinear` | Run the collinear solver. |
| `--target <SEW_FpL>` | Target SEW name (e.g. `SEW_3p1` for 2-loop, `SEW_5p1` for 3-loop). |
| `--rhs <rhs.wxf>` or `--rhs 0` | RHS path. Required — exits with code 1 if missing. `--rhs 0` means an all-zero RHS constructed in-memory. |
| `--projection <finite\|divergent>` | Which projection to apply. Required — no default. |
| `--letter-projection <file\|identity>` | Letter-slot projection matrix. Required. A path (e.g. `output/collinear/colprojdiv_w1.wxf`) projects each 11-dim letter slot to a lower-dim subspace; `identity` skips projection (solve in full letter space). Relative paths resolve against the executable directory. |
| `--basis <basis.wxf>` | Expansion basis file (repeatable; highest weight first). If omitted, auto-detected as `first_w{N}_basis.wxf` for weights `target_weight-1` down to 2. |
| `--data-dir <dir>` | Data directory with seed files (default: `<exec_dir>/data`). Resolved against the executable directory. |
| `--output-dir <dir>` | Output directory (default: `<exec_dir>/output`). Same resolution as `--data-dir`. |

## How it works

1. **Collinear projection chain** (`solve_collinear.hpp`):
   - Copy `data/colprojdiv.wxf` → `output/collinear/colprojdiv_w1.wxf`
     (and `colprojfin_w1.wxf`).
   - Recursively compute `colprojdiv_w{N}.wxf` and `colprojfin_w{N}.wxf`
     for `N = 2..F` by contracting the weight-`N-1` projection with the
     FEC weight-`N` basis.
   - Compute the SEW-level projection `colprojdiv_<sew_name>.wxf` (e.g.
     `colprojdiv_SEW_5p1.wxf`).

2. **Expand the SEW collinear basis** (`tensor_expand.hpp`):
   - Load `output/collinear/SEW_<name>_basis.wxf` (shape
     `{sew_dim, FEC_F_basis, 11}`).
   - Contract axis 1 with each `first_w{N}_basis.wxf` (highest weight
     first, down to weight 2).
   - Result: `A` with rank `2L+1`, dims `{sew_dim, 11, ..., 11}`.

3. **Project to the divergent subspace** (`apply_colprojdiv_slots`):
   - Apply the `--letter-projection` matrix to each 11-dim letter slot of
     both `A` and the boundary. The default example is `colprojdiv_w1`,
     which projects to the 2-dim divergent subspace.
   - `A_proj` becomes `{sew_dim, 2, ..., 2}`; `boundary_proj` becomes
     `{2, ..., 2}`.
   - If `--letter-projection identity`, this step is skipped — the solve
     happens in the full 11-dim letter space.
   - Projection to the divergent subspace is **required** at `L ≥ 3`
     when `E1` has divergent-letter entries (`E1[0,0] = -2`,
     `E1[1,1] = -2`), so the boundary has divergent components and
     solving in the full space fails.

4. **Match positions** (`A_match`, `b_match`):
   - Iterate `A_proj` directly (do **not** collapse into
     `std::map<key, T>` — that overwrites duplicate sew entries and
     breaks the multi-unknown system).
   - Preserve the sew axis in `A_match` so positions with both sew
     entries contribute two coefficients: `c0·A[0,key] + c1·A[1,key] = b[key]`.
   - `b_match` has one entry per matching key.

5. **Linear solve** (`linear_solve.hpp`):
   - Reshape `A_match` to `(n_unknowns, n_constraints)` and `b_match` to
     `(1, n_constraints)`.
   - If too many constraints, sample `3 × n_unknowns` constraints, solve,
     then verify against all.
   - Uses modular RREF over `Z / 2^61` with reconstruction to `rat_t`.

6. **Indicator-vector verification** (skipped if `--letter-projection identity`):
   - Compute `R* = c·A - boundary`.
   - Collect the distinct letter indices appearing in `R*`'s non-zero
     entries.
   - Build an 11-dim indicator vector (1 at those indices, 0 elsewhere).
   - Project with the `--letter-projection` matrix. If the result is
     zero, `R*` is divergent-free (no entries at letters `{0, 1}`), so
     the collinear constraint is satisfied.

## Inputs

- `data/colprojdiv.wxf`, `data/colprojfin.wxf` (weight-1 seeds).
- `output/collinear/SEW_<name>_basis.wxf` (from `--project --symmetry collinear`).
- `output/collinear/first_w{N}_basis.wxf` (from the projection chain).
- The RHS file (passed via `--rhs`).

## Outputs

### To `output/collinear/`

| File | Description |
|------|-------------|
| `colprojfin_w1.wxf`, `colprojdiv_w1.wxf` | Copies of the `data/` seeds (skipped if present) |
| `colprojfin_w{N}.wxf`, `colprojdiv_w{N}.wxf` (`N ≥ 2`) | FEC-level projections |
| `colprojfin_<sew_name>.wxf`, `colprojdiv_<sew_name>.wxf` | SEW-level projections (e.g. `colprojdiv_SEW_5p1.wxf`) |

Empty projections (0 rows) are not written.

### To `output/<L>loop/` (SEW targets only)

| File | Description |
|------|-------------|
| `solMHV_LL.wxf` | Solution coefficient vector `(1 × n_unknowns)`, written only if the system is consistent. `L = target_weight / 2`. |

## Key files

- `solve_collinear.hpp` — `run_collinear_proj_chain`,
  `run_collinear_solver`, `apply_colprojdiv_slots`, `detect_chain_base_paths`.
- `linear_solve.hpp` — `solve_linear_system` (modular RREF + reconstruct).
- `tensor_expand.hpp` — `expand_tensor` (contracts with a basis chain).

## Conventions

- **`--rhs` is required.** Missing `--rhs` → exit code 1 with a helpful
  message. `--rhs 0` means an empty (all-zero) RHS tensor.
- **`--projection` is required** — there is no default. Pass either
  `finite` or `divergent` explicitly.
- **`--letter-projection` is required** — there is no default. Pass
  either a file path (e.g. `output/collinear/colprojdiv_w1.wxf`) or
  the literal `identity` to skip projection (solve in full letter space).
- **`--data-dir` / `--output-dir`**: default to `<exec_dir>/data` and
  `<exec_dir>/output`. Relative paths resolve against the executable
  directory (same convention as `bootstrap`).
- **SEW-level naming**: `colprojdiv_<sew_name>.wxf` where `sew_name`
  already includes the `SEW_` prefix (e.g.
  `colprojdiv_SEW_5p1.wxf`). The old name `colprojdiv_w6.wxf` was
  renamed to this.
- **Slot order**: must remain in natural ascending weight order
  `[w1, w2, w3, w4]`. The `expansion_perm` function places new letters
  immediately after the FEC axis to maintain this; placing them at the
  end reverses to `[w1, w4, w3, w2]` and breaks verification.
- **`colprojdiv_w1` swaps** the two divergent letters: for `E6`,
  `colprojdiv[0,1] = colprojdiv[1,0] = 1`. Do not assume it is a
  projection onto the first coordinate.
- **Every run generates a log file** in `logs/` recording tensor
  dimensions and file CRC32.

## Smoke test

```bash
# L=2: solve for c such that c·A = E1²/2 (the boundary)
./bootstrap --solve-collinear --target SEW_3p1 --rhs output/2loop/boundary_2L.wxf \
    --projection divergent --letter-projection output/collinear/colprojdiv_w1.wxf

# L=3: boundary = E1³/6 + E1·R2 (requires L=2 outputs to exist)
./bootstrap --solve-collinear --target SEW_5p1 --rhs output/3loop/boundary_3L.wxf \
    --projection divergent --letter-projection output/collinear/colprojdiv_w1.wxf

# Identity case (no divergent projection; only works when E1 has no divergent entries
# or when solving at L=2 where the boundary already lies in the finite subspace)
./bootstrap --solve-collinear --target SEW_3p1 --rhs output/2loop/boundary_2L.wxf \
    --projection divergent --letter-projection identity
```

For the full recursive workflow (computing the boundary too), use
`compute_rhs` instead — see [05_compute_rhs.md](05_compute_rhs.md).
`compute_rhs` internally invokes `--solve-collinear` for each loop order.

## Verified status

- **L=2 (`SEW_3p1`)**: unique solution `c[0] = 8`, 8 constraints
  verified. `R2` is divergent-free.
- **L=3 (`SEW_5p1`)**: unique solution `c[0] = -24, c[1] = 2`, all 32
  constraints verified. `R3` contains only letters
  `{2,3,4,5,6,7,8,9,10}` (no divergent letters `{0,1}`), so
  `R3 = R*` is divergent-free.

## Pitfalls

1. **Map overwrite bug**: never collapse `A` into `std::map<key, T>`
   when `sew_dim > 1`. Iterate `A` directly and preserve the sew axis.
2. **Divergent-subspace projection is required at `L ≥ 3`**: because
   `E1` has divergent-letter entries, the boundary has divergent
   components, and solving in the full space produces inconsistent
   systems (e.g. at L=3: 19 distinct ratio values).
3. **`expansion_perm` slot order**: new letters go immediately after
   the FEC axis, not at the end.
