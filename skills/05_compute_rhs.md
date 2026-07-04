# Skill 05 — Recursive RHS (boundary) computation

## Purpose

Recursively compute the collinear boundary `E[L]`, the remainder `R[L]`,
and the boundary tensor `boundary_LL` from loop 2 up to a target loop
order `L`. This is a standalone executable that **couples with** the main
`bootstrap` module: it invokes `./bootstrap --extend`, `--sew`, and
`--project` as needed to generate missing SEW basis files.

The boundary at loop `L` is computed from the lower-loop results:

```
L=2: boundary = E1^2 / 2
L=3: boundary = E1^3 / 6 + E1 · R2
L=4: boundary = -E1^4 / 12 + E2^2 / 2 + E1 · R3
L=5: boundary = E1^5 / 20 - E1·E2^2 / 2 + E2·R3 + E1 · R4
```

where `E_L` is the expanded collinear projection of `hepMHV_LL`, and
`R_L = E_L - boundary_L` is the remainder (which must be divergent-free).

## CLI entry point

```bash
./compute_rhs --target <SEW_FpL> [--data-dir <dir>] [--output-dir <dir>]
```

## Flags

| Flag | Description |
|------|-------------|
| `--target <SEW_FpL>` | Target SEW name (required). Loop order `L = (F+L)/2`. Supported `L = 2..5`. |
| `--data-dir <dir>` | Data directory with seed files. Default: `<exec_dir>/data`. |
| `--output-dir <dir>` | Output directory. Default: `<exec_dir>/output`. |
| `-h` / `--help` | Print usage. |

## How it works (per loop order `L`)

1. **Load `E1`** from `data/E1.wxf` (rank 2, `11 × 11`, 5 nnz).
2. **Recursively compute `E[2..L-1]` and `R[2..L-1]`** by loading from
   `output/{l}loop/` if they already exist, or by invoking itself (this
   is the recursive part — it ensures the lower-loop prerequisites are
   present).
3. **Compute the boundary** (`compute_rhs.hpp::compute_boundary`):
   - `E1^n` via `shuffle_power` (sequential shuffle product).
   - `E1 · R_L` via `tensor_shuffle_product_parallel` (sequential).
   - Weighted sum via `tensor_add_weighted`.
4. **Expand the SEW collinear basis** (`A`) by contracting
   `SEW_<name>_basis.wxf` with `first_w{N}_basis.wxf` (highest weight
   first).
5. **Project `A` and the boundary to the divergent subspace** via
   `apply_colprojdiv_slots` (apply `colprojdiv_w1` to each 11-dim
   letter slot). This is required at `L ≥ 3`.
6. **Match positions** and **solve** `c · A_match = b_match` via
   `solve_linear_system` (modular RREF + reconstruct). See
   [04_collinear_solving.md](04_collinear_solving.md) for the matching
   logic.
7. **Compute `hepMHV_LL`** = `contract(solMHV_LL, SEW_basis, axis 1, 0)`.
   Per the design decision (Q3 in the original spec): use the SEW basis
   directly — do **not** apply `colprojdiv` to `hepMHV` or `E_L`.
8. **Expand `hepMHV_LL` to `E_L`** (rank `2L`, dims `11^2L`) via
   `expand_hepmhv`.
9. **Compute `R_L = E_L - boundary`** and save.
10. **Verify `R_L` is divergent-free** via the indicator-vector method:
    collect distinct letter indices, build an 11-dim indicator, project
    with `colprojdiv_w1`. If zero, `R_L = R*` is divergent-free.

## Inputs

- `data/E1.wxf` (the one-loop seed).
- All other `data/` seeds (transitively, via `--project`).
- Lower-loop outputs from `output/{l}loop/` (if they exist; otherwise
  computed recursively).

## Outputs (to `output/<L>loop/`)

Note: directory names use a digit prefix (`2loop`, `3loop`, `4loop`,
`5loop`), **not** `twoloop`/`threeloop`.

| File | Dims | Description |
|------|------|-------------|
| `solMHV_LL.wxf` | `1 × n_unknowns` | Solution coefficient vector |
| `hepMHV_LL.wxf` | `FEC_F_basis × 11` | Contracted solution (rank 2) |
| `E_LL.wxf` | `11^2L` | Expanded collinear projection (rank `2L`) |
| `R_LL.wxf` | `11^2L` | Remainder `E_L - boundary` (rank `2L`) |
| `boundary_LL.wxf` | `11^2L` | The boundary tensor (rank `2L`) |

Also writes `output/oneloop/E1.wxf` (a copy of `data/E1.wxf`) on the
first run, and triggers writes to `output/collinear/` via `--project`.

## Key files

- `compute_rhs.cpp` — CLI parsing, loop-order derivation, path resolution.
- `compute_rhs.hpp` — `compute_rhs_for_loop`, `compute_boundary`,
  `shuffle_power`, `expand_hepmhv`, `apply_colprojdiv_slots`,
  `ensure_fec_tensors`, `ensure_sew_basis`.
- `tensor_shuffle.h` — `tensor_shuffle_product_parallel` (sequential
  variant used for boundary computation; the parallel variant is
  incorrect).

## Conventions

- **`--data-dir` / `--output-dir` are exposed** (unlike `bootstrap`'s
  `--project` / `--solve-*` modes). Default to `<exec_dir>/data` and
  `<exec_dir>/output`.
- **Subprocess invocation**: `compute_rhs` shells out to `./bootstrap
  --extend`, `--sew`, `--project` to generate missing prerequisites.
  **Caveat**: it does NOT currently propagate `--data-dir` /
  `--output-dir` to these subprocesses — they will use their own defaults
  (`<exec_dir>/data`, `<exec_dir>/output`). This is a known limitation
  for the multi-project layout; see `README.md` §"Multi-Project Layout".
- **Sequential shuffle product**: always pass `pool = nullptr` to
  `tensor_shuffle_product_parallel` for boundary computation. The
  parallel variant produces incorrect results.
- **Boundary formulas are hardcoded** for `L = 2..5` in
  `compute_boundary`. Adding `L = 6` requires extending this function.

## Smoke test

```bash
# L=2 (SEW_3p1): computes E2, R2, boundary_2L
./compute_rhs --target SEW_3p1

# L=3 (SEW_5p1): computes E3, R3, boundary_3L (requires L=2 outputs)
./compute_rhs --target SEW_5p1
```

## Verified status

- **L=2**: boundary `E1²/2` (68 nnz before scaling), `E2` (rank 4),
  `R2` (rank 4, divergent-free).
- **L=3**: boundary `E1³/6 + E1·R2` (5894 nnz), `E3` (rank 6, 11606 nnz),
  `R3` (rank 6, 10461 nnz, divergent-free — contains only letters
  `{2,3,4,5,6,7,8,9,10}`).
- Solution at L=3: `c[0] = -24, c[1] = 2` (unique, all 32 constraints
  verified).

## Pitfalls

1. **`B.dims()` returns by value**: calling `B.dims().begin()` and
   `B.dims().end()` separately creates dangling iterators. Capture into
   a local first: `auto bdims = B.dims();`. (Fixed in `tensor_shuffle.h`.)
2. **`insert_add` corruption**: the COO `insert_add` (ordered insert)
   corrupts tensor dims for large results. Use `unordered_map` +
   `push_back` instead. (Fixed in `tensor_shuffle.h`.)
3. **Divergent-subspace projection is required at `L ≥ 3`**: because
   `E1` has divergent-letter entries. See
   [04_collinear_solving.md](04_collinear_solving.md).
4. **Matching must preserve the sew axis**: do not collapse `A` into
   `std::map<key, T>`. See [04_collinear_solving.md](04_collinear_solving.md).
