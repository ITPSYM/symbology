# Skills

Concise, model-agnostic reference documents for the `symbology` symbol-space
bootstrap codebase. Each file describes one module: its purpose, the CLI
entry point, the inputs/outputs, the key functions, and the conventions an
AI agent (or new contributor) must respect when modifying or extending it.

These documents describe *how to use and modify* the code; they are not
tutorials on the underlying physics. Read the relevant skill file before
working on a module.

## Index

| File | Module | Entry point |
|------|--------|-------------|
| [01_bootstrap.md](01_bootstrap.md) | Core bootstrap (extension + sewing) | `./bootstrap --extend` / `--sew` |
| [02_projection.md](02_projection.md) | Universal projection matrix module | `./bootstrap --project` |
| [03_symmetry_solving.md](03_symmetry_solving.md) | Symmetry invariant subspace solver | `./bootstrap --solve-symmetry` |
| [04_collinear_solving.md](04_collinear_solving.md) | Collinear constraint solver | `./bootstrap --solve-collinear` |
| [05_compute_rhs.md](05_compute_rhs.md) | Recursive RHS (boundary) computation | `./compute_rhs --target` |
| [06_sparserref_tensors.md](06_sparserref_tensors.md) | SparseRREF tensor library essentials | (library reference) |

## How to read these

- **Build first**: `make CXX=g++-14` on macOS (Homebrew GCC); plain `make`
  on Linux. See `README.md` for dependency setup.
- **Conventions**: every module returns CSR tensors; COO is used only
  internally for reshaping or for the shuffle product. WXF files are
  CSR-native. Logs go to `logs/`; summaries go to
  `output/<symmetry>/summary.txt`.
- **Multi-project**: see `README.md` §"Multi-Project Layout" for the
  `data_<PROJECT>/` + `output_<PROJECT>/` convention.

## Verified status (as of 2026-07-04)

| Module | Target | Status |
|--------|--------|--------|
| Bootstrap | `FEC_2..FEC_6`, `LEC_2..LEC_4`, `SEW_2p2..SEW_5p1` | Verified (CRC32 matches archive after Mathematica roundtrip) |
| Projection (collinear) | `SEW_5p1` | Verified (basis dims `2 × 546 × 11`, rank 7 expansion) |
| Projection (symmetries) | `SEW_5p1` cyclic/flip/parity | Verified (all three are identity `diag(1,1)` for `E6`) |
| Symmetry solving | `SEW_5p1` cyclic/flip/parity | Verified (invariant subspace == full space) |
| Collinear solving | `SEW_3p1` (L=2) | Verified (unique solution `c[0] = 8`, 8 constraints, `R2` divergent-free) |
| Collinear solving | `SEW_5p1` (L=3) | Verified (unique solution `c[0] = -24, c[1] = 2`, 32 constraints, `R3` divergent-free) |
| Compute RHS | L=2, L=3 | Verified end-to-end |

## Pitfalls (cross-cutting)

These are documented in detail in the per-module files, but listed here
because they affect every module:

1. **`B.dims()` returns by value.** Calling `B.dims().begin()` and
   `B.dims().end()` separately creates iterators from two different
   temporaries — a dangling-iterator bug. Capture into a local first:
   `auto bdims = B.dims();`.
2. **CSR has no `reshape` method.** Only COO has `reshape`. Convert
   CSR→COO→reshape→CSR when reshaping is needed.
3. **`tensor_shuffle_product_parallel` parallel variant is incorrect for
   boundary computation.** Always pass `pool = nullptr` for the
   sequential variant when computing `E1^n` or `E1·R_L` terms.
4. **`std::map<key, T>` overwrites duplicate keys.** When iterating a
   tensor that has multiple entries per key (e.g. the sew axis in `A`),
   do not collapse into a `std::map<key, T>` — iterate directly and
   preserve the full index.
5. **`expansion_perm` slot ordering.** After `tensor_contract`, the new
   letter must be placed immediately after the FEC axis to maintain the
   natural ascending weight order `[w1, w2, w3, w4]`. Placing it at the
   end reverses the slot order to `[w1, w4, w3, w2]` and breaks
   verification.
6. **`colprojdiv_w1` swaps the two divergent letters.** For `E6`:
   `colprojdiv[0,1] = colprojdiv[1,0] = 1`. Do not assume it is a
   projection onto the first coordinate.
7. **E1 has divergent-letter entries.** `E1[0,0] = -2` and
   `E1[1,1] = -2`. The boundary `E1^L / L!` therefore has divergent
   components, and solving `c·A = boundary` in the full space fails at
   `L ≥ 3`. Always project both `A` and the boundary to the divergent
   subspace (`apply_colprojdiv_slots`) before matching.
