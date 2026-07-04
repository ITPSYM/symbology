# Changelog

All notable changes to this repository are documented here. Dates use
ISO 8601 (YYYY-MM-DD) and the local timezone is Asia/Shanghai.

## [2026-07-04] Multi-project path support (--data-dir / --output-dir)

### Summary
Threaded `--data-dir` / `--output-dir` through every pipeline mode so a
single checkout can serve multiple bootstrap projects (e.g. `E6`, `E7`,
`D5`) without editing source files. Each project lives in its own
`data_<PROJECT>/` + `output_<PROJECT>/` pair; the default (`data/` +
`output/`) remains backward compatible.

### Modified
- **`bootstrap.cpp`** — added `--data-dir` / `--output-dir` parsing and
  dispatch for `--project`, `--solve-symmetry`, `--solve-collinear`.
  Relative paths are resolved against the executable directory (same
  convention as `--condition` / `--first` / `--last` / `--output`).
  `--extend` / `--sew` ignore the flags (they use explicit `-c`/`-f`/`-l`/`-o`).
- **`projection.hpp`** — `run_projection_pipeline` signature changed from
  `base_path` to `data_dir` + `output_dir`; removed the internal
  `base / "data"` / `base / "output"` resolution.
- **`solve_symmetry.hpp`** — `run_symmetry_solver` signature changed from
  `base_path` to `data_dir` + `output_dir`; forwards both to
  `run_projection_pipeline`.
- **`compute_rhs.hpp`** — added `find_dlogmat(data_dir)` helper that scans
  for `dlogmat_*.wxf` instead of hardcoding `dlogmat_E6.wxf`; added a
  project-aware `run_bootstrap_cmd(cmd, data_dir, output_dir)` overload
  that appends `--data-dir <abs>` / `--output-dir <abs>` to every subprocess
  invocation so `--extend` / `--sew` / `--project` auto-invocations use the
  same project directories.
- **`compute_rhs.cpp`** — relative `--data-dir` / `--output-dir` now
  resolve against the executable directory (matching `bootstrap.cpp` and
  `inspect_tensors.cpp`); previously they were left relative to cwd.
- **`inspect_tensors.cpp`** — added `--output-dir` (and `--data-dir` for
  symmetry); `main` now takes `argc`/`argv` and resolves paths against the
  executable directory.
- **`run_workflow.sh`** — parameterized with `PROJECT` env var
  (`data_<PROJECT>/` + `output_<PROJECT>/`); auto-detects
  `dlogmat_*.wxf` in the data dir instead of hardcoding `dlogmat_E6.wxf`.
- **`run_projection.sh`**, **`run_solve.sh`** — parameterized with
  `PROJECT` env var; pass `--data-dir` / `--output-dir` to `bootstrap`.
- **`solve_collinear.hpp`** — fixed a misleading log message that printed
  `colprojdiv_SEW_SEW_3p1.wxf` (the actual file is `colprojdiv_SEW_3p1.wxf`
  because `sew_name` already includes the `SEW_` prefix).
- **`.gitignore`** — added `output_*/` to cover multi-project output dirs.
- **`README.md`** — removed the "Current limitation" caveat; documented
  the now-functional multi-project layout, the `PROJECT` env var, and the
  `--data-dir` / `--output-dir` flags for all three executables.

### Verified
- Default paths (`./compute_rhs --target SEW_3p1`): L=2 solve succeeds,
  `c[0] = 8`, all 8 constraints verified, `R2` divergent-free.
- Custom paths (`cp -r data data_test && ./compute_rhs --target SEW_3p1
  --data-dir data_test --output-dir output_test`): full pipeline runs
  end to end (subprocess `--extend` / `--sew` / `--project` all receive
  absolute `--data-dir` / `--output-dir`); produces byte-identical
  outputs (CRC32 match for `SEW_3p1_basis.wxf`, `E2.wxf`, `R2.wxf`,
  `boundary_2L.wxf`, `hepMHV_2L.wxf`).

## [2026-07-04] Standalone RHS computation module and divergent-subspace solve

### Added
- **`compute_rhs.cpp` / `compute_rhs.hpp`** — a standalone executable that
  recursively computes the collinear boundary (`E[L]`, `R[L]`, `boundary_LL`)
  from loop 2 up to a target SEW (e.g. `SEW_5p1` = 3-loop). Couples with
  the main `bootstrap` module to generate any missing SEW basis files.
  Exposes `--target`, `--data-dir`, `--output-dir` flags.
- **`tensor_shuffle.h`** — shuffle product implementation for the boundary
  formulas (`E1^n`, `E1·R2`, `E2^2`, etc.). Includes both sequential and
  parallel variants; the sequential variant is used for boundary
  computation because the parallel variant was observed to produce
  incorrect results for `E1^2/2`.
- **`inspect_tensors.cpp`** — a small diagnostic tool that prints the
  contents of `output/oneloop/E1.wxf` and `output/2loop/boundary_2L.wxf`
  for quick sanity checks.
- **`data/E1.wxf`** — the one-loop collinear seed tensor (renamed from the
  archive's `coloneloop.wxf`). Rank-2, `11 × 11`, 5 non-zero entries.
- **`data/DESCRIPTION.md`** — describes every seed file under `data/`.
- **`skills/`** — per-module skill documentation, structured for reuse by
  any AI agent (see `skills/README.md`).
- **`CHANGELOG.md`** (this file).

### Modified
- **`Makefile`** — added build rules for `compute_rhs` and `inspect_tensors`
  (each depends on its own `.cpp` plus the shared headers).
- **`bootstrap.cpp`** — added `--project`, `--solve-symmetry`, and
  `--solve-collinear` modes dispatched to the new modules.
- **`solve_collinear.hpp`** — added the divergent-subspace projection step
  (`apply_colprojdiv_slots`) before matching `A` against the boundary, and
  rewrote the matching logic to preserve the sew axis in `A_match` (the
  prior `std::map<key,T>` overwrite dropped duplicate sew entries and broke
  the 2-unknown system at `L = 3`).
- **`linear_solve.hpp`** — adapted the linear solver to accept the
  sew-axis-preserving `A_match` / `b_match` tensors.
- **`tensor_expand.hpp`** — fixed the `expansion_perm` slot ordering so
  the new letter is placed immediately after the FEC axis (preserving the
  natural ascending weight order `[w1, w2, w3, w4]`); the prior code
  reversed the slot order to `[w1, w4, w3, w2]`.
- **`.gitignore`** — added `compute_rhs` and `inspect_tensors` to the
  ignored executables list.
- **`README.md`** — added a "Multi-Project Layout" section, per-module
  documentation, and updated the CLI reference.

### Verified
- L=2 (`SEW_3p1`): unique solution `c[0] = 8`, all 44 constraints
  satisfied; `R2` is divergent-free.
- L=3 (`SEW_5p1`): unique solution `c[0] = -24, c[1] = 2`, all 32
  constraints satisfied; `R3` contains only letters `{2,3,4,5,6,7,8,9,10}`
  (no divergent letters `{0,1}`), so `R3 = R*` is divergent-free as
  required by the collinear constraint.

## [2026-07-04] Collinear solver with projection chain, expansion, and linear solving

### Added
- **`solve_collinear.hpp`** — collinear solver implementing:
  - `colprojfin` and `colprojdiv` projection chain (FEC and SEW levels),
  - Tensor expansion by contracting with bases,
  - Non-homogeneous linear equation solving for `c·A = boundary`,
  - Indicator-vector verification that `R* = c·A - boundary` is
    divergent-free.
- **`linear_solve.hpp`** — generic sparse linear solver over `rat_t` using
  modular reconstruction (RREF over `Z / 2^61`, then reconstruct).
- **`tensor_expand.hpp`** — universal expansion function that contracts a
  tensor with a chain of bases (highest weight first).

### Modified
- **`Makefile`** — added the new headers to the `bootstrap` build rule.
- **`bootstrap.cpp`** — added the `--solve-collinear` mode, requiring
  `--target` and `--rhs` (exits with code 1 if `--rhs` is missing; `--rhs 0`
  means an empty RHS tensor).

## [2026-07-04] Symmetry solving module

### Added
- **`solve_symmetry.hpp`** — symmetry solver that computes the invariant
  subspace of a target's projection (cyclic, flip, parity). Writes
  `<target>_invariant.wxf` under `output/<symmetry>/`.

### Modified
- **`bootstrap.cpp`** — added the `--solve-symmetry` mode, requiring
  `--symmetry` and `--target`.

## [2026-07-04] Universal projection matrix module

### Added
- **`projection.hpp`** — universal projection pipeline supporting both
  collinear (dimension-shrinking) and symmetry (dimension-preserving)
  projections. Handles starting seeds (`colmat42`, `cycrepmat`, `flipmat`,
  `paritymat`) with `FEC_1`/`LEC_1` truncations. Stores the reduced basis
  for collinear projections (`first_w{N}_basis.wxf`,
  `last_w{N}_basis.wxf`, `SEW_{name}_basis.wxf`).
- **`run_projection.sh`** — driver script for the projection pipeline.
- **`run_solve.sh`** — driver script for the symmetry solver.

### Modified
- **`bootstrap.cpp`** — added the `--project` mode, requiring `--symmetry`
  and `--target`.
- **`run_workflow.sh`** — updated to use the new naming convention.

## [2026-07-04] CSR/COO consistency and projection logging

### Modified
- **`projection.hpp`** — all tensor functions now return CSR (matching the
  `bootstrap.cpp` and `bootstrap.hpp` convention); COO overloads removed.
- **`run_projection.sh`** — now uses `tee` to log to `logs/` and writes a
  `summary.txt` (with tensor dimensions, `nnz`, and CRC32) under each
  `output/<symmetry>/` directory.
- **`bootstrap.cpp`** — added CRC32 printing on every file read/write.

## [2026-07-04] Initial commit and macOS build setup

### Added
- Initial `bootstrap` executable with `--extend`, `--sew`, `--induce` modes.
- `bootstrap.hpp` — tensor layouts, forward/backward extension, sewing,
  timing, helper routines.
- `data/dlogmat_E6.wxf`, `data/FEC_1.wxf`, `data/LEC_1.wxf` — the `E6`
  seed tensors.
- `SparseRREF/` — included as a non-vendored dependency (clone separately).
- `README.md` — build and usage documentation.
- `.gitignore` — added `output/`, `logs/`, `temp/`, `tmp/`, `.trae/`,
  `.local/`, and the `bootstrap` executable.
- Updated `README.md` with macOS Homebrew GCC build instructions (the
  libc++ limitation on `std::execution::par` and `std::chrono::zoned_time`
  requires `g++-14`).
