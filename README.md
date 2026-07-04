
# Symbology

Working prototype for symbol-space bootstrap experiments on the `E6` heptagon. The codebase covers the full bootstrap pipeline: recursive first-entry/last-entry growth and sewing (`bootstrap --extend` / `--sew`), universal projection matrices (`--project`), symmetry invariant-subspace solving (`--solve-symmetry`), collinear constraint solving (`--solve-collinear`), and recursive RHS (collinear boundary) computation (`compute_rhs`). All sparse rational linear algebra goes through [`SparseRREF`](https://github.com/munuxi/SparseRREF).

## What You Need

The executable is a C++20 program. A local build needs:

- a C++20 compiler with `<format>` and chrono time-zone support;
- `make`;
- FLINT and GMP;
- TBB;
- mimalloc;
- Git, to fetch this repository and `SparseRREF`;
- optional: Wolfram/Mathematica, useful for inspecting or round-tripping WXF `SparseArray` data.

`SparseRREF` is not vendored as a submodule. Put a checkout of its current default branch at repository root, so the headers live under `SparseRREF/`.

```bash
git clone https://github.com/ITPSYM/symbology.git
cd symbology
git clone https://github.com/munuxi/SparseRREF.git
```

## Dependency Setup

### macOS

#### Why Homebrew GCC is required

`bootstrap` uses `std::chrono::zoned_time` (for timestamps) and `SparseRREF` uses `std::execution::par` (parallel algorithms) unconditionally. On macOS **every clang toolchain ships libc++**, and libc++ gates both of these behind build-time feature flags that no available distribution enables — Apple clang, Homebrew clang, and conda-forge clang all fail with `no member named 'par' in namespace 'std::execution'` and `no member named 'zoned_time' in namespace 'std::chrono'`. The `-D_LIBCPP_HAS_PARALLEL_ALGORITHMS` macro does not help (it is decided when libc++ itself is compiled, not at use-site). Conda-forge `gxx` on macOS is also a clang wrapper driving libc++, so it has the same problem.

The reliable route is **Homebrew GCC** (`g++-14`), which brings its own libstdc++ where `std::execution::par` and `zoned_time` are available unconditionally.

#### Setup

Install the libraries and Homebrew GCC:

```bash
xcode-select --install
brew install flint gmp tbb mimalloc gcc
```

Then build with `g++-14` (note: this bypasses the default `make` rule, which uses the system `g++`/`clang++`):

```bash
make CXX=g++-14
```

If `g++-14` is not on `PATH`, point at it explicitly, for example:

```bash
make CXX=/opt/homebrew/bin/g++-14
```

#### If `g++-14` cannot find the libraries

Pass the Homebrew include and library paths explicitly:

```bash
make CXX=g++-14 \
  CXXFLAGS="-O3 -std=c++20 -I. -I/opt/homebrew/include" \
  LDLIBS="-L/opt/homebrew/lib -Wl,-rpath,/opt/homebrew/lib -lflint -lgmp -lmimalloc -ltbb"
```

#### Linux note

On Linux (the archive's original target) libstdc++ is the default and the plain `make` rule works once FLINT, GMP, TBB, and mimalloc are installed — none of the macOS libc++ issues arise.

### Ubuntu/Debian

```bash
sudo apt-get update
sudo apt-get install build-essential make git libflint-dev libgmp-dev libtbb-dev libmimalloc-dev
make
```

Package names can vary slightly across distributions. The important libraries are FLINT, GMP, TBB, and mimalloc.

## Build

There are three executables, all built with `make`:

- `bootstrap` — the main dispatcher (extension, sewing, projection, symmetry solving, collinear solving). Built from `bootstrap.cpp` plus the shared headers.
- `compute_rhs` — standalone recursive RHS (collinear boundary) computation. Built from `compute_rhs.cpp`.
- `inspect_tensors` — small diagnostic tool that prints tensor contents. Built from `inspect_tensors.cpp`.

Each rule in the `Makefile` declares its own header dependencies, so `make` will only rebuild what has changed. The shared headers are: `bootstrap.hpp`, `projection.hpp`, `solve_symmetry.hpp`, `solve_collinear.hpp`, `linear_solve.hpp`, `tensor_expand.hpp`, `tensor_shuffle.h`.

Build from the repository root:

```bash
make               # builds bootstrap only (default target)
make compute_rhs   # builds the RHS computation module
make inspect_tensors  # builds the diagnostic tool
make all           # same as `make bootstrap`
```

Clean the binaries with:

```bash
make clean
```

## Multi-Project Layout

By default, every executable reads seed tensors from `data/` and writes outputs to `output/`. This is the `E6` problem that ships with the repository. For a different symmetry group or a different bootstrap project, use a per-project directory pair:

```
data_<PROJECT>/      # e.g. data_E7/, data_D5/
output_<PROJECT>/    # e.g. output_E7/, output_D5/
```

The default (no tag) is `data/` + `output/`, which is backward compatible.

### Conventions

- Inside `data_<PROJECT>/`, the seed files keep the same **roles** as in `data/` (see `data/DESCRIPTION.md`), but encode the symmetry group in the filename where it matters: `dlogmat_<group>.wxf` (e.g. `dlogmat_E7.wxf`), `FEC_1.wxf`, `LEC_1.wxf`, `colmat<N>.wxf` (where `<N>` is the FEC weight-1 dimension — `42` for `E6`), `colprojdiv.wxf`, `colprojfin.wxf`, `<group>repmat_*.wxf`, `E1.wxf`.
- All executables accept `--data-dir <dir>` and `--output-dir <dir>`. Relative paths are resolved against the executable directory. The flags are threaded through every pipeline mode: `bootstrap --project` / `--solve-symmetry` / `--solve-collinear`, `compute_rhs`, and `inspect_tensors`. When `compute_rhs` shells out to `./bootstrap --extend` / `--sew` / `--project` to generate missing prerequisites, it passes the absolute `--data-dir` / `--output-dir` to the subprocess so the same project directories are used end to end.
- The driver scripts (`run_workflow.sh`, `run_projection.sh`, `run_solve.sh`) honor a `PROJECT=<name>` environment variable: setting `PROJECT=E7` makes them use `data_E7/` + `output_E7/`. With `PROJECT` unset they default to `data/` + `output/` (backward compatible).
- `run_workflow.sh` auto-detects the condition tensor as `dlogmat_*.wxf` inside the data directory, so it generalizes to other symmetry groups without editing the script.

## Minimal Smoke Test

After building, run one forward step, one backward step, and one sewing step:

```bash
./bootstrap --extend -c data/dlogmat_E6.wxf -f data/FEC_1.wxf -o output/FEC_2.wxf
./bootstrap --extend -c data/dlogmat_E6.wxf -l data/LEC_1.wxf -o output/LEC_2.wxf
./bootstrap --sew -c data/dlogmat_E6.wxf -f output/FEC_2.wxf -l output/LEC_2.wxf -o output/SEW_2p2.wxf
```

A successful run prints tensor ranks, dimensions, nonzero counts, RREF timing information, and CRC32 values for the files it reads and writes. `output/` is generated locally and is ignored by git.

For a longer but still local workflow, use:

```bash
./run_workflow.sh smoke
```

This runs:

- `FEC_1 -> FEC_6`;
- `LEC_1 -> LEC_4`;
- `SEW_2p2`, `SEW_3p1`, `SEW_4p2`, `SEW_5p1`.

The full workflow is available as:

```bash
./run_workflow.sh full
```

It lists the intended first-stage path through `FEC_8`, `LEC_5`, and `SEW_8p2`, but it is not expected to be practical on a normal laptop.

### Projection and solving smoke test

After running the bootstrap smoke test above (or `run_workflow.sh smoke`), the SEW tensors exist. Run the projection, symmetry solving, and RHS computation modules:

```bash
./bootstrap --project --symmetry collinear --target SEW_5p1
./bootstrap --solve-symmetry --symmetry cyclic --target SEW_5p1
./compute_rhs --target SEW_3p1    # 2-loop (requires only E1)
./compute_rhs --target SEW_5p1    # 3-loop (requires L=2 outputs)
```

Or via the driver scripts:

```bash
./run_projection.sh SEW_5p1   # all four symmetries
./run_solve.sh SEW_5p1        # cyclic, flip, parity
```

A successful `compute_rhs --target SEW_5p1` run prints the unique solution (`c[0] = -24, c[1] = 2`), verifies all 32 constraints, and confirms `R3` is divergent-free (no entries at letters `{0, 1}`).

## CLI

### Core bootstrap (`./bootstrap`)

Forward extension:

```bash
./bootstrap --extend -c data/dlogmat_E6.wxf -f data/FEC_1.wxf -o output/FEC_2.wxf
```

Backward extension:

```bash
./bootstrap --extend -c data/dlogmat_E6.wxf -l data/LEC_1.wxf -o output/LEC_2.wxf
```

Sewing:

```bash
./bootstrap --sew -c data/dlogmat_E6.wxf -f output/FEC_4.wxf -l output/LEC_2.wxf -o output/SEW_4p2.wxf
```

Projection matrix (collinear or symmetry):

```bash
./bootstrap --project --symmetry collinear --target SEW_5p1
./bootstrap --project --symmetry cyclic    --target SEW_5p1
```

Symmetry invariant subspace solver:

```bash
./bootstrap --solve-symmetry --symmetry cyclic --target SEW_5p1
```

Collinear constraint solver:

```bash
./bootstrap --solve-collinear --target SEW_5p1 --rhs output/3loop/boundary_3L.wxf
./bootstrap --solve-collinear --target SEW_3p1 --rhs 0   # empty RHS
```

Options:

- `--extend`: grow either forward (`-f/--first`) or backward (`-l/--last`) data by one weight;
- `--sew`: combine a forward tensor and a backward tensor into a sewing matrix;
- `--induce`: reserved for future induced-transformation workflows;
- `--project`: run the universal projection pipeline (requires `--symmetry`, `--target`);
- `--solve-symmetry`: compute the invariant subspace of a target's projection (requires `--symmetry`, `--target`);
- `--solve-collinear`: finite/divergent split + expansion + linear solve (requires `--target`, `--rhs`, `--projection`; `--basis` optional);
- `--symmetry <collinear|cyclic|flip|parity>`: symmetry name for `--project` / `--solve-symmetry`;
- `--target <SEW_FpL|FEC_W|LEC_W>`: target name (e.g. `SEW_5p1`, `FEC_3`, `LEC_2`);
- `--rhs <rhs.wxf>` or `--rhs 0`: RHS path for `--solve-collinear`; `"0"` means an all-zero RHS constructed in-memory. Missing → exit code 1;
- `--projection <finite|divergent>`: which projection to apply in `--solve-collinear` (required — no default);
- `--basis <basis.wxf>`: expansion basis file (repeatable; highest weight first). Auto-detected as `first_w{N}_basis.wxf` if omitted;
- `--data-dir <dir>`: data directory with seed files (default: `<exec_dir>/data`). Used by `--project`, `--solve-symmetry`, `--solve-collinear`; ignored by `--extend` / `--sew` (which use explicit `-c`/`-f`/`-l`/`-o` paths);
- `--output-dir <dir>`: output directory (default: `<exec_dir>/output`). Same scope as `--data-dir`;
- `-c/--condition`: condition tensor, currently `data/dlogmat_E6.wxf`;
- `-f/--first`, `-l/--last`, `-o/--output`: input/output file paths;
- `-h/--help`: print usage.

Thread count is chosen automatically by `SparseRREF`.

### RHS computation (`./compute_rhs`)

```bash
./compute_rhs --target SEW_3p1    # 2-loop: computes E2, R2, boundary_2L
./compute_rhs --target SEW_5p1    # 3-loop: computes E3, R3, boundary_3L (requires L=2 outputs)
./compute_rhs --target SEW_5p1 --data-dir data_E7 --output-dir output_E7   # multi-project
```

Options:

- `--target <SEW_FpL>`: target SEW name (required). Loop order `L = (F+L)/2`; supported `L = 2..5`;
- `--data-dir <dir>`: data directory with seed files (default: `<exec_dir>/data`);
- `--output-dir <dir>`: output directory (default: `<exec_dir>/output`);
- `-h/--help`: print usage.

### Inspection (`./inspect_tensors`)

```bash
./inspect_tensors                                  # default: reads from ./output/
./inspect_tensors --output-dir output_E7           # multi-project
```

Options:

- `--output-dir <dir>`: output directory (default: `<exec_dir>/output`);
- `--data-dir <dir>`: accepted for symmetry with the other tools (not used by `inspect_tensors`);
- `-h/--help`: print usage.

Reads `<output-dir>/oneloop/E1.wxf` and `<output-dir>/2loop/boundary_2L.wxf`.

## Tensor Files

Tracked seed data (see `data/DESCRIPTION.md` for a complete listing with dimensions):

- `data/dlogmat_E6.wxf`: RREF-reduced `E6` adjacency/integrability condition tensor. It corresponds to `bootstrap_E6_archive/dlogmatE6RREF.wxf` and has dimensions `{42, 42, 1191}`.
- `data/FEC_1.wxf`: forward expansion coefficient seed, copied from `bootstrap_E6_archive/FCC_1.wxf`; layout `{basis_w, basis_{w-1}, letter}`.
- `data/LEC_1.wxf`: backward expansion coefficient seed, obtained from `bootstrap_E6_archive/LCC_1.wxf` by transposing to the new layout `{basis_w, letter, basis_{w-1}}`.
- `data/colmat42.wxf`: collinear seed on the FEC weight-1 space (`42 × 2`).
- `data/cycrepmat.wxf`, `data/fliprepmat.wxf`, `data/parityrepmat.wxf`: cyclic / flip / parity symmetry representation matrices on the FEC weight-1 space (`42 × 42`).
- `data/colprojdiv.wxf`: weight-1 colprojdiv seed (`11 × 2`); projects each letter slot to its 2-dim divergent subspace.
- `data/colprojfin.wxf`: weight-1 colprojfin seed (`11 × 9`); projects each letter slot to its 9-dim finite subspace.
- `data/E1.wxf`: one-loop collinear seed tensor (`11 × 11`, 5 nnz); renamed from the archive's `coloneloop.wxf`. Used by `compute_rhs`.

Generated files (under `output/`):

- `output/FEC_w.wxf`: forward expansion coefficients;
- `output/LEC_w.wxf`: backward expansion coefficients;
- `output/SEW_fpl.wxf`: sewing matrices with layout `{sew_basis, FEC_f_basis, LEC_l_basis}`;
- `output/collinear/`: collinear projection chain — `first_w{N}.wxf`, `last_w{N}.wxf`, `first_w{N}_basis.wxf`, `last_w{N}_basis.wxf`, `SEW_<name>_basis.wxf`, `colprojfin_w{N}.wxf`, `colprojdiv_w{N}.wxf`, `colprojfin_<sew_name>.wxf`, `colprojdiv_<sew_name>.wxf`, plus a `summary.txt`;
- `output/cyclic/`, `output/flip/`, `output/parity/`: symmetry projections — `first_w{N}.wxf`, `last_w{N}.wxf`, `SEW_<name>.wxf`, `<target>_invariant.wxf`, plus a `summary.txt`;
- `output/oneloop/E1.wxf`: copy of `data/E1.wxf` (written by `compute_rhs`);
- `output/{L}loop/` (digit prefix — `2loop`, `3loop`, `4loop`, `5loop`): per-loop results from `compute_rhs` — `solMHV_LL.wxf`, `hepMHV_LL.wxf`, `E_LL.wxf`, `R_LL.wxf`, `boundary_LL.wxf`;
- `logs/*.log`: stdout/stderr logs for each workflow step.

`output/`, `output_*/`, `logs/`, the compiled `bootstrap` / `compute_rhs` / `inspect_tensors` executables, `temp/`, and `tmp/` are ignored by git.

## Skills and Changelog

- `skills/` holds per-module reference documents (concise, model-agnostic) for AI agents and new contributors. Start at `skills/README.md`.
- `CHANGELOG.md` records all notable changes (new files, modified files, new functionality) grouped by date.
- `data/DESCRIPTION.md` describes every seed file under `data/`.

## Format Notes

`bootstrap` writes SparseRREF-native WXF. Some archive files were exported through Mathematica, so byte-level CRC32 values can differ even when tensor contents are identical.

For forward files, the current workflow has been checked as follows:

1. Generate `output/FEC_2.wxf` through `output/FEC_6.wxf`.
2. Roundtrip each file with the top-level `wxf_roundtrip.wls`.
3. Compare CRC32 with `bootstrap_E6_archive/FCC_2_rref.wxf` through `FCC_6_rref.wxf`.

After Mathematica roundtrip, all checked `FEC_2..FEC_6` CRC32 values match the archive. `LEC` files intentionally use a different axis order from the old `LCC` files.

## Troubleshooting

- If `make` cannot find `SparseRREF/sparse_mat.h`, clone `SparseRREF` into the repository root.
- If `git` fails on macOS with an `xcode-select` error, install the Apple command-line tools or use the conda-forge setup above.
- If compilation fails on macOS with `no member named 'par' in namespace 'std::execution'` or `no member named 'zoned_time' in namespace 'std::chrono'`, you are hitting the libc++ limitation described in the macOS section above. A newer clang will **not** fix it — build with Homebrew GCC instead: `make CXX=g++-14`.
- If the linker cannot find FLINT, GMP, TBB, or mimalloc, check that the matching include and library paths are visible to `make`.
- If byte-level WXF CRC32 values differ from archived Mathematica exports, compare after a Mathematica roundtrip rather than comparing raw SparseRREF-native WXF bytes.
