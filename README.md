
# Symbology

Working prototype for symbol-space bootstrap experiments. The current code focuses on recursive first-entry/last-entry growth and sewing for the `E6` heptagon data, using sparse rational linear algebra through [`SparseRREF`](https://github.com/munuxi/SparseRREF).

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

If you use Homebrew and the Apple command-line tools:

```bash
xcode-select --install
brew install flint gmp tbb mimalloc
```

Then build with:

```bash
make
```

On some macOS systems the system compiler or libc++ is too old for the C++20 features used here or for the latest `SparseRREF`. In that case, the most reproducible route is a conda-forge compiler environment:

```bash
micromamba create -n symbology -c conda-forge cxx-compiler make libflint gmp tbb mimalloc git
micromamba activate symbology
make CXX="$CXX"
```

If your compiler finds the wrong headers or libraries, pass the include and library paths explicitly, for example:

```bash
make CXX="$CXX" \
  CXXFLAGS="-O3 -std=c++20 -I. -I$CONDA_PREFIX/include" \
  LDLIBS="-L$CONDA_PREFIX/lib -Wl,-rpath,$CONDA_PREFIX/lib -lflint -lgmp -lmimalloc -ltbb"
```

### Ubuntu/Debian

```bash
sudo apt-get update
sudo apt-get install build-essential make git libflint-dev libgmp-dev libtbb-dev libmimalloc-dev
make
```

Package names can vary slightly across distributions. The important libraries are FLINT, GMP, TBB, and mimalloc.

## Build

The main executable is `bootstrap`, built from:

- `bootstrap.cpp`: command-line parsing, WXF input/output, CRC32 printing, and mode dispatch;
- `bootstrap.hpp`: tensor layouts, forward/backward extension, sewing, timing, and helper routines;
- `SparseRREF/`: sparse tensors, sparse RREF, rational arithmetic support, threading, and WXF support.

Build from the repository root:

```bash
make
```

Clean the binary with:

```bash
make clean
```

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

## CLI

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

Options:

- `--extend`: grow either forward (`-f/--first`) or backward (`-l/--last`) data by one weight;
- `--sew`: combine a forward tensor and a backward tensor into a sewing matrix;
- `--induce`: reserved for future induced-transformation workflows;
- `-c/--condition`: condition tensor, currently `data/dlogmat_E6.wxf`;
- `-o/--output`: output WXF file.

Thread count is chosen automatically by `SparseRREF`.

## Tensor Files

Tracked seed data:

- `data/dlogmat_E6.wxf`: RREF-reduced `E6` adjacency/integrability condition tensor. It corresponds to `bootstrap_E6_archive/dlogmatE6RREF.wxf` and has dimensions `{42, 42, 1191}`.
- `data/FEC_1.wxf`: forward expansion coefficient seed, copied from `bootstrap_E6_archive/FCC_1.wxf`; layout `{basis_w, basis_{w-1}, letter}`.
- `data/LEC_1.wxf`: backward expansion coefficient seed, obtained from `bootstrap_E6_archive/LCC_1.wxf` by transposing to the new layout `{basis_w, letter, basis_{w-1}}`.

Generated files:

- `output/FEC_w.wxf`: forward expansion coefficients;
- `output/LEC_w.wxf`: backward expansion coefficients;
- `output/SEW_fpl.wxf`: sewing matrices with layout `{sew_basis, FEC_f_basis, LEC_l_basis}`;
- `logs/*.log`: stdout/stderr logs for each workflow step.

`output/`, `logs/`, the compiled `bootstrap` executable, `temp/`, and `tmp/` are ignored by git.

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
- If compilation fails around `<format>`, `std::chrono::zoned_time`, or parallel algorithms in `SparseRREF`, use a newer C++20 compiler from conda-forge or another current toolchain.
- If the linker cannot find FLINT, GMP, TBB, or mimalloc, check that the matching include and library paths are visible to `make`.
- If byte-level WXF CRC32 values differ from archived Mathematica exports, compare after a Mathematica roundtrip rather than comparing raw SparseRREF-native WXF bytes.
