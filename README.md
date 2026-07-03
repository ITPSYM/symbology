# Symbology

AI assisted symbology study

This repository contains the working version of the symbol bootstrap workflow.
The first stage focuses on the `E6` heptagon bootstrap data that was previously
kept in `bootstrap_E6_archive/`.

## Current Workflow

The main executable is `bootstrap`, built from:

- `bootstrap.cpp`: command-line parsing, WXF input/output, CRC32 printing, and
  mode dispatch.
- `bootstrap.hpp`: tensor layouts, forward/backward extension, sewing, timing,
  and small helper routines.

Build with:

```bash
make
```

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

- `--extend`: grow either forward (`-f/--first`) or backward (`-l/--last`) data
  by one weight.
- `--sew`: combine a forward tensor and a backward tensor into a sewing matrix.
- `--induce`: reserved for future induced-transformation workflows.
- `-c/--condition`: condition tensor, currently `data/dlogmat_E6.wxf`.
- `-o/--output`: output WXF file.

Thread count is chosen automatically by `SparseRREF`.

## Tensor Files

Tracked seed data:

- `data/dlogmat_E6.wxf`: RREF-reduced `E6` adjacency/integrability condition
  tensor. It corresponds to `bootstrap_E6_archive/dlogmatE6RREF.wxf` and has
  dimensions `{42, 42, 1191}`.
- `data/FEC_1.wxf`: forward expansion coefficient seed, copied from
  `bootstrap_E6_archive/FCC_1.wxf`; layout `{basis_w, basis_{w-1}, letter}`.
- `data/LEC_1.wxf`: backward expansion coefficient seed, obtained from
  `bootstrap_E6_archive/LCC_1.wxf` by transposing to the new layout
  `{basis_w, letter, basis_{w-1}}`.

Generated files:

- `output/FEC_w.wxf`: forward expansion coefficients.
- `output/LEC_w.wxf`: backward expansion coefficients.
- `output/SEW_fp l.wxf`: sewing matrices with layout
  `{sew_basis, FEC_f_basis, LEC_l_basis}`.
- `logs/*.log`: stdout/stderr logs for each workflow step.

`output/`, `logs/`, the compiled `bootstrap` executable, `temp/`, and `tmp/`
are ignored by git.

## Workflow Script

`run_workflow.sh` is the main driver.

Smoke workflow:

```bash
./run_workflow.sh smoke
```

This runs locally tractable checks:

- `FEC_1 -> FEC_6`
- `LEC_1 -> LEC_4`
- `SEW_2p2`, `SEW_3p1`, `SEW_4p2`, `SEW_5p1`

Full workflow:

```bash
./run_workflow.sh full
```

This lists the intended first-stage path through `FEC_8`, `LEC_5`, and
`SEW_8p2`. It is not expected to be practical on a normal local machine.

## Format Notes

`bootstrap` writes SparseRREF-native WXF. Some archive files were exported
through Mathematica, so byte-level CRC32 values can differ even when tensor
contents are identical.

For forward files, the current workflow has been checked as follows:

1. Generate `output/FEC_2.wxf` through `output/FEC_6.wxf`.
2. Roundtrip each file with the top-level `wxf_roundtrip.wls`.
3. Compare CRC32 with `bootstrap_E6_archive/FCC_2_rref.wxf` through
   `FCC_6_rref.wxf`.

After Mathematica roundtrip, all checked `FEC_2..FEC_6` CRC32 values match the
archive. `LEC` files intentionally use a different axis order from the old
`LCC` files.
