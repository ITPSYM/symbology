# Skill 01 — Core bootstrap (extension + sewing)

## Purpose

Grow the forward (`FEC`) and backward (`LEC`) expansion coefficients one
weight at a time, and sew a forward/backward pair into a sewing matrix
(`SEW`). This is the foundational module — every other module depends on
the `FEC_*` / `LEC_*` / `SEW_*` tensors it produces.

## CLI entry points

```bash
# Forward extension: FEC_{w-1} -> FEC_w
./bootstrap --extend -c data/dlogmat_E6.wxf -f data/FEC_1.wxf -o output/FEC_2.wxf

# Backward extension: LEC_{w-1} -> LEC_w
./bootstrap --extend -c data/dlogmat_E6.wxf -l data/LEC_1.wxf -o output/LEC_2.wxf

# Sewing: FEC_F + LEC_L -> SEW_FpL
./bootstrap --sew -c data/dlogmat_E6.wxf -f output/FEC_4.wxf -l output/LEC_2.wxf -o output/SEW_4p2.wxf
```

## Flags

| Flag | Description |
|------|-------------|
| `--extend` | Extend a FEC (with `-f`) or LEC (with `-l`) one weight forward. |
| `--sew` | Sew a FEC and LEC into a SEW tensor (requires `-f` and `-l`). |
| `--induce` | Reserved for a future workflow stage; currently rejected. |
| `-c` / `--condition <file>` | Condition (dlog) matrix path. |
| `-f` / `--first <file>` | FEC input. |
| `-l` / `--last <file>` | LEC input. |
| `-o` / `--output <file>` | Output tensor path. |
| `-h` / `--help` | Print usage. |

## Inputs (from `data/`)

- `dlogmat_E6.wxf` — RREF-reduced `E6` condition tensor, dims
  `{42, 42, 1191}`. Read via `-c/--condition`.
- `FEC_1.wxf` — forward seed, dims `{7, 1, 42}`. Read via `-f` for the
  first `--extend` step.
- `LEC_1.wxf` — backward seed, dims `{14, 42, 1}`. Read via `-l` for the
  first `--extend` step.

## Outputs (to `output/`)

- `FEC_w.wxf` — forward expansion coefficients at weight `w`
  (`w = 2, 3, ...`).
- `LEC_w.wxf` — backward expansion coefficients at weight `w`.
- `SEW_FpL.wxf` — sewing matrix with layout `{sew_basis, FEC_F_basis,
  LEC_L_basis}`. Naming: `SEW_<F>p<L>.wxf` (e.g. `SEW_5p1.wxf` for
  FEC_5 + LEC_1).

## Key files

- `bootstrap.cpp` — CLI parsing, mode dispatch, WXF I/O, CRC32 printing.
- `bootstrap.hpp` — tensor layouts, forward/backward extension, sewing,
  timing helpers.

## Conventions

- **Path resolution**: relative paths are resolved against the executable's
  directory (`argv[0]` parent), not the current working directory. Pass
  absolute paths if in doubt.
- **Tensor layout**: every output is CSR-format WXF. CRC32 is printed on
  every read/write for traceability.
- **Naming**: `SEW_FpL` means `F` (forward weight) + `p` + `L` (last
  weight), e.g. `SEW_5p1` is the 3-loop sewing matrix (since
  `L = (F+L)/2 = (5+1)/2 = 3`).
- **Logging**: every workflow step should be tee'd to `logs/<name>.log`.

## Smoke test

```bash
./bootstrap --extend -c data/dlogmat_E6.wxf -f data/FEC_1.wxf -o output/FEC_2.wxf
./bootstrap --extend -c data/dlogmat_E6.wxf -l data/LEC_1.wxf -o output/LEC_2.wxf
./bootstrap --sew -c data/dlogmat_E6.wxf -f output/FEC_2.wxf -l output/LEC_2.wxf -o output/SEW_2p2.wxf
```

A successful run prints tensor ranks, dimensions, nonzero counts, RREF
timing, and CRC32 values for the files it reads and writes.

## Longer workflow

```bash
./run_workflow.sh smoke   # FEC_1 -> FEC_6, LEC_1 -> LEC_4, SEW_2p2..SEW_5p1
./run_workflow.sh full    # intended path through FEC_8, LEC_5, SEW_8p2 (not practical on a laptop)
```

## Verification

Forward files have been checked as follows:

1. Generate `output/FEC_2.wxf` through `output/FEC_6.wxf`.
2. Roundtrip each file with `wxf_roundtrip.wls` (Mathematica).
3. Compare CRC32 with `bootstrap_E6_archive/FCC_2_rref.wxf` through
   `FCC_6_rref.wxf`.

After Mathematica roundtrip, all checked `FEC_2..FEC_6` CRC32 values match
the archive. `LEC` files intentionally use a different axis order from the
old `LCC` files.
