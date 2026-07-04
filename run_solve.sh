#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

# Multi-project support: set PROJECT=<name> (e.g. PROJECT=E7) to use
# data_<name>/ and output_<name>/ instead of the default data/ and output/.
if [[ -n "${PROJECT:-}" ]]; then
  DATA_DIR="data_${PROJECT}"
  OUTPUT_DIR="output_${PROJECT}"
else
  DATA_DIR="data"
  OUTPUT_DIR="output"
fi

# Build with Homebrew GCC on macOS (libc++ cannot handle std::execution::par
# and std::chrono::zoned_time used by the project).
if [[ -x "/opt/homebrew/bin/g++-14" ]]; then
  make CXX=/opt/homebrew/bin/g++-14 \
    CXXFLAGS="-O3 -std=c++20 -I. -I/opt/homebrew/include" \
    LDLIBS="-L/opt/homebrew/lib -Wl,-rpath,/opt/homebrew/lib -lflint -lgmp -lmimalloc -ltbb"
else
  make
fi

BOOTSTRAP="./bootstrap"
if [[ -x "./bootstrap.exe" ]]; then
  BOOTSTRAP="./bootstrap.exe"
fi

TARGET="${1:-SEW_5p1}"
mkdir -p "${OUTPUT_DIR}" logs

# Run symmetry solving for all three symmetry types (cyclic, flip, parity).
# Collinear has its own solver (--solve-collinear, see run_collinear.sh).
for SYMMETRY in cyclic flip parity; do
  LOG="logs/solve_symmetry_${SYMMETRY}_${TARGET}.log"
  echo "==> ${SYMMETRY} (${LOG})"
  "$BOOTSTRAP" --solve-symmetry --symmetry "$SYMMETRY" --target "$TARGET" \
    --data-dir "$DATA_DIR" --output-dir "$OUTPUT_DIR" 2>&1 | tee "$LOG"
done

echo ""
echo "=== All symmetry solving complete ==="
echo "Data dir:   ${DATA_DIR}"
echo "Output dir: ${OUTPUT_DIR}"
echo "Logs:       logs/solve_symmetry_*_${TARGET}.log"
echo "Invariants: ${OUTPUT_DIR}/<symmetry>/${TARGET}_invariant.wxf"
