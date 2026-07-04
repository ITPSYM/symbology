#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

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
mkdir -p output logs

# Run all four symmetries, each with its own log file.
for SYMMETRY in collinear cyclic flip parity; do
  LOG="logs/project_${SYMMETRY}_${TARGET}.log"
  echo "==> ${SYMMETRY} (${LOG})"
  "$BOOTSTRAP" --project --symmetry "$SYMMETRY" --target "$TARGET" 2>&1 | tee "$LOG"
done

echo ""
echo "=== All projections complete ==="
echo "Logs:       logs/project_*_${TARGET}.log"
echo "Summaries:  output/<symmetry>/summary.txt"
