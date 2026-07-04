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

# Run symmetry solving for all three symmetry types (cyclic, flip, parity).
# Collinear has its own solver (solve_collinear, to be developed).
for SYMMETRY in cyclic flip parity; do
  LOG="logs/solve_symmetry_${SYMMETRY}_${TARGET}.log"
  echo "==> ${SYMMETRY} (${LOG})"
  "$BOOTSTRAP" --solve-symmetry --symmetry "$SYMMETRY" --target "$TARGET" 2>&1 | tee "$LOG"
done

echo ""
echo "=== All symmetry solving complete ==="
echo "Logs:       logs/solve_symmetry_*_${TARGET}.log"
echo "Invariants: output/<symmetry>/${TARGET}_invariant.wxf"
