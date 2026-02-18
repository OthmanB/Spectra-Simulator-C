#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT_DIR/build}"

echo "Building specsim (CMake)..."
cmake -S "$ROOT_DIR" -B "$BUILD_DIR" >/dev/null
cmake --build "$BUILD_DIR" -j >/dev/null

echo "Running smoke tests..."
python3 "$ROOT_DIR/scripts/smoke_test.py" \
  --specsim "$BUILD_DIR/specsim" \
  --examples-dir "$ROOT_DIR/Configurations/examples_cfg" \
  --noise-cfg "$ROOT_DIR/Configurations/noise_Kallinger2014.cfg"
