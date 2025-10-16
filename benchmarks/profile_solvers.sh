#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
BUILD_DIR="$REPO_ROOT/bin"
PROF_DIR="$REPO_ROOT/benchmarks/profiles"

mkdir -p "$PROF_DIR"
cd "$REPO_ROOT"

if [[ ! -x "$BUILD_DIR/ns_fd_solver" || ! -x "$BUILD_DIR/ns_spectral_solver" ]]; then
  echo "[profile] Building solvers with profiling flags (-pg)..."
  make clean
  make performance PROFILE_FLAGS="-pg"
fi

run_gprof() {
  local binary=$1
  local output=$2
  echo "[profile] Running $binary under gprof..."
  GMON_OUT="gmon.out" "$binary" >/dev/null 2>&1 || true
  if [[ -f gmon.out ]]; then
    gprof "$binary" gmon.out > "$output"
    mv gmon.out "$PROF_DIR/${binary##*/}.gmon"
  else
    echo "[profile] Warning: gmon.out not generated for $binary" >&2
  fi
}

run_gprof "$BUILD_DIR/ns_fd_solver" "$PROF_DIR/fd_gprof.txt"
run_gprof "$BUILD_DIR/ns_spectral_solver" "$PROF_DIR/spectral_gprof.txt"

generate_flamegraph() {
  local binary=$1
  local svg=$2
  local flamegraph_exec="${FLAMEGRAPH:-FlameGraph/flamegraph.pl}"
  if command -v perf >/dev/null 2>&1 && [[ -x "$flamegraph_exec" || "$(command -v flamegraph.pl 2>/dev/null)" != "" ]]; then
    echo "[profile] Capturing perf data for $binary..."
    perf record -F 999 -g "$binary" >/dev/null 2>&1 || true
    if [[ ! -x "$flamegraph_exec" ]]; then
      flamegraph_exec="$(command -v flamegraph.pl)"
    fi
    perf script | "$flamegraph_exec" > "$svg"
    mv perf.data "$PROF_DIR/${binary##*/}.perf"
  else
    echo "[profile] Skipping flame graph for $binary (perf or FlameGraph missing)."
  fi
}

generate_flamegraph "$BUILD_DIR/ns_fd_solver" "$PROF_DIR/fd_flamegraph.svg"
generate_flamegraph "$BUILD_DIR/ns_spectral_solver" "$PROF_DIR/spectral_flamegraph.svg"

echo "[profile] Profiles stored in $PROF_DIR"
echo "[profile] Review gprof outputs for hotspot analysis and update optimisation notes accordingly."
