#!/usr/bin/env bash
set -euo pipefail

# Opt-in host regressions.  This script is never called by the target build
# or by verify_release.sh.  CC must name a GCC/Clang executable with GNU ld
# compatible --wrap and --gc-sections support (for example, on Linux/WSL).
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
BUILD_DIR=$(mktemp -d "${TMPDIR:-/tmp}/sqisign-quaternion.XXXXXX")
trap 'rm -rf -- "$BUILD_DIR"' EXIT
CC=${CC:-cc}

levels=(1 3 5)
if (( $# > 1 )); then
  echo "usage: bash tests/quaternion/run.sh [1|3|5]" >&2
  exit 2
fi
if (( $# == 1 )); then
  case "$1" in
    1|3|5) levels=("$1") ;;
    *) echo "level must be 1, 3, or 5" >&2; exit 2 ;;
  esac
fi

for level in "${levels[@]}"; do
  source_root="$ROOT/backends/lvl$level/source"
  quaternion="$source_root/src/quaternion/ref/generic"
  flags=(
    -std=c11 -O1 -g -Wall -Wextra
    -ffunction-sections -fdata-sections -fno-strict-aliasing
    -DDISABLE_NAMESPACING -DRADIX_32 -DSQISIGN_BUILD_TYPE_REF
    "-DSQISIGN_VARIANT=lvl$level" -DSQISIGN_INTBIG_OVERFLOW_CHECK
    "-I$source_root/include" "-I$quaternion/include"
    "-I$quaternion/internal_quaternion_headers" "-I$quaternion"
    "-I$quaternion/lll"
  )
  common=("$quaternion/intbig.c" "$quaternion/finit.c"
          "$quaternion/dim2.c" "$quaternion/dim4.c" "$quaternion/algebra.c")

  # UBSan catches the original signed subtraction overflow in normeq.c.
  "$CC" "${flags[@]}" -fsanitize=undefined -fno-sanitize-recover=all \
    "$ROOT/tests/quaternion/normeq_regression.c" \
    "$quaternion/normeq.c" "$quaternion/lattice.c" "${common[@]}" \
    -Wl,--wrap=ibz_rand_interval -Wl,--wrap=ibz_probab_prime \
    -Wl,--wrap=ibz_cornacchia_prime -Wl,--gc-sections -lm \
    -o "$BUILD_DIR/normeq-lvl$level"
  "$BUILD_DIR/normeq-lvl$level"

  # ml2_regression.c includes the actual ml2.c with a checked DPE divisor.
  # Keep this independent of the separate DPE zero-operand/UBSan repair.
  "$CC" "${flags[@]}" "$ROOT/tests/quaternion/ml2_regression.c" \
    "$quaternion/lll/mlll.c" "$quaternion/dpe.c" "${common[@]}" \
    -Wl,--gc-sections -lm -o "$BUILD_DIR/ml2-lvl$level"
  "$BUILD_DIR/ml2-lvl$level"
  echo "Level $level quaternion regressions passed"
done
