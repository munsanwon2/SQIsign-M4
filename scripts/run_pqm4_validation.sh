#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "usage: $0 /path/to/pqm4 PLATFORM [1|3|5|all]" >&2
  exit 2
fi
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
PQM4=$(cd -- "$1" && pwd)
PLATFORM=$2
SELECT=${3:-all}
case "$SELECT" in
  all) LEVELS=(1 3 5);;
  1|3|5) LEVELS=("$SELECT");;
  *) echo "unsupported level selector: $SELECT" >&2; exit 2;;
esac
if [[ "$PLATFORM" == stm32f4discovery ]]; then
  [[ -f "$PQM4/ldscripts/stm32f4discovery_sqisign_qlapoti.ld" ]] || {
    echo "missing SQIsign-specific F4 linker profile; rerun install_into_pqm4.sh" >&2; exit 1;
  }
fi
cd "$PQM4"
# Clean once so an all-level validation retains every Level-I/III/V ELF and
# host executable together for final hashing and cross-level audit.
make PLATFORM="$PLATFORM" clean
for LEVEL in "${LEVELS[@]}"; do
  SCHEME="sqisign-qlapoti-lvl$LEVEL"
  TARGET="$PQM4/crypto_sign/$SCHEME/m4"
  HOST="$PQM4/mupq/crypto_sign/$SCHEME/ref"
  for impl in "$TARGET" "$HOST"; do
    python3 "$ROOT/tools/check_pqm4_layout.py" --require-backend "$impl"
    python3 "$ROOT/tools/check_no_host_dependencies.py" \
      "$impl/backend" --library-manifest "$impl/BACKEND_LIBRARY_MANIFEST.json"
    python3 "$ROOT/tools/audit_no_mock_api.py" --allow-backend-definitions "$impl"
  done
  make -j2 PLATFORM="$PLATFORM" IMPLEMENTATION_PATH="crypto_sign/$SCHEME/m4" \
    "elf/crypto_sign_${SCHEME}_m4_test.elf" \
    "elf/crypto_sign_${SCHEME}_m4_speed.elf" \
    "elf/crypto_sign_${SCHEME}_m4_hashing.elf" \
    "elf/crypto_sign_${SCHEME}_m4_stack.elf" \
    "elf/crypto_sign_${SCHEME}_m4_testvectors.elf"
  make -j2 PLATFORM="$PLATFORM" IMPLEMENTATION_PATH="mupq/crypto_sign/$SCHEME/ref" \
    "bin-host/mupq_crypto_sign_${SCHEME}_ref_testvectors"
  cat <<EOF2
Level $LEVEL build gates passed. Mandatory runtime/measurement commands remain:
  python3 test.py --platform $PLATFORM --opt size $SCHEME
  python3 testvectors.py --platform $PLATFORM --opt size $SCHEME
  python3 benchmarks.py --platform $PLATFORM --opt size $SCHEME
Use the board measurements, ELF/map files, and stack canary results to confirm
the required margin for the intended interrupt configuration before deployment.
EOF2
done
