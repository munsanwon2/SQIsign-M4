#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "usage: $0 /path/to/pqm4 [1,3,5|LEVEL]" >&2
  exit 2
fi
PQM4=$(cd -- "$1" && pwd)
selector=${2:-1,3,5}
selector=${selector//,/ }
LEVELS=()
for level in $selector; do
  case "$level" in 1|3|5) LEVELS+=("$level");; *) echo "unsupported level: $level" >&2; exit 2;; esac
done
for level in "${LEVELS[@]}"; do
  for d in "$PQM4/crypto_sign/sqisign-qlapoti-lvl$level" \
           "$PQM4/mupq/crypto_sign/sqisign-qlapoti-lvl$level"; do
    [[ -f "$d/PORTKIT_MARKER" ]] || { echo "refusing to remove unmarked path: $d" >&2; exit 1; }
  done
done
for level in "${LEVELS[@]}"; do
  rm -rf "$PQM4/crypto_sign/sqisign-qlapoti-lvl$level" \
         "$PQM4/mupq/crypto_sign/sqisign-qlapoti-lvl$level"
done
# The two profiles are shared by all levels. Remove them only when no marked
# level implementation remains.
remaining=0
for level in 1 3 5; do
  [[ -e "$PQM4/crypto_sign/sqisign-qlapoti-lvl$level" ]] && remaining=1
  [[ -e "$PQM4/mupq/crypto_sign/sqisign-qlapoti-lvl$level" ]] && remaining=1
done
if [[ $remaining -eq 0 ]]; then
  rm -f "$PQM4/ldscripts/stm32f4discovery_sqisign_qlapoti.ld" \
        "$PQM4/ldscripts/stm32f4discovery_fullram_sqisign_qlapoti.ld"
fi
echo "Removed selected marked scheme directories; existing pqm4 core files were untouched."
