#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 || $# -gt 5 ]]; then
  echo "usage: $0 TARGET_BACKEND_SOURCE /path/to/pqm4 [PORTABLE_HOST_BACKEND_SOURCE] [LEVEL] [BACKEND_BINDINGS_C]" >&2
  exit 2
fi
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
TARGET_SRC=$(cd -- "$1" && pwd)
PQM4=$(cd -- "$2" && pwd)
HOST_SRC="$TARGET_SRC"
if [[ $# -ge 3 && -n ${3:-} ]]; then HOST_SRC=$(cd -- "$3" && pwd); fi
LEVEL=${4:-${SQISIGN_LEVEL:-1}}
BIND=${5:-$ROOT/pqm4-overlay/backend_bindings.c}
case "$LEVEL" in 1|3|5) ;; *) echo "unsupported level: $LEVEL" >&2; exit 2;; esac
[[ -f "$BIND" ]] || { echo "binding file not found: $BIND" >&2; exit 1; }
TARGET="$PQM4/crypto_sign/sqisign-qlapoti-lvl$LEVEL/m4"
HOST="$PQM4/mupq/crypto_sign/sqisign-qlapoti-lvl$LEVEL/ref"
[[ -d "$TARGET" && -d "$HOST" ]] || {
  echo "run install_into_pqm4.sh for Level $LEVEL first" >&2; exit 1;
}
"$ROOT/scripts/inject_backend.sh" "$TARGET_SRC" "$TARGET" "$BIND"
"$ROOT/scripts/inject_backend.sh" "$HOST_SRC" "$HOST" "$BIND"
echo "Imported Level $LEVEL target and portable host backend trees."
echo "Installed and validated the three backend contract providers in both implementations."
