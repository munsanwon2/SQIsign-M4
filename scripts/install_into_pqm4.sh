#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "usage: $0 /path/to/pqm4 [1,3,5|LEVEL]" >&2
  exit 2
fi
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
PQM4=$(cd -- "$1" && pwd)
case "$PQM4" in *[[:space:]]*) echo "pqm4 path containing whitespace is unsupported" >&2; exit 1;; esac
[[ -f "$PQM4/Makefile" && -d "$PQM4/crypto_sign" && -d "$PQM4/mupq/crypto_sign" && -d "$PQM4/ldscripts" ]] || {
  echo "not a recursive pqm4 checkout (mupq submodule or ldscripts missing?): $PQM4" >&2; exit 1;
}
selector=${2:-1,3,5}
selector=${selector//,/ }
LEVELS=()
for level in $selector; do
  case "$level" in 1|3|5) LEVELS+=("$level");; *) echo "unsupported level: $level" >&2; exit 2;; esac
done
[[ ${#LEVELS[@]} -gt 0 ]] || { echo "no levels selected" >&2; exit 2; }

LDS_NAMES=(stm32f4discovery_sqisign_qlapoti.ld stm32f4discovery_fullram_sqisign_qlapoti.ld)
created_dirs=()
created_ld=()
for level in "${LEVELS[@]}"; do
  target="$PQM4/crypto_sign/sqisign-qlapoti-lvl$level"
  host="$PQM4/mupq/crypto_sign/sqisign-qlapoti-lvl$level"
  for path in "$target" "$host"; do
    [[ ! -e "$path" ]] || { echo "refusing to overwrite $path" >&2; exit 1; }
  done
  created_dirs+=("$target" "$host")
done
# Shared linker profiles may already exist from a previously installed level,
# but only byte-identical port-kit copies are accepted.
for name in "${LDS_NAMES[@]}"; do
  src="$ROOT/pqm4-overlay/ldscripts/$name"
  dst="$PQM4/ldscripts/$name"
  if [[ -e "$dst" ]]; then
    cmp -s "$src" "$dst" || { echo "refusing conflicting linker profile: $dst" >&2; exit 1; }
  fi
done
cleanup() {
  for path in "${created_dirs[@]}"; do rm -rf "$path"; done
  for path in "${created_ld[@]}"; do rm -f "$path"; done
}
trap cleanup ERR INT TERM
for level in "${LEVELS[@]}"; do
  target="$PQM4/crypto_sign/sqisign-qlapoti-lvl$level"
  host="$PQM4/mupq/crypto_sign/sqisign-qlapoti-lvl$level"
  cp -a "$ROOT/pqm4-overlay/crypto_sign/sqisign-qlapoti-lvl$level" "$target"
  cp -a "$ROOT/pqm4-overlay/mupq/crypto_sign/sqisign-qlapoti-lvl$level" "$host"
  printf 'Installed Level %s target: %s\n' "$level" "$target/m4"
  printf 'Installed Level %s host ref: %s\n' "$level" "$host/ref"
done
for name in "${LDS_NAMES[@]}"; do
  src="$ROOT/pqm4-overlay/ldscripts/$name"
  dst="$PQM4/ldscripts/$name"
  if [[ ! -e "$dst" ]]; then
    cp "$src" "$dst"
    created_ld+=("$dst")
  fi
done
trap - ERR INT TERM
printf 'Linker profiles available under %s/ldscripts\n' "$PQM4"
echo "No existing pqm4 Makefile, mk/config, platform file, linker script, or submodule metadata was modified."
echo "Expected state: each final cryptographic link fails until its genuine level-specific backend is imported."
