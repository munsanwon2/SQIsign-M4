#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "usage: $0 /path/to/pqm4 [1,3,5|LEVEL]" >&2
  exit 2
fi

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
PQM4=$(cd -- "$1" && pwd)
SELECTOR=${2:-1,3,5}
SELECTOR=${SELECTOR//,/ }
LEVELS=()
for level in $SELECTOR; do
  case "$level" in
    1|3|5) LEVELS+=("$level");;
    *) echo "unsupported level: $level" >&2; exit 2;;
  esac
done
[[ ${#LEVELS[@]} -gt 0 ]] || { echo "no levels selected" >&2; exit 2; }

CSV=$(IFS=,; echo "${LEVELS[*]}")
"$ROOT/scripts/install_into_pqm4.sh" "$PQM4" "$CSV"

for level in "${LEVELS[@]}"; do
  export_root="$ROOT/backends/lvl$level"
  [[ -f "$export_root/library-manifest.json" && -d "$export_root/source" ]] || {
    echo "missing curated Level-$level backend: $export_root" >&2
    exit 1
  }
  "$ROOT/scripts/inject_scheme_backends.sh" \
    "$export_root" "$PQM4" "$export_root" "$level"
done

echo "Installed the selected self-contained SQIsign-Qlapoti target and host implementations."
