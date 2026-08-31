#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "usage: $0 BACKEND_SOURCE IMPLEMENTATION_DIR [BACKEND_BINDINGS_C]" >&2
  exit 2
fi
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
SRC=$(cd -- "$1" && pwd)
IMPL=$(cd -- "$2" && pwd)
BIND=${3:-}
# The source is copied into the implementation and is safe when quoted, so a
# workspace path containing spaces must remain usable.  GNU make cannot safely
# consume an implementation path containing whitespace because MAKEFILE_LIST
# and source lists are whitespace-delimited; retain that narrower restriction.
case "$IMPL" in *[[:space:]]*) echo "pqm4 implementation paths containing whitespace are unsupported" >&2; exit 1;; esac
[[ -f "$IMPL/api.h" && -f "$IMPL/config.mk" && -f "$IMPL/sqisign_m4_backend.h" ]] || {
  echo "destination is not an installed sqisign-qlapoti implementation" >&2; exit 1;
}
# A curated export must be passed at the level root so its manifest controls
# hashing, active TARGET_ARM audit, and generated compile flags.  Treating the
# source/ payload as a legacy tree would silently bypass all three gates.
if [[ $(basename -- "$SRC") == source && -f "$SRC/../library-manifest.json" ]]; then
  echo "curated source/ payload passed directly; pass its parent export root containing library-manifest.json" >&2
  exit 1
fi
CURATED=0
PAYLOAD="$SRC"
if [[ -f "$SRC/library-manifest.json" ]]; then
  [[ -d "$SRC/source" ]] || {
    echo "curated backend export has library-manifest.json but no source/ directory" >&2; exit 1;
  }
  CURATED=1
  PAYLOAD=$(cd -- "$SRC/source" && pwd)
fi
if find "$SRC" -type l -print -quit | grep -q .; then
  echo "materialize backend symlinks before import" >&2; exit 1
fi
if find "$SRC" -type f -printf '%P\n' | grep -Eq '[[:space:]]'; then
  echo "backend file names containing whitespace are unsupported by pqm4 make paths" >&2; exit 1
fi
assert_no_backend_artifacts() {
  local x
  if find "$IMPL" -mindepth 1 -maxdepth 1 -name 'backend_tu_*' \
      -print -quit | grep -q .; then
    echo "refusing pre-existing backend wrapper artifact under: $IMPL" >&2
    return 1
  fi
  for x in "$IMPL/backend" "$IMPL/BACKEND_TRANSLATION_UNITS.txt" \
           "$IMPL/BACKEND_LIBRARY_MANIFEST.json" "$IMPL/backend_compile_flags.mk" \
           "$IMPL/backend_bindings.c" "$IMPL/backend_bindings.c.todo"; do
    [[ ! -e "$x" ]] || {
      echo "refusing to overwrite existing backend artifact: $x" >&2; return 1;
    }
  done
}
audit_backend_tree() {
  local tree=$1
  local manifest=${2:-}
  local args=("$tree")
  if [[ -n "$manifest" ]]; then
    args+=(--library-manifest "$manifest")
  fi
  python3 "$ROOT/tools/check_no_host_dependencies.py" "${args[@]}"
  if find "$tree" -type f \( -iname 'randombytes.c' -o -iname 'randombytes_*.c' -o -iname 'randombytes.s' -o -iname 'randombytes.S' \) | grep -q .; then
    echo "refusing backend containing a randombytes implementation file" >&2; return 1
  fi
  if grep -RIE --include='*.c' \
    '(^|[^A-Za-z0-9_])int[[:space:]]+crypto_sign_(keypair|open)[[:space:]]*\(|(^|[^A-Za-z0-9_])int[[:space:]]+crypto_sign[[:space:]]*\(' \
    "$tree" >/dev/null; then
    echo "backend contains a competing NIST API wrapper; prepare a library-only tree" >&2; return 1
  fi
  if grep -RIE --include='*.c' '(^|[^A-Za-z0-9_])int[[:space:]]+main[[:space:]]*\(' "$tree" >/dev/null; then
    echo "backend contains application/test main(); prepare a library-only tree" >&2; return 1
  fi
}
assert_no_backend_artifacts
if [[ "$CURATED" == 1 ]]; then
  audit_backend_tree "$PAYLOAD" "$SRC/library-manifest.json"
else
  audit_backend_tree "$PAYLOAD"
fi

STAGE=$(mktemp -d "$IMPL/.backend-import.XXXXXX")
trap 'rm -rf "$STAGE"' EXIT
if [[ "$CURATED" == 1 ]]; then
  python3 "$ROOT/full-backend-tools/consume_pqm4_manifest.py" \
    --export-root "$SRC" --implementation "$IMPL" --stage "$STAGE"
else
  mkdir -p "$STAGE/backend"
  tar -C "$PAYLOAD" --exclude=.git --exclude=build --exclude='*.o' --exclude='*.a' \
    --exclude='*.so' --exclude='*.dylib' --exclude='__pycache__' \
    --exclude='PQCgenKAT*' --exclude='test' --exclude='tests' \
    -cf - . | tar -C "$STAGE/backend" -xf -
fi
# Audit the immutable staged snapshot as well as the input tree.  This closes
# the check/copy race without weakening legacy imports.
if [[ "$CURATED" == 1 ]]; then
  audit_backend_tree "$STAGE/backend" "$STAGE/BACKEND_LIBRARY_MANIFEST.json"
else
  audit_backend_tree "$STAGE/backend"
fi

python3 "$ROOT/tools/render_pqm4_backend_wrappers.py" \
  --backend "$STAGE/backend" --destination "$STAGE"

if [[ -n "$BIND" ]]; then
  [[ -f "$BIND" ]] || { echo "binding file not found: $BIND" >&2; exit 1; }
  cp "$BIND" "$STAGE/backend_bindings.c"
fi
count=$(python3 - "$STAGE" <<'PY3'
from collections import Counter
from pathlib import Path
import re, sys
stage = Path(sys.argv[1])
rx = re.compile(r'^[ \t]*int[ \t\r\n]+sqisign_m4_backend_(keypair|sign|verify)[ \t]*\([^;]*?\)[ \t\r\n]*\{', re.M | re.S)
found = []
paths = list((stage / 'backend').rglob('*.c'))
binding = stage / 'backend_bindings.c'
if binding.is_file():
    paths.append(binding)
for path in paths:
    text = path.read_text(encoding='utf-8', errors='replace')
    found.extend(rx.findall(text))
if not found:
    print(0)
else:
    counts = Counter(found)
    expected = Counter({'keypair': 1, 'sign': 1, 'verify': 1})
    if counts != expected:
        print('backend contract definitions must contain exactly one each of '
              'keypair, sign, and verify; found ' + repr(dict(counts)),
              file=sys.stderr)
        raise SystemExit(1)
    print(3)
PY3
)
if [[ "$count" == 0 ]]; then
  cp "$ROOT/templates/backend_bindings.c.in" "$STAGE/backend_bindings.c.todo"
fi
if [[ "$CURATED" == 1 ]]; then
  expected_config_hash=$(tr -d '\r\n' < "$STAGE/BACKEND_CONFIG_ORIGINAL.sha256")
  actual_config_hash=$(sha256sum "$IMPL/config.mk")
  actual_config_hash=${actual_config_hash%% *}
  [[ "$actual_config_hash" == "$expected_config_hash" ]] || {
    echo "config.mk changed while curated backend was staged; refusing partial import" >&2; exit 1;
  }
fi
assert_no_backend_artifacts
mv "$STAGE/backend" "$IMPL/backend"
for f in "$STAGE"/backend_tu_* "$STAGE/BACKEND_TRANSLATION_UNITS.txt" \
         "$STAGE/BACKEND_LIBRARY_MANIFEST.json" "$STAGE/backend_compile_flags.mk" \
         "$STAGE/backend_bindings.c" "$STAGE/backend_bindings.c.todo"; do
  [[ -e "$f" ]] && mv "$f" "$IMPL/"
done
if [[ "$CURATED" == 1 ]]; then
  mv "$STAGE/config.mk" "$IMPL/config.mk"
  rm "$STAGE/BACKEND_CONFIG_ORIGINAL.sha256"
fi
trap - EXIT
rm -rf "$STAGE"
if [[ "$count" == 0 ]]; then
  echo "Backend imported; complete backend_bindings.c.todo and rename it backend_bindings.c."
else
  echo "Backend imported with exactly three backend contract definitions."
fi
if [[ "$CURATED" == 1 ]]; then
  echo "Curated manifest definitions, options, and include directories are active via backend_compile_flags.mk."
fi
