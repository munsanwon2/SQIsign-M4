#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

for forbidden in \
  work docs reports tmp upstreams staging .agents .codex \
  README_FIRST.md SQIsign_Qlapoti_CortexM4_SESSION_HANDOFF_20260822.md \
  STATUS.md FULL_BACKEND_STATUS.md; do
  [[ ! -e "$ROOT/$forbidden" ]] || {
    echo "forbidden release path is present: $forbidden" >&2
    exit 1
  }
done

if find "$ROOT" -type f \( -name '*SESSION_HANDOFF*.md' -o -name '*AUDIT*.md' \) \
    -print -quit | grep -q .; then
  echo "work-memory or audit report found in release tree" >&2
  exit 1
fi

for script in scripts/*.sh; do
  bash -n "$script"
done

python3 - <<'PY'
import ast
from pathlib import Path

for path in sorted(Path('.').rglob('*.py')):
    ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
print('PYTHON SYNTAX: PASS')
PY

for level in 1 3 5; do
  export_root="$ROOT/backends/lvl$level"
  [[ -f "$export_root/library-manifest.json" && -d "$export_root/source" ]] || {
    echo "missing curated Level-$level backend" >&2
    exit 1
  }
  python3 "$ROOT/tools/check_no_host_dependencies.py" \
    "$export_root/source" --library-manifest "$export_root/library-manifest.json"
done

echo "RELEASE TREE AUDIT: PASS"
