#!/usr/bin/env python3
"""Reject fixed-vector/mock paths in production implementation sources."""
from pathlib import Path
import argparse, re, sys

TEXT_SUFFIX={'.c','.h','.S','.s','.inc','.mk','.txt'}
PATTERNS=[
 (re.compile(r'\b(mock|dummy|fake)[-_ ]?(key|sign|signature|backend)',re.I),'mock/dummy/fake cryptographic path'),
 (re.compile(r'fixed[-_ ]?(key|signature)',re.I),'fixed key/signature marker'),
 (re.compile(r'PQCLEAN.*FIXED',re.I),'fixed test-vector marker'),
 (re.compile(r'static\s+(?:const\s+)?(?:unsigned\s+char|uint8_t)\s+[^;]*(?:SECRETKEYBYTES|CRYPTO_BYTES)[^;]*=\s*\{',re.S),'embedded key/signature-sized array'),
]
DEF=re.compile(r'\bint\s+sqisign_m4_backend_(?:keypair|sign|verify)\s*\([^;]*\)\s*\{',re.S)


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--allow-backend-definitions',action='store_true')
    ap.add_argument('root',type=Path)
    a=ap.parse_args(); errors=[]; definitions=0
    for p in a.root.rglob('*'):
        if not p.is_file() or p.suffix not in TEXT_SUFFIX:
            continue
        raw=p.read_text(errors='replace')
        text=re.sub(r'/\*.*?\*/','',raw,flags=re.S)
        text=re.sub(r'//.*','',text)
        # Backend source may legitimately mention test terminology in comments,
        # but comments were stripped. Test/application subtrees must not be imported.
        for rx,label in PATTERNS:
            if rx.search(text): errors.append(f'{p}: {label}')
        if p.suffix=='.c': definitions += len(DEF.findall(text))
    if not a.allow_backend_definitions and definitions:
        errors.append(f'production scaffold defines {definitions} backend function(s); genuine backend must be injected separately')
    if a.allow_backend_definitions and definitions!=3:
        errors.append(f'imported implementation defines {definitions} backend functions; expected exactly 3')
    if errors:
        for e in errors: print('ERROR:',e,file=sys.stderr)
        return 1
    print('NO MOCK/FIXED-VECTOR CRYPTOGRAPHIC PATH: PASS')
    return 0
if __name__=='__main__': raise SystemExit(main())
