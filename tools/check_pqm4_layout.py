#!/usr/bin/env python3
from pathlib import Path
import argparse, re, sys

REQUIRED={
 'api.h','pqm4_api.c','sqisign_m4_backend.h','qlapoti_workspace.h',
 'qlapoti_workspace.c','m4_operation_guard.h','m4_operation_guard.c',
 'malloc_compat.c','secure_bzero.h','secure_bzero.c','config.mk','LICENSE'
}
CONTRACT=re.compile(r'^\s*int\s+sqisign_m4_backend_(keypair|sign|verify)\s*\([^;]*\)\s*\{',re.M|re.S)

def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--require-backend',action='store_true')
    ap.add_argument('implementation',type=Path)
    a=ap.parse_args(); p=a.implementation.resolve(); errors=[]
    if not p.is_dir():
        errors.append(f'not a directory: {p}')
    else:
        if p.name not in {'m4','ref'}: errors.append(f'expected m4 or ref, got {p.name}')
        if p.name=='m4' and not (p.parent.parent.name=='crypto_sign' and '/mupq/crypto_sign/' not in ('/'+p.as_posix().lstrip('/'))):
            errors.append('m4 must be under pqm4/crypto_sign/<scheme>/m4')
        if p.name=='ref' and not (p.parent.parent.name=='crypto_sign' and '/mupq/crypto_sign/' in ('/'+p.as_posix().lstrip('/'))):
            errors.append('ref must be under pqm4/mupq/crypto_sign/<scheme>/ref')
        names={x.name for x in p.iterdir()}
        for name in sorted(REQUIRED-names): errors.append(f'missing {name}')
        for f in p.iterdir():
            if f.is_symlink(): errors.append(f'top-level source symlink forbidden; use wrapper TU: {f.name}')
        for f in p.rglob('*'):
            if f.is_file() and f.name.lower() in {'randombytes.c','randombytes_system.c','randombytes.s'}:
                errors.append(f'do not bundle randombytes implementation: {f}')
        if a.require_backend:
            if not (p/'backend').is_dir(): errors.append('backend/ missing')
            if not (p/'BACKEND_TRANSLATION_UNITS.txt').is_file(): errors.append('backend TU manifest missing')
            if (p/'backend_bindings.c.todo').exists(): errors.append('backend_bindings.c.todo remains unresolved')
            defs={}
            for f in p.rglob('*.c'):
                text=f.read_text(errors='replace')
                for n in CONTRACT.findall(text): defs.setdefault(n,[]).append(str(f))
            for n in ('keypair','sign','verify'):
                if len(defs.get(n,[]))!=1:
                    errors.append(f'backend contract {n}: expected 1 definition, got {len(defs.get(n,[]))}')
    if errors:
        for e in errors: print('ERROR:',e,file=sys.stderr)
        return 1
    print('PQM4 LAYOUT'+(' + GENUINE-BACKEND CONTRACT' if a.require_backend else '')+': PASS')
    return 0
if __name__=='__main__': raise SystemExit(main())
