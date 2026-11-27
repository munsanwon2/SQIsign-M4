#!/usr/bin/env python3
"""Run wire-format checks against the installed pqm4 host audit archives."""

import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def field_vectors(prime, size):
    rng = random.Random(723)
    operands = [(prime - 1, prime - 2, 1, 0), (1, 0, prime - 1, prime - 1)]
    operands += [tuple(rng.randrange(1, prime) for _ in range(4)) for _ in range(4)]

    def c_bytes(*components):
        data = b"".join((x % prime).to_bytes(size, "little") for x in components)
        return '"' + "".join(f"\\x{x:02x}" for x in data) + '"'

    cases = []
    for a, b, c, d in operands:
        inverse_norm = pow((a * a + b * b) % prime, -1, prime)
        fields = [c_bytes(a, b), c_bytes(c, d), c_bytes(a + c), c_bytes(a - c),
                  c_bytes(a * c), c_bytes(pow(a, -1, prime)), c_bytes(a * a),
                  c_bytes(a * c - b * d, a * d + b * c),
                  c_bytes(a * inverse_norm, -b * inverse_norm), c_bytes(a * a - b * b, 2 * a * b)]
        cases.append("    {\n        " + ",\n        ".join(fields) + "\n    },\n")
    return "static const struct field_vector field_vectors[] = {\n" + "".join(cases) + "};\n"


def run(pqm4):
    source = Path(__file__).with_name("encoding_regression.c").resolve()
    with tempfile.TemporaryDirectory(prefix="sqisign-encoding-") as temporary:
        for level, cofactor, exponent in ((1, 3, 324), (3, 27, 500), (5, 17, 664)):
            scheme = f"sqisign-qlapoti-lvl{level}"
            path = pqm4 / f"mupq/crypto_sign/{scheme}/ref"
            archive = pqm4 / f"obj-host/libmupq_crypto_sign_{scheme}_ref_vla_audit.a"
            manifest = json.loads((path / "BACKEND_LIBRARY_MANIFEST.json").read_bytes())
            backend = path / "backend"
            defines = [f"-D{item}" for item in manifest["pqm4_compile_definitions"]
                       if not item.startswith(("TARGET_", "SQISIGN_TARGET_OS_"))]
            includes = [f"-I{backend / directory}" for directory in manifest["include_directories"]]
            output = Path(temporary) / f"encoding-lvl{level}"
            prime = (cofactor << exponent) - 1
            wire_bytes = {1: 48, 3: 64, 5: 88}[level]
            (Path(temporary) / "field_vectors.h").write_text(field_vectors(prime, wire_bytes), encoding="ascii")
            subprocess.run([os.environ.get("CC", "cc"), "-std=c11", "-O2",
                            "-Werror=implicit-function-declaration", "-ffunction-sections", "-fdata-sections",
                            *defines, *includes,
                            f"-I{temporary}",
                            f"-I{backend / 'src/signature/ref/lvlx'}", f"-I{backend / 'src/hd/ref/lvlx'}",
                            str(source), str(archive),
                            "-Wl,--gc-sections", "-lm", "-o", str(output)], check=True)
            print(f"Level {level}", flush=True)
            subprocess.run([str(output)], check=True)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} /path/to/built/pqm4")
    run(Path(sys.argv[1]).resolve())
