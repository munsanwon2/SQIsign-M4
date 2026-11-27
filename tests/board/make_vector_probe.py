#!/usr/bin/env python3
"""Instrument the first vector of the pinned pqm4 signature harness."""

from pathlib import Path
import sys

if len(sys.argv) != 2:
    sys.exit(f"Usage: {sys.argv[0]} /path/to/pinned/pqm4")
pqm4 = Path(sys.argv[1])
source = pqm4 / "mupq/crypto_sign/testvectors.c"
probe = pqm4 / "mupq/crypto_sign/measure_testvectors.c"
text = source.read_text(encoding="utf-8")
if "#define MAXMLEN 2048" not in text:
    sys.exit("Unexpected pqm4 testvector harness")
text = text.replace("#define MAXMLEN 2048", "#define MAXMLEN 1", 1)
text = text.replace("#include <stdint.h>", '''#include <stdint.h>
#include "m4_operation_guard.h"

extern char end, _stack;
static int report(const char *name, int rc, size_t peak, size_t available) {
    const qlw_stats_t arena = sqisign_m4_last_workspace_stats;
    char line[220];
    snprintf(line, sizeof(line), "%s status=%d stack=%u stack_bank=%u stack_available=%u capacity=%u peak=%u current=%u live=%u arena=%u",
             name, rc, (unsigned)peak, (unsigned)((uintptr_t)&_stack - (uintptr_t)&end),
             (unsigned)available, (unsigned)arena.capacity_bytes, (unsigned)arena.peak_bytes,
             (unsigned)arena.current_bytes, (unsigned)arena.live_blocks, arena.status_flags);
    hal_send_str(line);
    return peak < available && arena.capacity_bytes > 0 &&
           arena.capacity_bytes <= QLW_WORKSPACE_BYTES && arena.peak_bytes <= arena.capacity_bytes &&
           arena.current_bytes == 0 && arena.live_blocks == 0 && arena.status_flags == 0;
}
''')
text = text.replace("  size_t smlen;", "  size_t smlen, stack_available;")
for label, call in [
    ("keypair", "MUPQ_crypto_sign_keypair(pk, sk)"),
    ("sign", "MUPQ_crypto_sign(sm, &smlen, mi, i, sk)"),
    ("verify", "MUPQ_crypto_sign_open(sm, &mlen, sm, smlen, pk)"),
]:
    old = (("r = " if label == "verify" else "") + call + ";")
    assert text.count(old) == 1
    text = text.replace(old, f'''stack_available = hal_get_stack_size();
    hal_spraystack();
    r = {call};
    if (!report("{label}", r, hal_checkstack(), stack_available) || r != 0) {{
      hal_send_str("RESULT FAIL");
      hal_send_str("#");
      return 1;
    }}''')
text = text.replace('  hal_send_str("#");\n  return 0;', '  hal_send_str("RESULT PASS");\n  hal_send_str("#");\n  return 0;')
if probe.exists() and probe.read_text(encoding="utf-8") != text:
    sys.exit(f"Existing probe differs from generated content: {probe}")
probe.write_text(text, encoding="utf-8", newline="\n")
