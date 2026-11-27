"""Run inside GDB to record hal_send_str calls through a local st-util server."""

import os
import gdb

gdb.execute("set pagination off")
gdb.execute("set confirm off")
port = int(os.environ.get("SQISIGN_SWD_PORT", "42425"))
gdb.execute(f"target remote 127.0.0.1:{port}")
gdb.execute("monitor reset")

transcript = open(os.environ["SQI_TRANSCRIPT"], "x", encoding="utf-8", buffering=1)
finished = False

class Capture(gdb.Breakpoint):
    def stop(self):
        global finished
        value = gdb.parse_and_eval("(char *)$r0")
        line = value.string(encoding="ascii", errors="replace")
        transcript.write(line + "\n")
        finished = line == "#"
        return finished

Capture("*hal_send_str", type=gdb.BP_HARDWARE_BREAKPOINT, internal=True)
gdb.Breakpoint("hard_fault_handler", type=gdb.BP_HARDWARE_BREAKPOINT, internal=True)
gdb.execute("continue")
transcript.close()
if not finished:
    gdb.execute("info registers pc sp lr")
    gdb.execute("x/3wx 0xE000ED28")
    gdb.execute("quit 1")
gdb.execute("disconnect")
