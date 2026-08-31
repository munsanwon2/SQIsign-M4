# STM32F407 linker profiles

For `PLATFORM=stm32f4discovery`, every SQIsign-Qlapoti target and portable
reference `config.mk` automatically selects
`stm32f4discovery_sqisign_qlapoti.ld` for this scheme's ELF rules. Installing
the overlay copies the new script beside pqm4's existing linker scripts; it
does not overwrite pqm4 core linker files or change unrelated schemes.

The selected split profile preserves pqm4's physical layout:

- 112 KiB main SRAM for ordinary writable state and stack;
- 16 KiB `ram2` for decoded binding state;
- 64 KiB CCM for the 61,440-byte arena.

It adds explicit assertions for the exact arena size and the binding-state
ceiling. The linker MEMORY/section assignments enforce region capacity and
placement. Inspect the final ELF/map file and measure the runtime high-water
mark before deployment. Current L1/L3/L5 binding sizes are 6,024 / 8,992 /
11,664 bytes.

`stm32f4discovery_fullram_sqisign_qlapoti.ld` is the optional diagnostic
variant. It models 128 KiB of main SRAM plus 64 KiB CCM and moves binding state
into main RAM. Select it only deliberately and repeat all map, symbol,
callgraph, and runtime audits. It is not the default F407 scheme layout and
does not establish support for a 128-KiB-total MCU.
