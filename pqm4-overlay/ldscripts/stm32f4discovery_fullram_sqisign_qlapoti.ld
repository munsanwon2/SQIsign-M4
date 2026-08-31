EXTERN(vector_table)
ENTRY(reset_handler)

MEMORY
{
    ram (rwx) : ORIGIN = 0x20000000, LENGTH = 128K
    rom (rx) : ORIGIN = 0x08000000, LENGTH = 1024K
    ccm (rwx) : ORIGIN = 0x10000000, LENGTH = 64K
}

SECTIONS
{
    .text : {
        *(.vectors)
        *(.text*)
        . = ALIGN(4);
        *(.rodata*)
        . = ALIGN(4);
    } >rom

    .preinit_array : {
        . = ALIGN(4);
        __preinit_array_start = .;
        KEEP (*(.preinit_array))
        __preinit_array_end = .;
    } >rom

    .init_array : {
        . = ALIGN(4);
        __init_array_start = .;
        KEEP (*(SORT(.init_array.*)))
        KEEP (*(.init_array))
        __init_array_end = .;
    } >rom

    .fini_array : {
        . = ALIGN(4);
        __fini_array_start = .;
        KEEP (*(.fini_array))
        KEEP (*(SORT(.fini_array.*)))
        __fini_array_end = .;
    } >rom

    .ARM.extab : { *(.ARM.extab*) } >rom
    .ARM.exidx : {
        __exidx_start = .;
        *(.ARM.exidx*)
        __exidx_end = .;
    } >rom

    . = ALIGN(4);
    _etext = .;
    .data : {
        _data = .;
        *(.data*)
        . = ALIGN(4);
        _edata = .;
    } >ram AT >rom
    _data_loadaddr = LOADADDR(.data);

    .bss : {
        *(.bss*)
        *(COMMON)
        . = ALIGN(4);
        _ebss = .;
    } >ram

    /* The 128 KiB fallback treats SRAM1+SRAM2 as one contiguous region, but
       preserves the dedicated section so its cost remains visible in maps. */
    .ram2 (NOLOAD) : {
        . = ALIGN(8);
        __sqisign_binding_start__ = .;
        KEEP(*(.ram2.sqisign_binding))
        __sqisign_binding_end__ = .;
        *(.ram2*)
        . = ALIGN(4);
        /* pqm4's stack spray and _sbrk start at `end`; it must account for
           this reservation in the contiguous 128-KiB fallback profile. */
        end = .;
    } >ram

    /* qlw_reset()/initialize_if_needed() zero the arena in software.  Keep
       this uninitialised storage out of the flash image. */
    .ccm (NOLOAD) : {
        . = ALIGN(16);
        __qlw_workspace_start__ = .;
        KEEP(*(.ccmram.qlapoti_workspace))
        __qlw_workspace_end__ = .;
        *(.ccmram*)
        . = ALIGN(4);
    } >ccm

    /DISCARD/ : { *(.eh_frame) }
    . = ALIGN(4);
}

PROVIDE(_stack = ORIGIN(ram) + LENGTH(ram));
ASSERT((__qlw_workspace_end__ - __qlw_workspace_start__) == 61440,
       "Qlapoti workspace missing or not exactly 61440 bytes");
ASSERT((__sqisign_binding_end__ - __sqisign_binding_start__) <= 16K,
       "SQIsign binding state exceeds 16 KiB");
