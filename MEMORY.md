# STM32F407 memory layout

The default linker profile reserves SRAM1 for ordinary data and the main
stack, SRAM2 for decoded keys and signatures, and CCM for the transient
allocator. The signing workspace is zero-initialized in the SRAM1 `.bss`
and securely cleared after each operation. The guard serializes the public
operations. Host calls use automatic signing and ML2 storage.

| Bank | Start | Physical size | Contents |
| --- | --- | ---: | --- |
| SRAM1 | `0x20000000` | 112 KiB | `.data`, `.bss`, persistent signing state, main stack |
| SRAM2 | `0x2001c000` | 16 KiB | binding state in `.ram2` |
| CCM | `0x10000000` | 64 KiB | 61,440-byte `.ccm` arena with 61,416 bytes of initial payload capacity |
| Flash | `0x08000000` | 1 MiB | code and immutable precomputation |

The `.ram2` binding and `.ccm` arena use `NOLOAD`. The binding initializes
its own state and wipes it at release. The operation guard initializes the
arena and checks its state before returning. The linker profiles retain
their arena and binding assertions. The alternative full-RAM profile
places `end` after the binding reservation so stack spray excludes it.

## Target sizes

The following sizes are from the release `stm32f4discovery` test ELFs built
with GCC 13.2.1, `-mcpu=cortex-m4 -mthumb -mfloat-abi=hard
-mfpu=fpv4-sp-d16` and `-Os`. Struct sizes come from DWARF and sections from
`arm-none-eabi-size -A`. All sizes are in bytes.

| Object or section | I | III | V |
| --- | ---: | ---: | ---: |
| `ibz_t` | 272 | 424 | 552 |
| `fp_t` | 48 | 72 | 96 |
| ML2 fixed prefix | 6,816 | 9,096 | 11,016 |
| ML2 with 16 generators | 24,224 | 36,232 | 46,344 |
| Persistent signing state | 16,504 | 25,456 | 33,328 |
| Theta workspace | 17,280 | 25,920 | 34,560 |
| Binding (SRAM2) | 7,696 | 11,840 | 15,496 |
| `.data` + `.bss` | 17,228 | 26,172 | 34,052 |
| `[end, _stack)` reservation, before harness frames | 97,456 | 88,512 | 80,632 |
| Flash load (`.text`, `.ARM.exidx`, `.data`) | 137,180 | 156,076 | 172,604 |

The compiler uses 4-byte pointers and field digits, 8-byte `double`,
and 16-byte `dpe_t`. Target-compiled probes report 8-byte alignment for
`ibz_t`, `double`, `dpe_t`, the ML2 state, the signing state and the
binding. Pointers, field digits, `fp_t` and the theta workspace use
4-byte alignment. DPE retains a separate integer exponent. The hard
float ABI uses software binary64 operations in libgcc and the actual
Newlib `fegetround` and `fesetround` implementations. The target compiles
the portable 32-bit-half integer path.

## Allocation lifetimes

| Phase | Live CCM payload |
| --- | --- |
| Known-norm ideal construction | Eight incoming generators and ML2 with eight generator slots |
| Same-order intersection | Eight incoming generators and both reduced lattices during basis preparation. The second lattice is released before ML2 acquires eight generator slots |
| Ordered product reduction | ML2 with 16 slots. The shortened bases remain in the caller and each ordered product is generated directly into its destination slot |
| Prime-norm equivalent search | Two small-rank matrices. The final ideal multiplication can nest ML2 with four slots |
| Qlapoti shortest generator | Four IdealToIsogeny quaternion elements and two 4×4 reduction matrices |
| Qlapoti CVP search | Four IdealToIsogeny quaternion elements and five search candidates |
| IdealToIsogeny theta chain | Four IdealToIsogeny quaternion elements, eight Jacobian slots and nine theta slots |
| Verification and signing theta chains | Eight Jacobian slots and nine theta slots. Splits outside the offline table use a right subtree that fits the available slots |

The ML2 state uses `2736 + (15 + 4*d) * sizeof(ibz_t)` bytes for `d`
generators. The live payload for eight inputs plus its external generator
buffer is equal to the 16-generator ML2 payload. At Level V this is
46,344 bytes. The additional 9,384-byte reduced lattice brings the
intersection phase to 55,728 payload bytes, leaving 5,712 bytes in the
arena before allocator metadata and alignment. The IdealToIsogeny phase
has 11,040 bytes of quaternion elements alongside up to 34,560 bytes
for theta. The shortest-generator matrices, CVP candidates and theta
workspace are acquired and released in separate phases. The persistent
signing workspace remains live in SRAM1 across these phases. Each CCM
release clears the block, including on error exits from those phases.

The response computation can nest the intersection or ordered-product
frame, each of which can call basis reduction. KeyGen and commitment
conversion can nest Qlapoti and lattice membership under the IdealToIsogeny
frame. Lattice membership can in turn call the 4×4 adjugate routine. The
runtime canary includes the nested arithmetic, runtime helpers and
interrupts exercised by the measurement firmware.

| Level V compiler `.su` frame | Bytes |
| --- | ---: |
| IdealToIsogeny core | 13,520 |
| Qlapoti wrapper | 12,768 |
| Shortest generator | 1,136 |
| Qlapoti CVP search | 9,504 |
| Lattice membership | 13,824 |
| 4×4 adjugate | 16,656 |
| Response element computation | 23,760 |
| Ordered-product reduction | 18,240 |
| Same-order intersection | 12,784 |
| 4×4 basis reduction | 17,680 |

These are single-frame sizes. For example, the IdealToIsogeny core,
Qlapoti wrapper, CVP search, membership and adjugate frames account for
66,272 bytes before callers and integer helpers. The response computation,
intersection, basis reduction and its checked Gram frame account for
63,656 bytes before callers and integer helpers.

See [VALIDATION.md](VALIDATION.md) for the measured paths and test scope.
