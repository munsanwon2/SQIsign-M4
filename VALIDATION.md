# Validation

The Cortex-M4 import starts from `a83101a14380e9c05ac878a31b7effa3c89056f9`.
The compact source and pqm4 revisions are in [SOURCE_URLS.lock](SOURCE_URLS.lock).
The host builds used x86-64 Ubuntu 24.04 under WSL, GCC 13.3.0 and CMake 3.28.3.
The STM32F407 firmware used GNU Arm Embedded GCC 13.2.1 and Newlib.

## Reproduction and host tests

```sh
python3 full-backend-tools/rebuild_compact_exports.py \
    --upstream /path/to/compact-SQIsign --output /tmp/sqisign-exports
diff -ru backends /tmp/sqisign-exports
bash scripts/verify_release.sh
bash scripts/install_all_into_pqm4.sh /path/to/pinned/pqm4
OPT_SIZE=1 bash scripts/run_pqm4_validation.sh /path/to/pinned/pqm4 stm32f4discovery all
bash tests/quaternion/run.sh
python3 tests/encoding/run.py /path/to/pinned/pqm4
```

- The export and installer verified the source closures, included fragments,
  hashes, active ARM dependencies and one provider for each backend API.
- The clean pqm4 build linked the `test`, `speed`, `hashing`, `stack` and
  `testvectors` ELF for each level, and the three host testvector binaries.
- All 12 deterministic pqm4 host cases passed for each level. The pinned
  harness uses its SURF seed and tests message lengths 0, 1, 2, 4 through
  1024 with full KeyGen, Sign, Verify and message recovery.
- The arithmetic, encoding and workspace suites passed. The tests cover the
  original normeq and dependent-row regressions, exact ML2 reducedness
  and Gram updates, the fallback limb path, DPE rounding, direct sampler
  support and failure statuses, Fp and Fp2 against Python integers, scalar
  ideal products and intersections against exact expected lattices, field
  canonicality, signed SK widths, widened exponents, scalar inversion,
  the bounded theta strategy and workspace recovery from busy, OOM and
  live-allocation failures.

Both the AIO and backend-only GCC commands carry `-fno-fast-math`,
`-ffp-contract=off`, `-frounding-math` and `-fexcess-precision=standard`.
The release profile uses `NDEBUG`. For an assertion and overflow-checking
build, start with a clean pqm4 build directory and set `SQISIGN_CHECKED=1`.
Its backend and AIO commands end the `NDEBUG` settings with `-UNDEBUG`
and define `SQISIGN_INTBIG_OVERFLOW_CHECK`. All three checked ARM test
ELFs linked with an assertion trap. All 36 checked host vectors matched
their release output. The board measurements use the release build.

With a trace on the host IdealToIsogeny return path, the same vectors
exercised both Qlapoti outcomes and retained byte-for-byte host parity.

| Level | Normal calls / theta exponent | Transformed calls / theta exponent |
| --- | ---: | ---: |
| I | 28 / 322 | 8 / 321 |
| III | 25 / 498 | 11 / 497 |
| V | 31 / 662 | 5 / 661 |

| Level | Library manifest SHA-256 | STM32F407 `test.elf` SHA-256 in the validation checkout |
| --- | --- | --- |
| I | `92c577823eb22ad0657a6c63e4e4cae001b38ff4350910f9015690a6e697e3ee` | `499d04c45faf7a22cc594671b7dbe8231b4f3256edd9c2d9183c0409a3875b46` |
| III | `52ba8e7cce89a14fb142f8e5f9d2eed0d0eeb1489fb21903fbee09cbc4908470` | `1dc055bb4d5755f807649f8edfb63d9a6f8085056c5b88b7627c4ef866a8ddbf` |
| V | `373de8ec587ba83a36cb403907c83a876fbc1144160858a87992b0a172cf7702` | `d0630b4fa7b4cff00b52902d98fa26313e1f39c435735249bb5f03782db8774f` |

## Board measurements

The connected STM32F4 ran the `stm32f4discovery` STM32F407G-DISC1
profile with the banked linker script and pqm4's 168 MHz `CLOCK_FAST`
setup. The measurements use pqm4's stack spray
and the workspace allocator's per-operation statistics. Test firmware
has its own data and stack footprint. The SRAM1 sizes and statically
allocated objects in [MEMORY.md](MEMORY.md) are from the standard test ELFs.

The board probe in `tests/board/sqisign_m4_probe.c` exercises hardware
entropy, an overlapping signed-message call, successful verification and
message tampering. It checks the linked binary64 rounding environment
and arena cleanup after each operation. The UART profile uses 38400 baud.
SWD can collect the same `hal_send_str` calls when the serial adapter is
unavailable.

To build the hardware-entropy probe in a validation checkout of pqm4, copy
it to `mupq/crypto_sign/sqisign_m4_probe.c`. For deterministic parity, the
generator instruments the pinned pqm4 harness and retains its first vector
and SURF seed.

```sh
python3 tests/board/make_vector_probe.py /path/to/pinned/pqm4
cd /path/to/pinned/pqm4
scheme=sqisign-qlapoti-lvl5
make PLATFORM=stm32f4discovery OPT_SIZE=1 \
    IMPLEMENTATION_PATH="crypto_sign/$scheme/m4" \
    "elf/crypto_sign_${scheme}_m4_measure_testvectors.elf"
```

The hardware-entropy probe uses the target suffix `_m4_sqisign_m4_probe.elf`.
Flash the resulting binary and collect the output from the UART or SWD.
For SWD, run `st-util -p 42425` and invoke
`SQI_TRANSCRIPT=/tmp/sqisign-probe.txt gdb-multiarch -batch <elf> -x <portkit>/tests/board/capture_swd.py`.
Compare the three hexadecimal records of the deterministic probe with the
first three lines of the corresponding host testvector binary. The reports
below store hashes and resource measurements.

### First deterministic vector on the board

The first pqm4 vector at each level completed on the connected board.
The public key, secret key and signed-message bytes all matched the host
vector. SWD captured the output. The hashes use the three hexadecimal lines
joined by linefeeds, with no trailing linefeed.

| Level | First vector SHA3-256 | Flashed probe binary SHA-256 |
| --- | --- | --- |
| I | `2c73cf6aa856beb3da3f3f3b1de821bed8352975de07ca5437b4edf3ef6a76f7` | `886ca2b9e6d73080a0017de7ade4a49ecc28e808c74e3dba23f8fe2dc5fd321f` |
| III | `1baf5b6fbd40bf28eb04193558e6e568a3be967225cd8db6cb4b25cf3a0abb09` | `c3ceaa40c679ff82babec9d709c2b1fff5a46cbe62ad20dd5c9621c5b1143d8d` |
| V | `2a59e672cee68cc2f9d6468fea3875539cf11bd2f0954501a2be75c2b8983acd` | `0197defbe22425cf041eaa950cf1dbdcaa56349784e9804de4f2cd9c37c09f02` |

The stack available column is measured immediately before the stack spray
and accounts for the live harness frame and any HAL reservation. Headroom
is the difference between that value and the observed high-water mark.
CCM peaks count payload bytes. All nine operations reported a capacity of
61,416 bytes with zero current bytes, live allocations and status flags on
return. These measurements cover the paths and interrupts exercised by the
probe at `CLOCK_FAST`.

| Level | Operation | Stack available | Stack peak | Stack headroom | CCM peak |
| --- | --- | ---: | ---: | ---: | ---: |
| I | keypair | 95,028 | 35,160 | 59,868 | 24,224 |
| I | sign | 95,028 | 35,008 | 60,020 | 28,848 |
| I | verify | 95,028 | 21,200 | 73,828 | 17,280 |
| III | keypair | 85,764 | 54,064 | 31,700 | 36,232 |
| III | sign | 85,764 | 54,312 | 31,452 | 43,440 |
| III | verify | 85,764 | 31,500 | 54,264 | 25,920 |
| V | keypair | 77,492 | 70,456 | 7,036 | 46,344 |
| V | sign | 77,492 | 70,568 | 6,924 | 55,728 |
| V | verify | 77,492 | 41,800 | 35,692 | 34,560 |

The Level I hardware-entropy probe also passed on the same board. It
exercised a nonempty message, overlapping input and output, recovery of the
message and rejection after a message byte was changed. The Newlib binary64
round-to-nearest and halfway addition check passed. The guard reported a
clean arena on all four calls.

## ARM emulation

The Cortex-M4 `mps2-an386` QEMU run used a custom linker layout with
112 KiB for ordinary data and the main stack, 16 KiB for the binding,
64 KiB for the CCM arena and 1 MiB for code. The arena used an emulated
RAM region at `0x20100000`. QEMU uses its own runtime and interrupt model.

All 12 vectors per level matched the corresponding portable host bytes.
Level V ran nine cases sequentially and the remaining three from the
corresponding deterministic SURF states. Each case executed KeyGen, Sign
and Verify. All operations completed with no live arena blocks or error
flags.

| Level | Vectors | Largest observed stack use | Largest observed CCM payload |
| --- | ---: | ---: | ---: |
| I | 12 | 35,160 | 28,848 |
| III | 12 | 54,312 | 43,440 |
| V | 12 | 70,568 | 55,728 |
