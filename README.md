# SQIsign-Qlapoti for Cortex-M4

This repository contains a fixed-precision SQIsign port for the Arm Cortex-M4,
packaged as three `pqm4` signature schemes:

- `sqisign-qlapoti-lvl1` (NIST Level I)
- `sqisign-qlapoti-lvl3` (NIST Level III)
- `sqisign-qlapoti-lvl5` (NIST Level V)

The implementation combines the compact-SQIsign algorithms with Qlapoti and a
bounded-memory Cortex-M4 integration. It avoids GMP and hosted heap services in
the target backend. The STM32F407 linker profile uses the device's split SRAM
banks so that the complete 192 KiB physical SRAM is available to the program.

This is research software, not a production cryptographic library. Validate
functional correctness, timing, stack headroom, and side-channel requirements
for the exact toolchain and board configuration before deployment. The current
target is the STM32F407G-DISC1. Operation on a 128 KiB-only device is not claimed.

## Repository layout

- `backends/lvl{1,3,5}`: content-addressed, level-specific backend source exports
- `pqm4-overlay`: pqm4 API wrappers, workspace allocator, fixed-precision helpers,
  configuration, and STM32F407 linker profiles
- `scripts`: installation, removal, validation, and release-audit helpers
- `tools` and `full-backend-tools`: source/manifest validation and wrapper generation

## Prerequisites

The commands below assume a Linux host or WSL2 with:

- Git, GNU Make, Bash, and Python 3
- `arm-none-eabi-gcc`, binutils, and Newlib
- ST-LINK tools supported by pqm4
- an external 3.3 V USB-to-UART adapter for pqm4 serial output

Use a pqm4 checkout whose absolute path contains no whitespace. This release is
validated against pqm4 commit `cc2c1b992b602d285bd15991a566d5f17b34c1fa`.

```sh
git clone --recursive https://github.com/mupq/pqm4.git /tmp/pqm4-sqisign
git -C /tmp/pqm4-sqisign checkout cc2c1b992b602d285bd15991a566d5f17b34c1fa
git -C /tmp/pqm4-sqisign submodule update --init --recursive

git clone https://github.com/munsanwon2/SQIsign-M4.git
cd SQIsign-M4
./scripts/install_all_into_pqm4.sh /tmp/pqm4-sqisign
```

The installer adds only the three scheme directories and two SQIsign-specific
linker scripts. It refuses to overwrite an existing scheme or a conflicting
linker profile. To install one level, pass `1`, `3`, or `5` as the second
argument. To remove marked files installed by this repository:

```sh
./scripts/uninstall_from_pqm4.sh /tmp/pqm4-sqisign
```

## Build and host-side validation

Run the release checks before installation when preparing a new checkout:

```sh
./scripts/verify_release.sh
```

The following command audits the installed source, builds the target test,
speed, hashing, stack, and test-vector ELFs for all levels, and builds the
portable host test-vector executables:

```sh
./scripts/run_pqm4_validation.sh /tmp/pqm4-sqisign stm32f4discovery all
```

To build one target manually:

```sh
cd /tmp/pqm4-sqisign
make PLATFORM=stm32f4discovery OPT_SIZE=1 \
  IMPLEMENTATION_PATH=crypto_sign/sqisign-qlapoti-lvl1/m4 \
  elf/crypto_sign_sqisign-qlapoti-lvl1_m4_test.elf \
  elf/crypto_sign_sqisign-qlapoti-lvl1_m4_speed.elf
```

Replace `lvl1` with `lvl3` or `lvl5` for the other parameter sets. The supplied
configuration selects size optimization and the SQIsign-specific STM32F407
linker layout.

## Connect the STM32F407G-DISC1

Connect the board's integrated ST-LINK USB port for programming. Stock pqm4
uses a separate UART at 38,400 baud for test output. Connect a 3.3 V USB-UART
adapter as follows:

| STM32F407G-DISC1 | USB-UART adapter |
| --- | --- |
| PA2 / USART2 TX | RX |
| PA3 / USART2 RX | TX (optional for output-only runs) |
| GND | GND |

Do not connect a 5 V UART signal to the MCU. On Linux, the adapter commonly
appears as `/dev/ttyUSB0`; grant the current user serial-device access as
required by the distribution.

For WSL2, expose both ST-LINK and the UART adapter from an elevated Windows
PowerShell prompt, then attach them from a normal prompt:

```powershell
usbipd list
usbipd bind --busid <BUSID>
usbipd attach --wsl --busid <BUSID>
```

Confirm inside WSL that the devices are visible before running pqm4.

## Functional tests on the board

From the pqm4 checkout, run one scheme at a time:

```sh
cd /tmp/pqm4-sqisign
python3 test.py -p stm32f4discovery -o size -u /dev/ttyUSB0 \
  sqisign-qlapoti-lvl1
python3 testvectors.py -p stm32f4discovery -o size -u /dev/ttyUSB0 \
  sqisign-qlapoti-lvl1
```

Repeat with `sqisign-qlapoti-lvl3` and `sqisign-qlapoti-lvl5`.

## Benchmarking

The following uses pqm4's standard speed harness for 100 iterations while
disabling the separate stack, hashing, and code-size benchmark passes:

```sh
cd /tmp/pqm4-sqisign
python3 benchmarks.py -p stm32f4discovery -o size -i 100 \
  -u /dev/ttyUSB0 --nostack --nohashing --nosize \
  sqisign-qlapoti-lvl1
```

Repeat for Levels III and V. pqm4 writes raw results beneath
`benchmarks/speed/crypto_sign/<scheme>/<implementation>/`; preserve those files
with the toolchain version, pqm4 commit, repository commit, board revision,
clock configuration, compiler flags, message length, and iteration count.

For comparable measurements, use the same board, power source, firmware,
toolchain, optimization mode, UART setup, and pqm4 commit for every run. At the
pinned revision, the stock `stm32f4discovery` speed harness runs the Cortex-M4F
at 24 MHz with zero flash wait states, measures SysTick cycle deltas, obtains
fresh randomness from the STM32 hardware RNG, and uses a 59-byte message. Each
iteration performs KeyGen, then Sign, then Verify; the command above selects
the port's `-Os` configuration. Record all 100 raw cycle counts and compute
summary statistics from those values. Do not compare them directly with results
from a custom clock profile or a host wall-clock benchmark.

To collect stack measurements as well, omit `--nostack`. Stack results remain a
required validation gate for a memory-constrained deployment even when the
reported performance table contains speed figures only.

## Memory and integration notes

The wrappers route transient allocations through a bounded workspace and do not
depend on a conventional target heap. The custom STM32F407 scripts account for
the device's physically split SRAM. Link success alone is not a proof of runtime
headroom: measure the high-water mark for KeyGen, Sign, and Verify on the final
firmware and reserve space for active interrupt frames and application state.

## Licensing and upstreams

Repository glue and tools are provided under the Apache License 2.0; see
`LICENSE`. Imported source retains its upstream terms and notices under
`third_party_licenses`. The implementation was developed from:

- [compact-SQIsign](https://github.com/pqc-lab-ku/compact-SQIsign)
- [Qlapoti](https://github.com/KULeuven-COSIC/Qlapoti)
- [SQIsign vectorization](https://github.com/LeeJ-art/the-sqisign-vectorization)
- [pqm4](https://github.com/mupq/pqm4)
