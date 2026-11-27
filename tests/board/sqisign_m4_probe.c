/* SPDX-License-Identifier: Apache-2.0 */
#include "api.h"
#include "hal.h"
#include "m4_operation_guard.h"
#include "qlapoti_workspace.h"
#include "dpe.h"
#include "fp.h"

#include <fenv.h>
#include <float.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

/* The backend headers also name the internal NIST API. This harness calls
 * the pqm4 entry points declared by api.h. */
#undef crypto_sign_keypair
#undef crypto_sign
#undef crypto_sign_open

_Static_assert(FLT_RADIX == 2 && DBL_MANT_DIG == 53 && FLT_EVAL_METHOD == 0,
               "expected binary64 evaluation on Cortex-M4");

static const unsigned char message[] = "SQIsign-M4 board probe";
static unsigned char pk[CRYPTO_PUBLICKEYBYTES];
static unsigned char sk[CRYPTO_SECRETKEYBYTES];
static unsigned char signed_message[CRYPTO_BYTES + sizeof(message)];
static unsigned char opened[sizeof(message)];

extern char end, _stack;

static int
report_operation(const char *name, int status, size_t stack)
{
    const qlw_stats_t arena = sqisign_m4_last_workspace_stats;
    char line[160];
    snprintf(line, sizeof(line), "%s status=%d stack=%u capacity=%u current=%u peak=%u live=%u arena=%u",
             name, status, (unsigned)stack, (unsigned)arena.capacity_bytes,
             (unsigned)arena.current_bytes, (unsigned)arena.peak_bytes,
             (unsigned)arena.live_blocks, arena.status_flags);
    hal_send_str(line);
    return arena.current_bytes == 0 && arena.live_blocks == 0 && arena.status_flags == 0;
}

int
main(void)
{
    char line[120];
    size_t signed_length = 0, opened_length = 0;
    size_t stack;
    int status;
    int passed = 0;

    hal_setup(CLOCK_FAST);
    hal_send_str("==========================");
    snprintf(line, sizeof(line), "sizes ibz=%u fp=%u digit=%u pointer=%u double=%u dpe=%u stack_bank=%u",
             (unsigned)sizeof(ibz_t), (unsigned)sizeof(fp_t), (unsigned)sizeof(digit_t),
             (unsigned)sizeof(void *), (unsigned)sizeof(double), (unsigned)sizeof(dpe_t),
             (unsigned)((uintptr_t)&_stack - (uintptr_t)&end));
    hal_send_str(line);

    /* Use volatile operands so this exercises the linked binary64 runtime. */
    volatile double one = 1.0, halfway = 0x1p-53;
    if (fesetround(FE_TONEAREST) != 0 || fegetround() != FE_TONEAREST || one + halfway != 1.0) {
        hal_send_str("rounding FAIL");
        goto done;
    }
    hal_send_str("rounding PASS");

    hal_spraystack();
    status = crypto_sign_keypair(pk, sk);
    stack = hal_checkstack();
    if (!report_operation("keypair", status, stack) || status != 0)
        goto done;

    memcpy(signed_message, message, sizeof(message));
    hal_spraystack();
    status = crypto_sign(signed_message, &signed_length, signed_message, sizeof(message), sk);
    stack = hal_checkstack();
    if (!report_operation("sign", status, stack) || status != 0 || signed_length != sizeof(signed_message))
        goto done;

    hal_spraystack();
    status = crypto_sign_open(opened, &opened_length, signed_message, signed_length, pk);
    stack = hal_checkstack();
    if (!report_operation("verify", status, stack) || status != 0 ||
        opened_length != sizeof(message) || memcmp(opened, message, sizeof(message)) != 0)
        goto done;

    signed_message[signed_length - 1] ^= 1;
    hal_spraystack();
    status = crypto_sign_open(opened, &opened_length, signed_message, signed_length, pk);
    stack = hal_checkstack();
    passed = report_operation("verify-tampered", status, stack) && status != 0 && opened_length == 0;

done:
    for (size_t i = 0; i < sizeof(sk); ++i)
        ((volatile unsigned char *)sk)[i] = 0;
    hal_send_str(passed ? "RESULT PASS" : "RESULT FAIL");
    hal_send_str("#");
    return passed ? 0 : 1;
}
