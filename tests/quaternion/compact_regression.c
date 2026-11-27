/* SPDX-License-Identifier: Apache-2.0 */
#include "dpe.h"
#include "lll_internals.h"
#include "m4_operation_guard.h"
#include "quaternion.h"

#include <fenv.h>
#include <stdio.h>
#include <stdlib.h>

static void
require(int test, const char *name)
{
    if (!test) {
        fprintf(stderr, "compact regression: %s\n", name);
        exit(EXIT_FAILURE);
    }
}

static void
check_integer_boundaries(void)
{
    ibz_t max, min = {0}, result, expected;
    min[IBZ_LIMBS - 1] = UINT64_C(1) << 63;
    for (size_t i = 0; i < IBZ_LIMBS; ++i)
        max[i] = ~min[i];
    ibz_set(&result, 91);
    require(!ibz_add_checked(&result, &max, &ibz_const_one) && ibz_get(&result) == 91,
            "maximum plus one publishes an overflow");
    require(!ibz_sub_checked(&result, &min, &ibz_const_one) && ibz_get(&result) == 91,
            "minimum minus one publishes an overflow");
    require(!ibz_mul_checked(&result, &max, &ibz_const_two) && ibz_get(&result) == 91,
            "multiply overflow in a release-style checked call");
    require(ibz_mul_checked(&result, &min, &ibz_const_one) && ibz_cmp(&result, &min) == 0,
            "minimum times one fits exactly");

    uint32_t words[2 * IBZ_LIMBS];
    uint64_t limbs[IBZ_LIMBS];
    require(!ibz_to_u32_digits_checked(words, 2 * IBZ_LIMBS - 1, &max),
            "32-bit destination capacity enforced");
    require(ibz_to_u32_digits_checked(words, 2 * IBZ_LIMBS, &max) &&
            ibz_to_u64_digits_checked(limbs, IBZ_LIMBS, &max), "both digit widths export capacity");
    require(limbs[0] == UINT64_MAX && words[2 * IBZ_LIMBS - 1] == UINT32_MAX / 2,
            "64-bit limb and 32-bit digit representation");
    ibz_copy_u32_digits(&result, words, 2 * IBZ_LIMBS);
    ibz_copy_u64_digits(&expected, limbs, IBZ_LIMBS);
    require(ibz_cmp(&result, &max) == 0 && ibz_cmp(&expected, &max) == 0, "digit round trip");
}

static void
check_dpe(void)
{
    ibz_t input, rounded, expected;
    dpe_t approximate;
    const uint64_t half = UINT64_C(1) << 53;
    for (int negative = 0; negative <= 1; ++negative) {
        ibz_set_u64(&input, half + 1);
        ibz_set_u64(&expected, half);
        if (negative) {
            ibz_neg(&input, &input);
            ibz_neg(&expected, &expected);
        }
        dpe_set_z(approximate, &input);
        require(dpe_get_z(&rounded, approximate) && ibz_cmp(&rounded, &expected) == 0,
                "positive and negative halfway round to the even significand");
    }
    ibz_set_u64(&input, 2 * half - 1);
    ibz_set_u64(&expected, 2 * half);
    dpe_set_z(approximate, &input);
    require(dpe_get_z(&rounded, approximate) && ibz_cmp(&rounded, &expected) == 0,
            "significand carry adjusts the exponent");
    ibz_mul_2exp(&input, &ibz_const_one, IBZ_BITS - 80);
    dpe_set_z(approximate, &input);
    require(dpe_get_z(&rounded, approximate) && ibz_cmp(&rounded, &input) == 0,
            "large exponents survive a DPE conversion");
    dpe_set_z(approximate, &ibz_const_zero);
    require(dpe_get_z(&rounded, approximate) && ibz_is_zero(&rounded), "zero sentinel round trip");
    dpe_set_d(approximate, 0.5);
    DPE_EXP(approximate) = IBZ_BITS;
    require(!dpe_get_z(&rounded, approximate) && ibz_is_zero(&rounded), "DPE capacity rejection");
}

static void
check_ml2_admission(void)
{
    ibz_vec_4_t input = {{0}}, output;
    const int exponent = IBZ_BITS / 2 - 6;
    ibz_mul_2exp(&input[0], &ibz_const_one, exponent);
    const int rounding = fegetround();
    require(rounding >= 0 && fesetround(FE_UPWARD) == 0, "caller fenv setup");
    require(quat_ml2(&output, 1, &input, 1, NULL) == 1 && fegetround() == FE_UPWARD,
            "512-bound admits the boundary and restores the caller rounding mode");
    require(fesetround(rounding) == 0, "caller fenv restore");
    ibz_mul_2exp(&input[0], &ibz_const_one, exponent + 1);
    ibz_set(&output[0], 91);
    require(quat_ml2(&output, 1, &input, 1, NULL) == QUAT_ML2_ERR_PRECISION && ibz_get(&output[0]) == 91,
            "512-bound failure preserves the output");
}

static unsigned char rng_choices[2];
static unsigned rng_call;
static int rng_fail;

int
randombytes(unsigned char *output, unsigned long long length)
{
    if (rng_fail)
        return -1;
    memset(output, rng_choices[rng_call++ % 2], (size_t)length);
    return 0;
}

static void
check_direct_sampler(void)
{
    quat_alg_t alg;
    quat_alg_init_set_ui(&alg, 17);
    quat_lattice_t lattice = {0};
    quat_alg_elem_t sample;
    ibz_set(&lattice.denom, 1);
    for (int i = 0; i < 4; ++i)
        ibz_set(&lattice.basis[i][i], 1);

    /* For B=1 and the weighted identity, the coefficient boxes are
     * [-3,3], [-2,2], {0}, {0}. The norm-one boundary has four points. */
    for (int x = -1; x <= 1; ++x) {
        for (int y = -1; y <= 1; ++y) {
            if (x * x + y * y > 1)
                continue;
            rng_choices[0] = (unsigned char)(3 + x);
            rng_choices[1] = (unsigned char)(2 + y);
            rng_call = 0;
            require(quat_lattice_sample_from_ball_direct(&sample, &lattice, &alg, &ibz_const_one) ==
                        QUAT_BALL_SUCCESS && ibz_cmp_int32(&sample.coord[0], x) == 0 &&
                        ibz_cmp_int32(&sample.coord[1], y) == 0 && ibz_is_zero(&sample.coord[2]) &&
                        ibz_is_zero(&sample.coord[3]) && rng_call == 2,
                    "direct box covers zero and every boundary point");
        }
    }
    rng_fail = 1;
    require(quat_lattice_sample_from_ball_direct(&sample, &lattice, &alg, &ibz_const_one) ==
                QUAT_BALL_RNG_FAILURE, "direct RNG failure is distinct");
    rng_fail = 0;
    rng_choices[0] = rng_choices[1] = 0;
    require(quat_lattice_sample_from_ball_direct(&sample, &lattice, &alg, &ibz_const_one) ==
                QUAT_BALL_EXHAUSTED, "coordinate rejection has a finite proposal limit");
    require(quat_lattice_sample_from_ball_direct(&sample, &lattice, &alg, &ibz_const_zero) ==
                QUAT_BALL_INVALID, "zero norm bound is invalid");
}

static int
unexpected_core(void *context)
{
    (void)context;
    require(0, "entry ran while another operation was live");
    return 0;
}

static int
check_nested_guard(void *context)
{
    (void)context;
    unsigned char *block = qlw_malloc(128);
    require(block != NULL, "workspace allocation succeeds");
    for (size_t i = 0; i < 128; ++i)
        require(block[i] == 0, "operation clears previous workspace data");
    block[127] = 0x5a;
    require(sqisign_m4_guard_run(unexpected_core, NULL) == SQISIGN_M4_GUARD_BUSY &&
                qlw_get_stats().live_blocks == 1 && block[127] == 0x5a,
            "busy return preserves the active allocation");
    qlw_free(block);
    return 0;
}

static int
force_workspace_oom(void *context)
{
    (void)context;
    require(qlw_malloc(QLW_WORKSPACE_BYTES) == NULL, "oversized allocation is refused");
    return 0;
}

static int
leave_workspace_allocation(void *context)
{
    (void)context;
    unsigned char *block = qlw_malloc(128);
    require(block != NULL, "allocation for the cleanup test");
    block[127] = 0xa5;
    return 0;
}

static void
check_operation_guard(void)
{
    require(sqisign_m4_guard_run(check_nested_guard, NULL) == SQISIGN_M4_GUARD_OK &&
                sqisign_m4_last_workspace_stats.live_blocks == 0,
            "successful outer operation releases its workspace");
    require(sqisign_m4_guard_run(force_workspace_oom, NULL) == SQISIGN_M4_GUARD_WORKSPACE_FAILED &&
                sqisign_m4_last_workspace_stats.status_flags == QLW_STATUS_OOM,
            "allocation exhaustion reaches the public guard");
    require(sqisign_m4_guard_run(leave_workspace_allocation, NULL) == SQISIGN_M4_GUARD_WORKSPACE_FAILED &&
                sqisign_m4_last_workspace_stats.status_flags == QLW_STATUS_LEAK,
            "live allocation makes the operation fail");
    require(sqisign_m4_guard_run(check_nested_guard, NULL) == SQISIGN_M4_GUARD_OK &&
                sqisign_m4_last_workspace_stats.status_flags == QLW_STATUS_OK,
            "operation guard recovers and clears the arena");
}

int
main(void)
{
    check_integer_boundaries();
    check_dpe();
    check_ml2_admission();
    check_direct_sampler();
    check_operation_guard();
    puts("compact arithmetic regression: PASS");
    return EXIT_SUCCESS;
}
