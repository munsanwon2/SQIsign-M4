/* SPDX-License-Identifier: Apache-2.0 */
#include "internal.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* Link-time wrappers control only candidate selection.  The production
 * normeq.c, fixed-width arithmetic, parity check, primitive conversion and
 * output publication all run unchanged.  No test RNG enters a target build. */
static ibz_vec_4_t candidate;
static ibz_t cornacchia_target;
static unsigned samples, primality_calls, cornacchia_calls;

static void
require(int ok, const char *message)
{
    if (!ok) {
        fprintf(stderr, "normeq regression: %s\n", message);
        exit(EXIT_FAILURE);
    }
}

int
__wrap_ibz_rand_interval(ibz_t *out, const ibz_t *lower, const ibz_t *upper)
{
    /* Reject a second attempt: a wrongly rejected candidate must not be
     * hidden by a later random success or by a million-iteration search. */
    if (samples >= 2)
        return 0;
    const ibz_t *value = &candidate[2 + samples++];
    require(ibz_cmp(value, lower) >= 0 && ibz_cmp(value, upper) <= 0,
            "fixture outside the production sampling interval");
    ibz_copy(out, value);
    return 1;
}

int
__wrap_ibz_probab_prime(const ibz_t *n, int reps)
{
    ++primality_calls;
    require(reps > 0 && ibz_cmp(n, &cornacchia_target) == 0,
            "unexpected primality-test input");
    /* This isolates the congruence gate; primality is not under test. */
    return 1;
}

int
__wrap_ibz_cornacchia_prime(ibz_t *x, ibz_t *y,
                           const ibz_t *q, const ibz_t *n)
{
    ++cornacchia_calls;
    require(ibz_is_one(q) && ibz_cmp(n, &cornacchia_target) == 0,
            "unexpected Cornacchia input");
    /* The fixture satisfies x^2 + y^2 = n exactly, even when n is composite. */
    ibz_copy(x, &candidate[0]);
    ibz_copy(y, &candidate[1]);
    return 1;
}

/* An unexpected call into the actual entropy path is a test failure. */
int
randombytes(unsigned char *out, unsigned long long length)
{
    (void)out;
    (void)length;
    require(0, "unexpected entropy request");
    return -1;
}

static void
check_candidate(const char *name, const uint64_t values[4],
                int wide, int accepted)
{
    quat_alg_t alg;
    quat_p_extremal_maximal_order_t order = {0};
    quat_alg_elem_t gamma, sentinel;
    ibz_t norm, term, remainder, four, numerator, denominator;
    quat_alg_init_set_ui(&alg, 7);
    quat_alg_elem_init(&gamma);
    quat_alg_elem_init(&sentinel);
    ibz_init(&norm);
    ibz_init(&term);
    ibz_init(&remainder);
    ibz_init(&four);
    ibz_init(&numerator);
    ibz_init(&denominator);
    ibz_set(&four, 4);
    quat_lattice_O0_set_extremal(&order);
    const quat_represent_integer_params_t params = {
        .algebra = &alg, .order = &order, .primality_test_iterations = 1
    };

    for (int i = 0; i < 4; ++i)
        ibz_set_u64(&candidate[i], values[i]);
    if (wide) {
        ibz_set(&term, 1);
        ibz_mul_2exp(&term, &term, 96);
        ibz_add(&candidate[0], &candidate[0], &term);
        ibz_add(&candidate[3], &candidate[3], &term);
    }
    ibz_set(&norm, 0);
    ibz_set(&cornacchia_target, 0);
    for (int i = 0; i < 4; ++i) {
        ibz_mul(&term, &candidate[i], &candidate[i]);
        if (i < 2)
            ibz_add(&cornacchia_target, &cornacchia_target, &term);
        else
            ibz_mul(&term, &term, &alg.p);
        ibz_add(&norm, &norm, &term);
    }
    ibz_div(&norm, &remainder, &norm, &four);
    require(ibz_is_zero(&remainder) && ibz_is_odd(&norm),
            "fixture must have positive odd integral norm");
    quat_alg_elem_set(&sentinel, 1, 101, 103, 107, 109);
    quat_alg_elem_copy(&gamma, &sentinel);
    samples = primality_calls = cornacchia_calls = 0;
    int result = quat_represent_integer(&gamma, &norm, 1, &params);
    if (result != accepted) {
        fprintf(stderr, "normeq fixture failed: %s (got %d, expected %d)\n",
                name, result, accepted);
        exit(EXIT_FAILURE);
    }
    require(samples == 2 && primality_calls == 1 && cornacchia_calls == 1,
            "fixture did not reach the production congruence gate exactly once");
    if (accepted) {
        require(ibz_cmp_int32(&gamma.denom, 2) == 0,
                "unexpected output denominator");
        int swapped = ibz_is_odd(&candidate[0]) != ibz_is_odd(&candidate[3]);
        for (int i = 0; i < 4; ++i) {
            int source = i < 2 && swapped ? 1 - i : i;
            ibz_copy(&term, &candidate[source]);
            if (i == 3)
                ibz_neg(&term, &term);
            require(ibz_cmp(&gamma.coord[i], &term) == 0,
                    "accepted candidate was not published correctly");
        }
        require(quat_alg_norm(&numerator, &denominator, &gamma, &alg),
                "output norm calculation failed");
        ibz_mul(&term, &norm, &denominator);
        require(ibz_cmp(&numerator, &term) == 0, "output norm changed");
    } else {
        require(ibz_cmp(&gamma.denom, &sentinel.denom) == 0,
                "rejected candidate modified the denominator");
        for (int i = 0; i < 4; ++i)
            require(ibz_cmp(&gamma.coord[i], &sentinel.coord[i]) == 0,
                    "rejected candidate modified an output coordinate");
    }
    ibz_finalize(&denominator);
    ibz_finalize(&numerator);
    ibz_finalize(&four);
    ibz_finalize(&remainder);
    ibz_finalize(&term);
    ibz_finalize(&norm);
    quat_alg_elem_finalize(&sentinel);
    quat_alg_elem_finalize(&gamma);
    quat_alg_finalize(&alg);
}

static void
check_negative_residues(void)
{
    static const int32_t values[] = {0, -1, -2, -3, -4, INT32_MIN, INT32_MAX};
    static const unsigned long residues[] = {0, 3, 2, 1, 0, 0, 3};
    ibz_t value;
    ibz_init(&value);
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); ++i) {
        ibz_set(&value, values[i]);
        require(ibz_mod_ui(&value, 4) == residues[i], "incorrect signed residue");
    }
    ibz_set(&value, 1);
    ibz_mul_2exp(&value, &value, 96);
    ibz_add(&value, &value, &ibz_const_one);
    ibz_neg(&value, &value);
    require(ibz_mod_ui(&value, 4) == 3, "incorrect multi-limb negative residue");
    ibz_finalize(&value);
}

int
main(void)
{
    static const struct {
        const char *name;
        uint64_t values[4]; /* x,y,z,t before the production parity swap */
        int wide, accepted;
    } fixtures[] = {
        {"both negative differences", {1, 2, 4, 3}, 0, 1},
        {"first negative difference", {1, 4, 2, 3}, 0, 1},
        {"second negative difference", {3, 2, 4, 1}, 0, 1},
        {"first subtraction underflows int32", {UINT64_C(2147483648), 3, 1, 2}, 0, 1},
        {"first subtraction overflows int32", {UINT64_C(2147483647), 2, 4, UINT64_C(2147483649)}, 0, 1},
        {"second subtraction underflows int32", {3, UINT64_C(2147483648), 2, 1}, 0, 1},
        {"second subtraction overflows int32", {2, UINT64_C(2147483647), UINT64_C(2147483649), 4}, 0, 1},
        {"multi-limb negative difference", {1, 2, 4, 3}, 1, 1},
        {"parity swap", {2, 1, 4, 3}, 0, 1},
        {"reject first congruence", {1, 2, 4, 1}, 0, 0},
        {"reject second congruence", {2, 1, 1, 4}, 0, 0},
    };
    ibz_vec_4_init(&candidate);
    ibz_init(&cornacchia_target);
    check_negative_residues();
    for (size_t i = 0; i < sizeof(fixtures) / sizeof(fixtures[0]); ++i)
        check_candidate(fixtures[i].name, fixtures[i].values,
                        fixtures[i].wide, fixtures[i].accepted);
    ibz_finalize(&cornacchia_target);
    ibz_vec_4_finalize(&candidate);
    puts("normeq regression: PASS");
    return EXIT_SUCCESS;
}
