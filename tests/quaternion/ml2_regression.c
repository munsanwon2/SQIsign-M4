/* SPDX-License-Identifier: Apache-2.0 */
#include "lll_internals.h"
#include "mlll_internals.h"
#include "dpe.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void
require(int ok, const char *message)
{
    if (!ok) {
        fprintf(stderr, "MLLL regression: %s\n", message);
        exit(EXIT_FAILURE);
    }
}

/* Instrument divisions in the actual ML2 source, including calls made by
 * the public MLLL wrapper.  Catch a bad intermediate even if a retry would
 * later hide it.  No production source or DPE implementation is replaced. */
static void
checked_dpe_div(dpe_t out, dpe_t numerator, dpe_t denominator)
{
    require(isfinite(DPE_MANT(numerator)) &&
            isfinite(DPE_MANT(denominator)) && DPE_MANT(denominator) > 0.0,
            "nonpositive or nonfinite Gram--Schmidt divisor");
    dpe_div(out, numerator, denominator);
    require(isfinite(DPE_MANT(out)), "nonfinite Gram--Schmidt quotient");
}

#define dpe_div checked_dpe_div
#include "lll/ml2.c"
#undef dpe_div

static void
check_gram_update(void)
{
    static const int coordinates[3][4] = {{1, -2, 3, 1}, {3, 0, -1, 2}, {7, -5, 2, 4}};
    static const int coefficients[2] = {3, -2};
    static const int metric[4] = {1, 1, 17, 17};
    int gram[3][3] = {{0}}, products[3] = {0};
    quat_alg_t alg;
    quat_alg_init_set_ui(&alg, 17);
    ml2_state_t state;
    ml2_state_init(&state, 3, &alg);
    for (int i = 0; i < 3; ++i) {
        for (int c = 0; c < 4; ++c)
            ibz_set(&state.b[i][c], coordinates[i][c]);
        for (int j = 0; j < 3; ++j) {
            for (int c = 0; c < 4; ++c)
                gram[i][j] += metric[c] * coordinates[i][c] * coordinates[j][c];
        }
    }
    for (int i = 0; i < 2; ++i) {
        ibz_set(&state.coefficients[i], coefficients[i]);
        for (int j = 0; j < 3; ++j)
            products[j] += coefficients[i] * gram[i][j];
    }
    require(ml2_apply_size_reduction(&state, 2), "Babai exact vector update");
    for (int j = 0; j < 3; ++j) {
        const int expected = j < 2 ? gram[2][j] - products[j] : gram[2][2] - 2 * products[2] +
            coefficients[0] * products[0] + coefficients[1] * products[1];
        ibz_t recomputed;
        require(ml2_gram_entry(&recomputed, &state, 2, j) && ibz_cmp_int32(&recomputed, expected) == 0,
                "on-demand Gram and exact snapshot update differ");
    }
    ml2_state_finalize(&state);
    quat_alg_finalize(&alg);
}

/* Bareiss elimination gives each Gram-Schmidt pivot and coefficient as
 * a ratio of exact integer minors. The small fixture limits bound every
 * product here well within 128 bits. */
static void
check_reduced(const ibz_mat_4x4_t *basis, int rank, unsigned p)
{
    __int128 gram[4][4] = {{0}};
    require(p <= 367, "fixture metric exceeds the exact test bound");
    for (int i = 0; i < rank; ++i) {
        for (int j = 0; j < rank; ++j) {
            for (int c = 0; c < 4; ++c) {
                require(ibz_bitsize(&(*basis)[c][i]) <= 5 && ibz_bitsize(&(*basis)[c][j]) <= 5,
                        "fixture output exceeds the exact test bound");
                gram[i][j] += (__int128)(c < 2 ? 1 : p) * ibz_get(&(*basis)[c][i]) *
                              ibz_get(&(*basis)[c][j]);
            }
        }
    }
    __int128 previous = 1;
    for (int k = 0; k < rank; ++k) {
        const __int128 pivot = gram[k][k];
        require(pivot > 0, "output basis is dependent");
        for (int i = k + 1; i < rank; ++i) {
            const __int128 value = gram[i][k];
            require(20 * (value < 0 ? -value : value) <= 11 * pivot,
                    "output exceeds the exact size-reduction bound");
            for (int j = k + 1; j < rank; ++j) {
                const __int128 numerator = pivot * gram[i][j] - value * gram[k][j];
                require(numerator % previous == 0, "non-exact Gram minor division");
                gram[i][j] = numerator / previous;
            }
        }
        if (k > 0) {
            const __int128 lambda = gram[k][k - 1];
            const __int128 earlier = k == 1 ? 1 : gram[k - 2][k - 2];
            require(43 * previous * previous <= 50 * (pivot * earlier + lambda * lambda),
                    "output violates the exact Lovasz inequality");
        }
        previous = pivot;
    }
}

static void
check_span(const ibz_mat_4x4_t *basis, int rank, int identity, unsigned p)
{
    ibz_mat_4x4_t completed;
    ibz_t determinant, difference;
    ibz_mat_4x4_init(&completed);
    ibz_init(&determinant);
    ibz_init(&difference);
    for (int col = 0; col < rank; ++col) {
        for (int row = 0; row < 4; ++row)
            ibz_copy(&completed[row][col], &(*basis)[row][col]);
        if (!identity) {
            /* Input lattice: x0 == x1 (mod 2), x2 == 0 (mod 2),
             * and x3 == 0 at rank three.  These are membership checks. */
            ibz_sub(&difference, &(*basis)[0][col], &(*basis)[1][col]);
            require(ibz_is_even(&difference) && ibz_is_even(&(*basis)[2][col]),
                    "output left the input integer lattice");
            if (rank == 3)
                require(ibz_is_zero(&(*basis)[3][col]),
                        "output left the rank-three coordinate space");
        }
    }
    for (int col = rank; col < 4; ++col)
        ibz_set(&completed[col][col], 1);
    ibz_mat_4x4_inv_with_det_as_denom(NULL, &determinant, &completed);
    ibz_abs(&determinant, &determinant);
    /* Membership plus the same index proves equality over Z, not just Q. */
    require(ibz_cmp_int32(&determinant, identity ? 1 : 4) == 0,
            "output changed the input lattice index");
    check_reduced(basis, rank, p);
    ibz_finalize(&difference);
    ibz_finalize(&determinant);
    ibz_mat_4x4_finalize(&completed);
}

static void
check_dimension(int d, unsigned p)
{
    ibz_vec_4_t generators[16], original[16], output[4];
    ibz_mat_4x4_t basis;
    quat_alg_t alg;
    const int expected_rank = d <= 5 ? 3 : 4;
    const int start = d == 4 ? 0 : 1;
    quat_alg_init_set_ui(&alg, p);
    ibz_mat_4x4_init(&basis);
    for (int i = 0; i < d; ++i) {
        ibz_vec_4_init(&generators[i]);
        ibz_vec_4_init(&original[i]);
    }
    for (int i = 0; i < 4; ++i) {
        ibz_vec_4_init(&output[i]);
        for (int c = 0; c < 4; ++c)
            ibz_set(&output[i][c], 7001 + 4 * i + c);
    }

    /* Minimal eager-CFA failure: [2e0,2e1,2e2,e0+e1].  Inserting
     * e0+e1 first makes the next two old rows dependent; the old eager
     * sweep divides by that zero diagonal while rebuilding the last row.
     * Larger fixtures add an initial zero, e3, and more dependencies. */
    ibz_set(&generators[start][0], 2);
    ibz_set(&generators[start + 1][1], 2);
    ibz_set(&generators[start + 2][2], 2);
    ibz_set(&generators[start + 3][0], 1);
    ibz_set(&generators[start + 3][1], 1);
    if (d > 5) {
        ibz_set(&generators[5][3], 1);
        for (int i = 6; i < d - 1; ++i) {
            ibz_set(&generators[i][0], (i & 1) ? -1 : 1);
            ibz_set(&generators[i][1], 1);
            ibz_set(&generators[i][2], 2);
        }
        /* The last generator stays zero. */
    }
    for (int i = 0; i < d; ++i)
        ibz_vec_4_copy(&original[i], &generators[i]);

    int rank = quat_ml2(output, 4, generators, d, &alg);
    require(rank == expected_rank, "incorrect ML2 rank");
    for (int col = 0; col < rank; ++col)
        for (int row = 0; row < 4; ++row)
            ibz_copy(&basis[row][col], &output[col][row]);
    check_span(&basis, rank, 0, p);
    for (int i = rank; i < 4; ++i)
        for (int c = 0; c < 4; ++c)
            require(ibz_cmp_int32(&output[i][c], 7001 + 4 * i + c) == 0,
                    "ML2 wrote beyond the returned rank");

    quat_mlll(&basis, &rank, generators, d, &alg);
    require(rank == expected_rank, "incorrect public MLLL rank");
    check_span(&basis, rank, 0, p);
    for (int col = rank; col < 4; ++col)
        for (int row = 0; row < 4; ++row)
            require(ibz_is_zero(&basis[row][col]), "MLLL left a nonzero unused column");
    for (int i = 0; i < d; ++i)
        for (int c = 0; c < 4; ++c)
            require(ibz_cmp(&generators[i][c], &original[i][c]) == 0,
                    "reduction changed its input");

    /* All-zero generators: no factorization or output publication. */
    for (int i = 0; i < d; ++i)
        for (int c = 0; c < 4; ++c)
            ibz_set(&generators[i][c], 0);
    for (int i = 0; i < 4; ++i)
        for (int c = 0; c < 4; ++c)
            ibz_set(&output[i][c], 9001 + 4 * i + c);
    require(quat_ml2(output, 4, generators, d, &alg) == 0,
            "all-zero ML2 input has nonzero rank");
    for (int i = 0; i < 4; ++i)
        for (int c = 0; c < 4; ++c)
            require(ibz_cmp_int32(&output[i][c], 9001 + 4 * i + c) == 0,
                    "rank-zero ML2 call modified output");
    quat_mlll(&basis, &rank, generators, d, &alg);
    require(rank == 0, "all-zero MLLL input has nonzero rank");
    for (int row = 0; row < 4; ++row)
        for (int col = 0; col < 4; ++col)
            require(ibz_is_zero(&basis[row][col]), "rank-zero MLLL output is not zero");

    /* Ordinary full-rank control with the same dimensions and extra zeros. */
    for (int i = 0; i < 4; ++i)
        ibz_set(&generators[i][i], 1);
    quat_mlll(&basis, &rank, generators, d, &alg);
    require(rank == 4, "full-rank control lost rank");
    check_span(&basis, rank, 1, p);

    for (int i = 0; i < 4; ++i)
        ibz_vec_4_finalize(&output[i]);
    for (int i = 0; i < d; ++i) {
        ibz_vec_4_finalize(&original[i]);
        ibz_vec_4_finalize(&generators[i]);
    }
    ibz_mat_4x4_finalize(&basis);
    quat_alg_finalize(&alg);
}

int
main(void)
{
    static const int dimensions[] = {4, 5, 8, 16};
    static const unsigned weights[] = {7, 367};
    for (size_t i = 0; i < sizeof(dimensions) / sizeof(dimensions[0]); ++i)
        for (size_t j = 0; j < sizeof(weights) / sizeof(weights[0]); ++j)
            check_dimension(dimensions[i], weights[j]);
    check_gram_update();
    puts("MLLL regression: PASS");
    return EXIT_SUCCESS;
}
