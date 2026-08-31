/**
 * @file mlll.c
 * @brief Compact modified LLL entry points backed by ML2 (NS09)
 *
 * Implementation of the MLLL algorithm that takes a generating set
 * (possibly linearly dependent) and produces an LLL-reduced basis.
 *
 * The paper calls the NS09 floating-point modified LLL routine ML2/MLLL.
 * Exact vectors and Gram entries remain ibz_t; dpe_t is used only for the
 * approximate reduction decisions inside ml2.c.
 *
 * Paper mapping ("Compact Quaternion Algorithms for SQIsign"):
 *   Algorithm 1 (MLLL) — entry point quat_mlll
 *   Algorithm 2 line 6-7 — quat_lattice_mul_mlll
 *   Lemma 3 (mlll-bound) — bound applies to b/h/G integer entries.
 */

#include <quaternion.h>
#include "internal.h"
#include "lll_internals.h"
#include "mlll_internals.h"

#define N MLLL_MAX_GENERATORS

static int
mlll_lattice_add_products_fit(const quat_lattice_t *lat1,
                              const quat_lattice_t *lat2)
{
    if (lat1 == NULL || lat2 == NULL || ibz_is_zero(&lat1->denom) ||
        ibz_is_zero(&lat2->denom))
        return 0;

    int basis1_bits = 0;
    int basis2_bits = 0;
    for (int row = 0; row < 4; row++) {
        for (int column = 0; column < 4; column++) {
            int bits1 = ibz_bitsize(&lat1->basis[row][column]);
            int bits2 = ibz_bitsize(&lat2->basis[row][column]);
            if (bits1 > basis1_bits)
                basis1_bits = bits1;
            if (bits2 > basis2_bits)
                basis2_bits = bits2;
        }
    }
    int product1_bits = ibz_bitsize(&lat1->denom) + basis2_bits;
    int product2_bits = ibz_bitsize(&lat2->denom) + basis1_bits;
    int generator_bits = product1_bits > product2_bits ? product1_bits : product2_bits;
    return product1_bits <= IBZ_BITS - 1 &&
           product2_bits <= IBZ_BITS - 1 &&
           ibz_bitsize(&lat1->denom) + ibz_bitsize(&lat2->denom) <= IBZ_BITS - 1 &&
           2 * generator_bits + 2 <= IBZ_BITS - 1;
}

static int
mlll_positive_sum_bound(int left, int right)
{
    if (left == 0)
        return right;
    if (right == 0)
        return left;
    return (left > right ? left : right) + 1;
}

static int
mlll_product_bound(const ibz_t *left, const ibz_t *right)
{
    int left_bits = ibz_bitsize(left);
    int right_bits = ibz_bitsize(right);
    return left_bits == 0 || right_bits == 0 ? 0 : left_bits + right_bits;
}

static int
mlll_weight_bound(int bound, const quat_alg_t *alg)
{
    if (bound == 0 || ibz_is_one(&alg->p))
        return bound;
    return bound + ibz_bitsize(&alg->p);
}

/* Bound every transient in one quat_alg_coord_mul exactly according to the
 * source evaluation order.  In particular, p is charged only when one of the
 * coordinate-2/3 products is nonzero, avoiding rejection of large purely
 * unweighted coordinates. */
static int
mlll_quaternion_product_fits(const ibz_vec_4_t *left,
                             const ibz_vec_4_t *right,
                             const quat_alg_t *alg)
{
    int bound0 = mlll_positive_sum_bound(
        mlll_product_bound(&(*left)[2], &(*right)[2]),
        mlll_product_bound(&(*left)[3], &(*right)[3]));
    bound0 = mlll_weight_bound(bound0, alg);
    bound0 = mlll_positive_sum_bound(
        bound0, mlll_product_bound(&(*left)[0], &(*right)[0]));
    bound0 = mlll_positive_sum_bound(
        bound0, mlll_product_bound(&(*left)[1], &(*right)[1]));

    int bound1 = mlll_positive_sum_bound(
        mlll_product_bound(&(*left)[2], &(*right)[3]),
        mlll_product_bound(&(*left)[3], &(*right)[2]));
    bound1 = mlll_weight_bound(bound1, alg);
    bound1 = mlll_positive_sum_bound(
        bound1, mlll_product_bound(&(*left)[0], &(*right)[1]));
    bound1 = mlll_positive_sum_bound(
        bound1, mlll_product_bound(&(*left)[1], &(*right)[0]));

    int bound2 = mlll_positive_sum_bound(
        mlll_product_bound(&(*left)[0], &(*right)[2]),
        mlll_product_bound(&(*left)[2], &(*right)[0]));
    bound2 = mlll_positive_sum_bound(
        bound2, mlll_product_bound(&(*left)[1], &(*right)[3]));
    bound2 = mlll_positive_sum_bound(
        bound2, mlll_product_bound(&(*left)[3], &(*right)[1]));

    int bound3 = mlll_positive_sum_bound(
        mlll_product_bound(&(*left)[0], &(*right)[3]),
        mlll_product_bound(&(*left)[3], &(*right)[0]));
    bound3 = mlll_positive_sum_bound(
        bound3, mlll_product_bound(&(*left)[2], &(*right)[1]));
    bound3 = mlll_positive_sum_bound(
        bound3, mlll_product_bound(&(*left)[1], &(*right)[2]));

    return bound0 <= IBZ_BITS - 1 && bound1 <= IBZ_BITS - 1 &&
           bound2 <= IBZ_BITS - 1 && bound3 <= IBZ_BITS - 1;
}

static int
mlll_lattice_mul_products_fit(const quat_lattice_t *lat1,
                              const quat_lattice_t *lat2,
                              const quat_alg_t *alg)
{
    if (lat1 == NULL || lat2 == NULL || alg == NULL ||
        ibz_is_zero(&lat1->denom) || ibz_is_zero(&lat2->denom) ||
        ibz_cmp(&alg->p, &ibz_const_zero) <= 0)
        return 0;

    if (ibz_bitsize(&lat1->denom) + ibz_bitsize(&lat2->denom) >
        IBZ_BITS - 1)
        return 0;

    for (int left_column = 0; left_column < 4; left_column++) {
        ibz_vec_4_t left;
        ibz_vec_4_init(&left);
        for (int row = 0; row < 4; row++)
            ibz_copy(&left[row], &lat1->basis[row][left_column]);
        for (int right_column = 0; right_column < 4; right_column++) {
            ibz_vec_4_t right;
            int fits;
            ibz_vec_4_init(&right);
            for (int row = 0; row < 4; row++)
                ibz_copy(&right[row], &lat2->basis[row][right_column]);
            fits = mlll_quaternion_product_fits(&left, &right, alg);
            ibz_vec_4_finalize(&right);
            if (!fits) {
                ibz_vec_4_finalize(&left);
                return 0;
            }
        }
        ibz_vec_4_finalize(&left);
    }
    return 1;
}

static int
mlll_div_exact(ibz_t *quotient, const ibz_t *numerator, const ibz_t *denominator)
{
    ibz_t remainder;
    int exact;

    if (ibz_is_zero(denominator)) {
        ibz_set(quotient, 0);
        return 0;
    }
    ibz_init(&remainder);
    ibz_div(quotient, &remainder, numerator, denominator);
    exact = ibz_is_zero(&remainder);
    if (!exact)
        ibz_set(quotient, 0);
    ibz_finalize(&remainder);
    return exact;
}

/* Bound the exact numerator of trd(alpha * conjugate(beta)) before the
 * denominator division.  In coordinates this is
 * 2*(a0*b0 + a1*b1 + p*a2*b2 + p*a3*b3). */
static int
mlll_trace_pairing_fits(const ibz_vec_4_t *left,
                        const ibz_vec_4_t *right,
                        const quat_alg_t *alg)
{
    int bound = 0;
    for (int coordinate = 0; coordinate < 4; coordinate++) {
        int term = mlll_product_bound(
            &(*left)[coordinate], &(*right)[coordinate]);
        if (coordinate >= 2)
            term = mlll_weight_bound(term, alg);
        bound = mlll_positive_sum_bound(bound, term);
    }
    return bound == 0 || bound + 1 <= IBZ_BITS - 1;
}

static int
mlll_reduced_trace_pairing(ibz_t *trace,
                           const ibz_vec_4_t *left,
                           const ibz_vec_4_t *right,
                           const ibz_t *denominator,
                           const quat_alg_t *alg)
{
    ibz_t numerator, term;
    int ok = 0;

    if (!mlll_trace_pairing_fits(left, right, alg) ||
        ibz_is_zero(denominator))
        return 0;

    ibz_init(&numerator);
    ibz_init(&term);
    ibz_set(&numerator, 0);
    for (int coordinate = 0; coordinate < 4; coordinate++) {
        ibz_mul(&term, &(*left)[coordinate], &(*right)[coordinate]);
        if (coordinate >= 2)
            ibz_mul(&term, &term, &alg->p);
        ibz_add(&numerator, &numerator, &term);
    }
    ibz_mul(&numerator, &numerator, &ibz_const_two);
    ok = mlll_div_exact(trace, &numerator, denominator);
    ibz_finalize(&term);
    ibz_finalize(&numerator);
    return ok;
}

static int
mlll_scaled_coordinate_fits(const ibz_t *coordinate,
                            const ibz_t *scale,
                            const ibz_t *ideal_scalar)
{
    int coordinate_bits = ibz_bitsize(coordinate);
    if (coordinate_bits == 0)
        return 1;
    return coordinate_bits + ibz_bitsize(scale) +
               ibz_bitsize(ideal_scalar) <=
           IBZ_BITS - 1;
}

/* ========== integer helpers ========== */

/* The former Cohen fraction-free port below cannot represent a zero
 * Gram--Schmidt determinant followed by an independent vector: after a
 * dependent generator is swapped, it attempts an exact division by that zero
 * determinant.  Keep it excluded as implementation history; the compact
 * specification identifies MLLL with the NS09 ML2 reducer implemented in
 * ml2.c, which is used by the active entry point below. */
#if 0

/* Quaternion norm bilinear form: <a,b> = a0*b0 + a1*b1 + p*a2*b2 + p*a3*b3 */
static void
ibz_vec_4_dot_quat(ibz_t *dot, const ibz_vec_4_t *a, const ibz_vec_4_t *b,
                   const ibz_t *p)
{
    ibz_t tmp;
    ibz_init(&tmp);
    ibz_set(dot, 0);
    for (int i = 0; i < 4; i++) {
        ibz_mul(&tmp, &((*a)[i]), &((*b)[i]));
        if (i >= 2)
            ibz_mul(&tmp, &tmp, p);
        ibz_add(dot, dot, &tmp);
    }
    ibz_finalize(&tmp);
}

/* b = b - t * a  (integer vectors, t integer) */
static void
ibz_vec_4_sub_scalar_mul(ibz_vec_4_t *b, const ibz_t *t, const ibz_vec_4_t *a)
{
    ibz_t tmp;
    ibz_init(&tmp);
    for (int i = 0; i < 4; i++) {
        ibz_mul(&tmp, t, &((*a)[i]));
        ibz_sub(&((*b)[i]), &((*b)[i]), &tmp);
    }
    ibz_finalize(&tmp);
}

static void
ibz_vec_4_swap(ibz_vec_4_t *a, ibz_vec_4_t *b)
{
    for (int i = 0; i < 4; i++)
        ibz_swap(&((*a)[i]), &((*b)[i]));
}

/* Exact integer division n / d.
 * Used for Cohen 2.6.7 GS / swap updates where the division is provably exact.
 * Return failure instead of silently accepting a non-zero remainder: with
 * fixed-width arithmetic that condition can also signal an earlier overflow. */
static int
ibz_div_exact(ibz_t *q, const ibz_t *n, const ibz_t *d)
{
    if (ibz_is_zero(d)) {
        ibz_set(q, 0);
        return 0;
    }
    ibz_t r;
    ibz_init(&r);
    ibz_div_floor(q, &r, n, d);
    int exact = ibz_is_zero(&r);
    if (!exact)
        ibz_set(q, 0);
    ibz_finalize(&r);
    return exact;
}

/* Round-to-nearest integer division: q = round(n / d), ties away from zero.
 * Implementation: q = floor((2n + sign(d)*|d|) / (2d)) for n,d positive,
 * generalized via signed shifting. */
static void
ibz_div_round(ibz_t *q, const ibz_t *n, const ibz_t *d)
{
    ibz_t two_n, two_d, abs_d, n_adj, r;
    ibz_init(&two_n);
    ibz_init(&two_d);
    ibz_init(&abs_d);
    ibz_init(&n_adj);
    ibz_init(&r);

    /* two_n = 2*n, abs_d = |d|, two_d = 2*d */
    ibz_mul(&two_n, n, &((ibz_t *)&ibz_const_two)[0]);
    ibz_abs(&abs_d, d);
    ibz_mul(&two_d, d, &((ibz_t *)&ibz_const_two)[0]);

    /* n_adj = 2n + sign(n*d) * |d|  -- ties away from zero */
    int sign_match = (ibz_cmp(n, (const ibz_t *)&ibz_const_zero) >= 0)
                     == (ibz_cmp(d, (const ibz_t *)&ibz_const_zero) >= 0);
    if (sign_match) {
        ibz_add(&n_adj, &two_n, &abs_d);
    } else {
        ibz_sub(&n_adj, &two_n, &abs_d);
    }
    ibz_div(q, &r, &n_adj, &two_d);

    ibz_finalize(&two_n);
    ibz_finalize(&two_d);
    ibz_finalize(&abs_d);
    ibz_finalize(&n_adj);
    ibz_finalize(&r);
}

/* Test 2*|Lambda[k][l]| > d[l+1]  (size-reduce trigger condition). */
static int
needs_size_reduce(const ibz_t *Lambda_kl, const ibz_t *d_l1)
{
    ibz_t two_abs, abs_d;
    ibz_init(&two_abs);
    ibz_init(&abs_d);
    ibz_abs(&two_abs, Lambda_kl);
    ibz_mul(&two_abs, &two_abs, &((ibz_t *)&ibz_const_two)[0]);
    ibz_abs(&abs_d, d_l1);
    int res = (ibz_cmp(&two_abs, &abs_d) > 0);
    ibz_finalize(&two_abs);
    ibz_finalize(&abs_d);
    return res;
}

/* ========== Cohen 2.6.6 GS computation for a single position ========== */

/**
 * @brief Compute Lambda[idx][j] for j < idx and d[idx+1].
 *        Assumes positions 0..idx-1 are already computed.
 *
 * Cohen 2.6.6 integer-only GSO:
 *   For j = 0..idx-1:
 *     u = <b_idx, b_j>
 *     For i = 0..j-1:
 *       if d[i+1] != 0:
 *         u = (d[i+1]*u - Lambda[idx][i]*Lambda[j][i]) / d[i]   (exact)
 *     Lambda[idx][j] = u
 *
 *   u = <b_idx, b_idx>
 *   For j = 0..idx-1:
 *     if d[j+1] != 0:
 *       u = (d[j+1]*u - Lambda[idx][j]^2) / d[j]   (exact)
 *   d[idx+1] = u
 */
static int
compute_gs_single(int idx,
                  const ibz_vec_4_t *b,
                  ibz_t d[],
                  ibz_t Lambda[][N],
                  const ibz_t *p)
{
    ibz_t u, dot, t1, t2;
    ibz_init(&u);
    ibz_init(&dot);
    ibz_init(&t1);
    ibz_init(&t2);

    /* Lambda[idx][j] for j = 0..idx-1 */
    for (int j = 0; j < idx; j++) {
        ibz_vec_4_dot_quat(&dot, &b[idx], &b[j], p);
        ibz_copy(&u, &dot);
        for (int i = 0; i < j; i++) {
            if (ibz_is_zero(&d[i + 1])) {
                /* Skip: b_i* = 0, contributes nothing */
                continue;
            }
            ibz_mul(&t1, &d[i + 1], &u);
            ibz_mul(&t2, &Lambda[idx][i], &Lambda[j][i]);
            ibz_sub(&t1, &t1, &t2);
            if (!ibz_div_exact(&u, &t1, &d[i])) {
                goto failure;
            }
        }
        ibz_copy(&Lambda[idx][j], &u);
    }

    /* d[idx+1] */
    ibz_vec_4_dot_quat(&dot, &b[idx], &b[idx], p);
    ibz_copy(&u, &dot);
    for (int j = 0; j < idx; j++) {
        if (ibz_is_zero(&d[j + 1]))
            continue;
        ibz_mul(&t1, &d[j + 1], &u);
        ibz_mul(&t2, &Lambda[idx][j], &Lambda[idx][j]);
        ibz_sub(&t1, &t1, &t2);
        if (!ibz_div_exact(&u, &t1, &d[j])) {
            goto failure;
        }
    }
    ibz_copy(&d[idx + 1], &u);

    ibz_finalize(&u);
    ibz_finalize(&dot);
    ibz_finalize(&t1);
    ibz_finalize(&t2);
    return 1;

failure:
    ibz_finalize(&u);
    ibz_finalize(&dot);
    ibz_finalize(&t1);
    ibz_finalize(&t2);
    return 0;
}

/* ========== size-reduce step (Cohen 2.6.7 step 4) ========== */

/**
 * @brief Reduce(m, l): if 2|Lambda[m][l]| > d[l+1], subtract q*b_l from b_m
 *        and update Lambda[m][*] accordingly.
 *
 *   q = round(Lambda[m][l] / d[l+1])
 *   b[m] -= q * b[l]
 *   Lambda[m][l] -= q * d[l+1]
 *   for i = 0..l-1: Lambda[m][i] -= q * Lambda[l][i]
 */
static void
size_reduce_one(int m, int l,
                ibz_vec_4_t b[],
                ibz_t d[],
                ibz_t Lambda[][N])
{
    if (ibz_is_zero(&d[l + 1]))
        return;
    if (!needs_size_reduce(&Lambda[m][l], &d[l + 1]))
        return;

    ibz_t q, prod;
    ibz_init(&q);
    ibz_init(&prod);

    ibz_div_round(&q, &Lambda[m][l], &d[l + 1]);

    /* b[m] -= q * b[l] */
    ibz_vec_4_sub_scalar_mul(&b[m], &q, &b[l]);

    /* Lambda[m][l] -= q * d[l+1] */
    ibz_mul(&prod, &q, &d[l + 1]);
    ibz_sub(&Lambda[m][l], &Lambda[m][l], &prod);

    /* Lambda[m][i] -= q * Lambda[l][i] for i = 0..l-1 */
    for (int i = 0; i < l; i++) {
        ibz_mul(&prod, &q, &Lambda[l][i]);
        ibz_sub(&Lambda[m][i], &Lambda[m][i], &prod);
    }

    ibz_finalize(&q);
    ibz_finalize(&prod);
}

/* ========== Lovász condition (Cohen 2.6.7 step 3) ==========
 *
 * Standard form: ||b_m*||^2 + mu_{m,m-1}^2 ||b_{m-1}*||^2 < (3/4) ||b_{m-1}*||^2
 *
 * Integer form (multiply both sides by 4 * d[m-1] * d[m]):
 *   4 * d[m+1] * d[m-1] + 4 * Lambda[m][m-1]^2 / d[m]  <  3 * d[m]^2 ?
 * Cleaner integer-only form (Cohen 2.6.7):
 *   4 * d[m+1] * d[m-1]  <  3 * d[m]^2 - 4 * Lambda[m][m-1]^2
 *
 * Returns 1 if swap needed, 0 otherwise.
 *
 * Special case: if d[m+1] == 0 (b_m dependent on b_0..b_{m-1}),
 * always swap so b_m moves left and gets eliminated.
 */
static int
lovasz_test(int m, const ibz_t d[], const ibz_t Lambda[][N])
{
    if (m <= 0)
        return 0;
    if (ibz_is_zero(&d[m]))
        return 0; /* b_{m-1}* = 0, no meaningful swap */

    if (ibz_is_zero(&d[m + 1])) {
        /* b_m is dependent — swap to eliminate */
        return 1;
    }

    ibz_t lhs, rhs, t1, four;
    ibz_init(&lhs);
    ibz_init(&rhs);
    ibz_init(&t1);
    ibz_init(&four);
    ibz_set(&four, 4);

    /* lhs = 4 * d[m+1] * d[m-1] */
    ibz_mul(&lhs, &d[m + 1], &d[m - 1]);
    ibz_mul(&lhs, &lhs, &four);

    /* rhs = 3 * d[m]^2 - 4 * Lambda[m][m-1]^2 */
    ibz_mul(&rhs, &d[m], &d[m]);
    ibz_mul(&rhs, &rhs, &((ibz_t *)&ibz_const_three)[0]);
    ibz_mul(&t1, &Lambda[m][m - 1], &Lambda[m][m - 1]);
    ibz_mul(&t1, &t1, &four);
    ibz_sub(&rhs, &rhs, &t1);

    int need_swap = (ibz_cmp(&lhs, &rhs) < 0);

    ibz_finalize(&lhs);
    ibz_finalize(&rhs);
    ibz_finalize(&t1);
    ibz_finalize(&four);
    return need_swap;
}

/* ========== Cohen 2.6.7 SWAP step ==========
 *
 *   lambda = Lambda[m][m-1]
 *   B = (d[m+1] * d[m-1] + lambda^2) / d[m]   (exact)
 *   for i = m+1..beta-1:
 *     t = Lambda[i][m]
 *     Lambda[i][m]   = (d[m+1] * Lambda[i][m-1] - lambda * t) / d[m]   (exact)
 *     Lambda[i][m-1] = (B * t + lambda * Lambda[i][m]) / d[m+1]        (exact)
 *   swap b[m] <-> b[m-1]
 *   for j = 0..m-2: swap Lambda[m][j] <-> Lambda[m-1][j]
 *   d[m] = B
 *
 * Special-case (Cohen 2.6.8 / Pohst MLLL):
 *   - d[m+1] = 0 (b_m dependent): cannot do exact div by d[m+1].
 *     Fallback: just swap b/Lambda, recompute GS from m-1.
 *   - d[m] = 0 (b_{m-1}* = 0): caller filtered, won't arrive here.
 */
static int
swap_step(int m, int beta,
          ibz_vec_4_t b[],
          ibz_t d[],
          ibz_t Lambda[][N],
          const ibz_t *p)
{
    if (ibz_is_zero(&d[m + 1])) {
        /* Dependent path: swap and recompute GS from m-1.
         * After enough swaps + size_reduce, b_m becomes zero and is removed. */
        ibz_vec_4_swap(&b[m], &b[m - 1]);
        for (int j = 0; j < m - 1; j++)
            ibz_swap(&Lambda[m][j], &Lambda[m - 1][j]);
        /* Recompute GS rows m-1..beta-1. d[m] may become 0 too. */
        for (int i = (m - 1 < 0 ? 0 : m - 1); i < beta; i++)
            if (!compute_gs_single(i, b, d, Lambda, p))
                return 0;
        return 1;
    }

    ibz_t lambda, B, t, t1, t2, num;
    ibz_init(&lambda);
    ibz_init(&B);
    ibz_init(&t);
    ibz_init(&t1);
    ibz_init(&t2);
    ibz_init(&num);

    ibz_copy(&lambda, &Lambda[m][m - 1]);

    /* B = (d[m+1] * d[m-1] + lambda^2) / d[m] */
    ibz_mul(&t1, &d[m + 1], &d[m - 1]);
    ibz_mul(&t2, &lambda, &lambda);
    ibz_add(&num, &t1, &t2);
    if (!ibz_div_exact(&B, &num, &d[m]))
        goto failure;

    /* Update Lambda[i][m], Lambda[i][m-1] for i = m+1..beta-1 */
    for (int i = m + 1; i < beta; i++) {
        ibz_copy(&t, &Lambda[i][m]);

        /* Lambda[i][m] = (d[m+1] * Lambda[i][m-1] - lambda * t) / d[m] */
        ibz_mul(&t1, &d[m + 1], &Lambda[i][m - 1]);
        ibz_mul(&t2, &lambda, &t);
        ibz_sub(&num, &t1, &t2);
        if (!ibz_div_exact(&Lambda[i][m], &num, &d[m]))
            goto failure;

        /* Lambda[i][m-1] = (B * t + lambda * Lambda[i][m]) / d[m+1] */
        ibz_mul(&t1, &B, &t);
        ibz_mul(&t2, &lambda, &Lambda[i][m]);
        ibz_add(&num, &t1, &t2);
        if (!ibz_div_exact(&Lambda[i][m - 1], &num, &d[m + 1]))
            goto failure;
    }

    /* Swap b[m] <-> b[m-1] */
    ibz_vec_4_swap(&b[m], &b[m - 1]);
    /* Swap Lambda[m][j] <-> Lambda[m-1][j] for j = 0..m-2 */
    for (int j = 0; j < m - 1; j++)
        ibz_swap(&Lambda[m][j], &Lambda[m - 1][j]);

    /* d[m] = B  (d[m-1] and d[m+1] unchanged) */
    ibz_copy(&d[m], &B);

    ibz_finalize(&lambda);
    ibz_finalize(&B);
    ibz_finalize(&t);
    ibz_finalize(&t1);
    ibz_finalize(&t2);
    ibz_finalize(&num);
    return 1;

failure:
    ibz_finalize(&lambda);
    ibz_finalize(&B);
    ibz_finalize(&t);
    ibz_finalize(&t1);
    ibz_finalize(&t2);
    ibz_finalize(&num);
    return 0;
}

/* ========== MLLL core (Pohst Algorithm 1, 0-indexed) ========== */

void
quat_mlll(ibz_mat_4x4_t *basis,
          int *rank,
          const ibz_vec_4_t *generators,
          int g,
          const quat_alg_t *alg)
{
    if (g < 1 || g > MLLL_MAX_GENERATORS) {
        *rank = 0;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                ibz_set(&((*basis)[i][j]), 0);
        return;
    }

    ibz_vec_4_t b[N];
    ibz_t d[N + 1];
    ibz_t Lambda[N][N];

    for (int i = 0; i < g; i++)
        ibz_vec_4_init(&b[i]);
    for (int i = 0; i <= g; i++)
        ibz_init(&d[i]);
    for (int i = 0; i < g; i++)
        for (int j = 0; j < g; j++)
            ibz_init(&Lambda[i][j]);

    /* d[0] = 1 */
    ibz_set(&d[0], 1);

    int alpha = 0; /* next generator to load */
    int beta = 0;  /* count of working vectors */
    int tau = 1;   /* lower bound on m for reduction restart */
    int m;

    /* ===== LOAD ===== */
load:
    while (alpha < g) {
        if (ibz_vec_4_is_zero(&generators[alpha])) {
            alpha++;
            continue;
        }
        for (int i = 0; i < 4; i++)
            ibz_copy(&(b[beta][i]), &(generators[alpha][i]));
        alpha++;

        if (!compute_gs_single(beta, b, d, Lambda, &alg->p))
            goto failure;
        beta++;

        if (!ibz_is_zero(&d[beta]) && alpha < g)
            continue;
        break;
    }

    if (beta <= 1)
        goto done;

    m = tau;
    if (m >= beta) {
        if (alpha < g) {
            tau = beta;
            goto load;
        }
        goto done;
    }

    /* ===== REDUCTION ===== */
reduction:
    /* Size-reduce b[m] against b[m-1], then test Lovász. */

    /* Cohen 2.6.7 ordering: size_reduce against m-1 only first, then test. */
    size_reduce_one(m, m - 1, b, d, Lambda);

    /* If b[m] became zero, remove it. */
    if (ibz_vec_4_is_zero(&b[m]))
        goto remove_vector;

    /* Lovász test */
    if (lovasz_test(m, d, Lambda)) {
        if (!swap_step(m, beta, b, d, Lambda, &alg->p))
            goto failure;
        if (m > 1)
            m--;
        goto reduction;
    }

    /* Size-reduce b[m] against b[l] for l = m-2 down to 0 */
    for (int l = m - 2; l >= 0; l--)
        size_reduce_one(m, l, b, d, Lambda);

    if (ibz_vec_4_is_zero(&b[m]))
        goto remove_vector;

    m++;
    if (m >= beta) {
        if (alpha < g) {
            tau = beta;
            goto load;
        }
        goto done;
    }
    goto reduction;

    /* ===== REMOVE ===== */
remove_vector:
    /* Shift down from position m */
    for (int i = m; i < beta - 1; i++) {
        for (int c = 0; c < 4; c++)
            ibz_copy(&(b[i][c]), &(b[i + 1][c]));
    }
    /* Clear top vector */
    for (int c = 0; c < 4; c++)
        ibz_set(&(b[beta - 1][c]), 0);
    beta--;

    if (alpha >= g) {
        /* Recompute all GS for final pass */
        ibz_set(&d[0], 1);
        for (int i = 0; i < beta; i++)
            if (!compute_gs_single(i, b, d, Lambda, &alg->p))
                goto failure;
        goto done;
    }

    /* Recompute GS from position m */
    for (int i = m; i < beta; i++)
        if (!compute_gs_single(i, b, d, Lambda, &alg->p))
            goto failure;

    tau = m + 1;
    if (tau < 1)
        tau = 1;
    goto load;

done:
    /* Final full LLL pass: ensures positions 0..tau-1 (potentially unreduced
     * after removals) are also LLL-reduced. */
    if (beta > 1) {
        ibz_set(&d[0], 1);
        for (int i = 0; i < beta; i++)
            if (!compute_gs_single(i, b, d, Lambda, &alg->p))
                goto failure;

        int changed = 1;
        int safety = 0;
        const int SAFETY_MAX = 1000 * beta * beta;
        while (changed && safety++ < SAFETY_MAX) {
            changed = 0;
            int mm = 1;
            while (mm < beta) {
                /* Size-reduce b[mm] against b[l] for l = mm-1 down to 0 */
                int reduced_any = 0;
                for (int l = mm - 1; l >= 0; l--) {
                    if (ibz_is_zero(&d[l + 1]))
                        continue;
                    if (needs_size_reduce(&Lambda[mm][l], &d[l + 1])) {
                        size_reduce_one(mm, l, b, d, Lambda);
                        reduced_any = 1;
                    }
                }
                if (reduced_any)
                    changed = 1;

                if (ibz_vec_4_is_zero(&b[mm])) {
                    for (int i = mm; i < beta - 1; i++)
                        for (int c = 0; c < 4; c++)
                            ibz_copy(&(b[i][c]), &(b[i + 1][c]));
                    for (int c = 0; c < 4; c++)
                        ibz_set(&(b[beta - 1][c]), 0);
                    beta--;
                    ibz_set(&d[0], 1);
                    for (int i = 0; i < beta; i++)
                        if (!compute_gs_single(i, b, d, Lambda, &alg->p))
                            goto failure;
                    changed = 1;
                    continue;
                }

                if (lovasz_test(mm, d, Lambda)) {
                    if (!swap_step(mm, beta, b, d, Lambda, &alg->p))
                        goto failure;
                    if (mm > 1)
                        mm--;
                    changed = 1;
                    continue;
                }
                mm++;
            }
        }
        if (changed)
            goto failure;
    }

    goto publish;

failure:
    /* Do not publish a partially reduced lattice after an arithmetic or
     * convergence failure.  Callers already use rank as the status channel. */
    beta = 0;

publish:
    /* Copy output: basis columns = non-zero b vectors */
    *rank = 0;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ibz_set(&((*basis)[i][j]), 0);
    for (int i = 0; i < beta && *rank < 4; i++) {
        if (!ibz_vec_4_is_zero(&b[i])) {
            for (int j = 0; j < 4; j++)
                ibz_copy(&((*basis)[j][*rank]), &(b[i][j]));
            (*rank)++;
        }
    }

    /* Cleanup */
    for (int i = 0; i < g; i++)
        ibz_vec_4_finalize(&b[i]);
    for (int i = 0; i <= g; i++)
        ibz_finalize(&d[i]);
    for (int i = 0; i < g; i++)
        for (int j = 0; j < g; j++)
            ibz_finalize(&Lambda[i][j]);
}

#endif

/* ML2 returns its four basis vectors as an array of rows.  Lattice matrices
 * store basis vectors in columns.  Writing directly into private candidate
 * matrix storage and transposing in place avoids a second four-vector bank. */
static void
mlll_transpose_basis_in_place(ibz_mat_4x4_t *basis)
{
    for (int row = 0; row < 4; row++)
        for (int column = row + 1; column < 4; column++)
            ibz_swap(&(*basis)[row][column], &(*basis)[column][row]);
}

#if defined(TARGET_ARM) && (defined(__GNUC__) || defined(__clang__))
#define MLLL_TRANSIENT_NOINLINE __attribute__((noinline))
#else
#define MLLL_TRANSIENT_NOINLINE
#endif

/* Keep quaternion-product temporaries out of the caller frame that remains
 * live across ML2.  quat_alg_coord_mul computes into a private sum before
 * publishing, so using `right` as both its second input and output is safe. */
static MLLL_TRANSIENT_NOINLINE void
mlll_form_ordered_products(ibz_vec_4_t generators[16],
                           const quat_lattice_t *left_lattice,
                           const quat_lattice_t *right_lattice,
                           const quat_alg_t *alg)
{
    ibz_vec_4_t left, right;

    ibz_vec_4_init(&left);
    ibz_vec_4_init(&right);
    for (int left_column = 0; left_column < 4; left_column++) {
        ibz_vec_4_copy_ibz(&left,
                           &left_lattice->basis[0][left_column],
                           &left_lattice->basis[1][left_column],
                           &left_lattice->basis[2][left_column],
                           &left_lattice->basis[3][left_column]);
        for (int right_column = 0; right_column < 4; right_column++) {
            ibz_vec_4_copy_ibz(&right,
                               &right_lattice->basis[0][right_column],
                               &right_lattice->basis[1][right_column],
                               &right_lattice->basis[2][right_column],
                               &right_lattice->basis[3][right_column]);
            quat_alg_coord_mul(&right, &left, &right, alg);
            for (int coordinate = 0; coordinate < 4; coordinate++)
                ibz_copy(&generators[4 * left_column + right_column][coordinate],
                         &right[coordinate]);
        }
    }
    ibz_vec_4_finalize(&right);
    ibz_vec_4_finalize(&left);
}

/* Active compact MLLL entry point.  ML2 publishes its output only after the
 * reduction finishes.  A failed attempt is retried using span-preserving
 * permutations; no surrounding protocol state is regenerated. */
void
quat_mlll(ibz_mat_4x4_t *basis,
          int *rank,
          const ibz_vec_4_t *generators,
          int g,
          const quat_alg_t *alg)
{
    int rho = -1;

    if (basis == NULL || rank == NULL)
        return;
    *rank = 0;
    for (int row = 0; row < 4; row++)
        for (int column = 0; column < 4; column++)
            ibz_set(&(*basis)[row][column], 0);

    if (generators == NULL || g < 1 || g > MLLL_MAX_GENERATORS)
        return;

    rho = quat_ml2_mlll(&(*basis)[0], 4, generators, g, alg);

    if (rho >= 0 && rho <= 4) {
        mlll_transpose_basis_in_place(basis);
        *rank = rho;
    }
}

/* Full-rank ideal operations must not accept a lower-rank first ML2 attempt.
 * Retry the fixed span-preserving permutations, then publish only rank four. */
static int
mlll_full_rank_basis(ibz_mat_4x4_t *basis,
                     const ibz_vec_4_t *generators,
                     int generator_count,
                     const quat_alg_t *alg)
{
    int rank;

    rank = quat_ml2_retry(
        &(*basis)[0], 4, generators, generator_count, alg);
    if (rank == 4)
        mlll_transpose_basis_in_place(basis);
    return rank == 4;
}

/* ========== Lattice operations using MLLL ========== */

int
quat_lattice_mul_mlll(quat_lattice_t *res,
                      const quat_lattice_t *lat1,
                      const quat_lattice_t *lat2,
                      const quat_alg_t *alg)
{
    ibz_vec_4_t generators[16];
    quat_lattice_t lat1_red, lat2_red;
    int ok = 0;

    if (res == NULL || lat1 == NULL || lat2 == NULL || alg == NULL ||
        ibz_is_zero(&lat1->denom) || ibz_is_zero(&lat2->denom) ||
        ibz_cmp(&alg->p, &ibz_const_zero) <= 0)
        return 0;

    quat_lattice_init(&lat1_red);
    quat_lattice_init(&lat2_red);
    for (int i = 0; i < 16; i++)
        ibz_vec_4_init(&generators[i]);

    /* Algorithm 2, Lines 3-4: shorten both input bases before forming the
     * sixteen ordered products. */
    if (!quat_lattice_lll(&lat1_red.basis, lat1, alg) ||
        !quat_lattice_lll(&lat2_red.basis, lat2, alg))
        goto cleanup;
    ibz_copy(&lat1_red.denom, &lat1->denom);
    ibz_copy(&lat2_red.denom, &lat2->denom);
    if (!mlll_lattice_mul_products_fit(&lat1_red, &lat2_red, alg))
        goto cleanup;

    mlll_form_ordered_products(generators, &lat1_red, &lat2_red, alg);

    /* The shortened left basis is dead after product formation.  Reuse its
     * private storage for the unpublished ML2 candidate. */
    if (!mlll_full_rank_basis(&lat1_red.basis, generators, 16, alg))
        goto cleanup;
    ibz_mul(&lat1_red.denom, &lat1_red.denom, &lat2_red.denom);
    if (!quat_lattice_reduce_denom(&lat1_red, &lat1_red))
        goto cleanup;
    ibz_mat_4x4_copy(&res->basis, &lat1_red.basis);
    ibz_copy(&res->denom, &lat1_red.denom);
    ok = 1;

cleanup:
    quat_lattice_finalize(&lat2_red);
    quat_lattice_finalize(&lat1_red);
    for (int i = 0; i < 16; i++)
        ibz_vec_4_finalize(&generators[i]);
    return ok;
}

int
quat_lattice_intersect_mlll(quat_lattice_t *res,
                            ibz_t *intersection_norm,
                            const quat_lattice_t *lat1,
                            const ibz_t *norm1,
                            const quat_lattice_t *lat2,
                            const ibz_t *norm2,
                            const quat_alg_t *alg)
{
    /* Paper Algorithm 3, CompactIdealIntersection.
     *
     * For integral left O_0-ideals, the trace gcd below is nrd(I1 + I2).
     * With a = nrd(I1)/d and b = nrd(I2)/d, Lemma 14 gives
     * I1 intersect I2 = b*I1 + a*I2. ML2 reduces those eight generators
     * directly, avoiding both lattice duals and HNF. */
    quat_lattice_t lat1_red, lat2_red;
    ibz_vec_4_t generators[8];
    ibz_vec_4_t left, right;
    ibz_t denominator1, denominator2, denominator_product;
    ibz_t denominator_gcd, common_denominator, scale1, scale2;
    ibz_t trace_gcd, trace, a, b, candidate_norm, quotient_scratch;
    int ok = 0;

    if (res == NULL || intersection_norm == NULL || lat1 == NULL ||
        norm1 == NULL || lat2 == NULL ||
        norm2 == NULL || alg == NULL || ibz_is_zero(&lat1->denom) ||
        ibz_is_zero(&lat2->denom) ||
        ibz_cmp(norm1, &ibz_const_zero) <= 0 ||
        ibz_cmp(norm2, &ibz_const_zero) <= 0 ||
        ibz_cmp(&alg->p, &ibz_const_zero) <= 0)
        return 0;

    quat_lattice_init(&lat1_red);
    quat_lattice_init(&lat2_red);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_init(&generators[i]);
    ibz_vec_4_init(&left);
    ibz_vec_4_init(&right);
    ibz_init(&denominator1);
    ibz_init(&denominator2);
    ibz_init(&denominator_product);
    ibz_init(&denominator_gcd);
    ibz_init(&common_denominator);
    ibz_init(&scale1);
    ibz_init(&scale2);
    ibz_init(&trace_gcd);
    ibz_init(&trace);
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&candidate_norm);
    ibz_init(&quotient_scratch);

    /* Algorithm 3, Lines 1-2. */
    if (!quat_lattice_lll(&(lat1_red.basis), lat1, alg))
        goto cleanup;
    if (!quat_lattice_lll(&(lat2_red.basis), lat2, alg))
        goto cleanup;
    ibz_copy(&(lat1_red.denom), &(lat1->denom));
    ibz_copy(&(lat2_red.denom), &(lat2->denom));

    ibz_abs(&denominator1, &lat1_red.denom);
    ibz_abs(&denominator2, &lat2_red.denom);
    if (ibz_bitsize(&denominator1) + ibz_bitsize(&denominator2) >
        IBZ_BITS - 1)
        goto cleanup;
    ibz_mul(&denominator_product, &denominator1, &denominator2);

    /* Algorithm 3, Lines 3-4: d = gcd(N1, N2, all trace pairings). */
    ibz_gcd(&trace_gcd, norm1, norm2);
    for (int r = 0; r < 4; r++) {
        for (int coordinate = 0; coordinate < 4; coordinate++)
            ibz_copy(&left[coordinate], &lat1_red.basis[coordinate][r]);
        for (int s = 0; s < 4; s++) {
            for (int coordinate = 0; coordinate < 4; coordinate++)
                ibz_copy(&right[coordinate], &lat2_red.basis[coordinate][s]);
            if (!mlll_reduced_trace_pairing(
                    &trace, &left, &right, &denominator_product, alg))
                goto cleanup;
            ibz_gcd(&trace_gcd, &trace_gcd, &trace);
        }
    }
    ibz_abs(&trace_gcd, &trace_gcd);
    if (ibz_is_zero(&trace_gcd) ||
        !mlll_div_exact(&a, norm1, &trace_gcd) ||
        !mlll_div_exact(&b, norm2, &trace_gcd))
        goto cleanup;
    if (mlll_product_bound(norm1, &b) > IBZ_BITS - 1)
        goto cleanup;
    ibz_mul(&candidate_norm, norm1, &b);

    /* Put both scaled ideals over one common denominator. */
    ibz_gcd(&denominator_gcd, &denominator1, &denominator2);
    if (ibz_is_zero(&denominator_gcd) ||
        !mlll_div_exact(&scale1, &denominator2, &denominator_gcd) ||
        !mlll_div_exact(&scale2, &denominator1, &denominator_gcd))
        goto cleanup;
    ibz_mul(&common_denominator, &denominator1, &scale1);

    /* Algorithm 3, Lines 5-6: ML2(b*alpha_1,...,b*alpha_4,
     *                              a*beta_1,...,a*beta_4). */
    for (int column = 0; column < 4; column++) {
        for (int row = 0; row < 4; row++) {
            if (!mlll_scaled_coordinate_fits(
                    &lat1_red.basis[row][column], &scale1, &b) ||
                !mlll_scaled_coordinate_fits(
                    &lat2_red.basis[row][column], &scale2, &a))
                goto cleanup;
            ibz_mul(&quotient_scratch, &lat1_red.basis[row][column], &scale1);
            ibz_mul(&generators[column][row], &quotient_scratch, &b);
            ibz_mul(&quotient_scratch, &lat2_red.basis[row][column], &scale2);
            ibz_mul(&generators[4 + column][row], &quotient_scratch, &a);
        }
    }
    /* Both scaled generator sets are complete, so the old left basis is no
     * longer live and can hold the unpublished intersection candidate. */
    if (!mlll_full_rank_basis(&lat1_red.basis, generators, 8, alg))
        goto cleanup;
    ibz_copy(&lat1_red.denom, &common_denominator);
    if (!quat_lattice_reduce_denom(&lat1_red, &lat1_red))
        goto cleanup;

    /* Publish only after every exact division and reduction succeeded. */
    ibz_mat_4x4_copy(&res->basis, &lat1_red.basis);
    ibz_copy(&res->denom, &lat1_red.denom);
    ibz_copy(intersection_norm, &candidate_norm);
    ok = 1;

cleanup:
    ibz_finalize(&quotient_scratch);
    ibz_finalize(&candidate_norm);
    ibz_finalize(&b);
    ibz_finalize(&a);
    ibz_finalize(&trace);
    ibz_finalize(&trace_gcd);
    ibz_finalize(&scale2);
    ibz_finalize(&scale1);
    ibz_finalize(&common_denominator);
    ibz_finalize(&denominator_gcd);
    ibz_finalize(&denominator_product);
    ibz_finalize(&denominator2);
    ibz_finalize(&denominator1);
    ibz_vec_4_finalize(&right);
    ibz_vec_4_finalize(&left);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_finalize(&generators[i]);
    quat_lattice_finalize(&lat1_red);
    quat_lattice_finalize(&lat2_red);
    return ok;
}

int
quat_lattice_add_mlll_with_reducer(quat_lattice_t *res,
                                   const quat_lattice_t *lat1,
                                   const quat_lattice_t *lat2,
                                   const quat_alg_t *alg,
                                   quat_ml2_reducer_t reducer)
{
    /* paper 03Ideal.tex: "we use the modified LLL algorithm with
     * floating-point arithmetic [NS09], which we call ML2 or modified LLL"
     * "We also call LLL as the special case of MLLL when d=4 in Algorithm ML2"
    * So paper's MLLL == ML2 (NS09 Fig 9, 53-bit float). Production uses
     * quat_ml2_retry so a failed ordering is retried by pure permutations. */
    ibz_vec_4_t generators[8];
    quat_lattice_t candidate;
    int ok = 0;

    if (res == NULL || alg == NULL ||
        ibz_cmp(&alg->p, &ibz_const_zero) <= 0 || reducer == NULL ||
        !mlll_lattice_add_products_fit(lat1, lat2))
        return 0;

    for (int i = 0; i < 8; i++)
        ibz_vec_4_init(&generators[i]);
    quat_lattice_init(&candidate);

    /* candidate is private until rank four and denominator normalization have
     * both succeeded.  Reuse its basis while materializing the two scaled
     * generator sets instead of keeping a second 4-by-4 temporary. */
    ibz_mat_4x4_scalar_mul(
        &candidate.basis, &(lat1->denom), &(lat2->basis));
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ibz_copy(&(generators[j][i]), &(candidate.basis[i][j]));

    ibz_mat_4x4_scalar_mul(
        &candidate.basis, &(lat2->denom), &(lat1->basis));
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ibz_copy(&(generators[4 + j][i]), &(candidate.basis[i][j]));

    int rho;
#if defined(TARGET_ARM)
    /* Alternate reducers are host-only fault-injection instrumentation.  The
     * embedded production path is a direct edge to the native retry driver,
     * keeping an indirect call out of conservative Cortex-M4 stack audits. */
    if (reducer != quat_ml2_retry)
        goto cleanup;
    rho = quat_ml2_retry(&candidate.basis[0], 4, generators, 8, alg);
#else
    rho = reducer(&candidate.basis[0], 4, generators, 8, alg);
#endif
    if (rho == 4) {
        mlll_transpose_basis_in_place(&candidate.basis);

        ibz_mul(&candidate.denom, &lat1->denom, &lat2->denom);
        if (quat_lattice_reduce_denom(&candidate, &candidate)) {
            /* The reducer may fail or return a non-full-rank result. Do not
             * publish any part of the candidate until rank four and a valid
             * denominator normalization are confirmed. */
            ibz_mat_4x4_copy(&res->basis, &candidate.basis);
            ibz_copy(&res->denom, &candidate.denom);
            ok = 1;
        }
    }

#if defined(TARGET_ARM)
cleanup:
#endif
    quat_lattice_finalize(&candidate);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_finalize(&generators[i]);
    return ok;
}

int
quat_lattice_add_mlll(quat_lattice_t *res,
                      const quat_lattice_t *lat1,
                      const quat_lattice_t *lat2,
                      const quat_alg_t *alg)
{
    return quat_lattice_add_mlll_with_reducer(res, lat1, lat2, alg, quat_ml2_retry);
}
