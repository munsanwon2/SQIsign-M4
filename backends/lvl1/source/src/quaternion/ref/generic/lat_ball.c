#include <quaternion.h>
#include <rng.h>
#include "internal.h"
#include "lll_internals.h"
#include "dpe.h"

static int
add_nonnegative_checked(ibz_t *sum, const ibz_t *a, const ibz_t *b)
{
    return !ibz_is_negative(a) && !ibz_is_negative(b) && ibz_add_checked(sum, a, b);
}

static int
mul_nonnegative_checked(ibz_t *product, const ibz_t *a, const ibz_t *b)
{
    return !ibz_is_negative(a) && !ibz_is_negative(b) && ibz_mul_checked(product, a, b);
}

static int
abs_checked(ibz_t *absolute, const ibz_t *value)
{
    ibz_abs(absolute, value);
    return !ibz_is_negative(absolute);
}

/* Bound every output coordinate of mat*[-box,box] (or mat^t*box), while
 * preflighting exactly the same products and sums used by matrix evaluation. */
static int
matrix_box_bound(ibz_vec_4_t *out,
                 const ibz_mat_4x4_t *mat,
                 const ibz_vec_4_t *box,
                 int transpose)
{
    ibz_t coefficient, product;
    ibz_init(&coefficient);
    ibz_init(&product);
    for (int i = 0; i < 4; i++) {
        ibz_set(&(*out)[i], 0);
        for (int j = 0; j < 4; j++) {
            const ibz_t *entry = transpose ? &(*mat)[j][i] : &(*mat)[i][j];
            if (!abs_checked(&coefficient, entry) ||
                !mul_nonnegative_checked(&product, &coefficient, &(*box)[j]) ||
                !add_nonnegative_checked(&(*out)[i], &(*out)[i], &product)) {
                ibz_finalize(&coefficient);
                ibz_finalize(&product);
                return 0;
            }
        }
    }
    ibz_finalize(&coefficient);
    ibz_finalize(&product);
    return 1;
}

static int
qf_box_fits(const ibz_mat_4x4_t *G, const ibz_vec_4_t *box)
{
    int ok = 0;
    ibz_vec_4_t row_bounds;
    ibz_t product, total;
    ibz_vec_4_init(&row_bounds);
    ibz_init(&product);
    ibz_init(&total);
    if (!matrix_box_bound(&row_bounds, G, box, 0))
        goto cleanup;
    for (int i = 0; i < 4; i++) {
        if (!mul_nonnegative_checked(&product, &row_bounds[i], &(*box)[i]) ||
            !add_nonnegative_checked(&total, &total, &product))
            goto cleanup;
    }
    ok = 1;
cleanup:
    ibz_finalize(&total);
    ibz_finalize(&product);
    ibz_vec_4_finalize(&row_bounds);
    return ok;
}

/* quat_lattice_gram forms two ordinary products, optionally multiplies by p,
 * sums four terms, and doubles.  Reproduce absolute upper bounds first so a
 * Release build cannot wrap before the reduced-basis certificate sees G. */
static int
lattice_gram_fits(const quat_lattice_t *lattice, const quat_alg_t *alg)
{
    ibz_t a, b, p, product, sum;
    ibz_init(&a);
    ibz_init(&b);
    ibz_init(&p);
    ibz_init(&product);
    ibz_init(&sum);
    if (!abs_checked(&p, &alg->p))
        goto failure;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&sum, 0);
            for (int k = 0; k < 4; k++) {
                if (!abs_checked(&a, &lattice->basis[k][i]) ||
                    !abs_checked(&b, &lattice->basis[k][j]) ||
                    !mul_nonnegative_checked(&product, &a, &b) ||
                    (k >= 2 && !mul_nonnegative_checked(&product, &product, &p)) ||
                    !add_nonnegative_checked(&sum, &sum, &product))
                    goto failure;
            }
            if (!add_nonnegative_checked(&sum, &sum, &sum))
                goto failure;
        }
    }
    ibz_finalize(&a); ibz_finalize(&b); ibz_finalize(&p);
    ibz_finalize(&product); ibz_finalize(&sum);
    return 1;
failure:
    ibz_finalize(&a); ibz_finalize(&b); ibz_finalize(&p);
    ibz_finalize(&product); ibz_finalize(&sum);
    return 0;
}

/* Compute floor(c * numerator / denominator), where 1 <= c <= 5, without
 * ever forming c*numerator.  The modular-accumulation loop also avoids the
 * otherwise possible overflow in c*remainder when denominator is full width. */
static int
floor_ratio_small(ibz_t *quotient, const ibz_t *numerator, const ibz_t *denominator, unsigned c)
{
    int ok = 0;
    unsigned carry = 0;
    ibz_t q, r, acc, gap, scaled_q, small;
    ibz_init(&q);
    ibz_init(&r);
    ibz_init(&acc);
    ibz_init(&gap);
    ibz_init(&scaled_q);
    ibz_init(&small);

    if (c == 0 || c > 5 || ibz_is_negative(numerator) ||
        ibz_cmp(denominator, &ibz_const_zero) <= 0)
        goto cleanup;

    ibz_div(&q, &r, numerator, denominator);
    /* c <= 5.  This conservative bit check guarantees the signed product
     * fits even in builds where arithmetic overflow checks are disabled. */
    if (!ibz_is_zero(&q)) {
        const unsigned headroom = (c == 1) ? 0 : (c == 2) ? 1 : (c <= 4) ? 2 : 3;
        if (ibz_bitsize(&q) > IBZ_BITS - 1 - (int)headroom)
            goto cleanup;
    }
    ibz_set(&small, (int32_t)c);
    ibz_mul(&scaled_q, &q, &small);

    for (unsigned i = 0; i < c; i++) {
        ibz_sub(&gap, denominator, &acc);
        if (ibz_cmp(&r, &gap) >= 0) {
            ibz_sub(&acc, &r, &gap);
            carry++;
        } else {
            ibz_add(&acc, &acc, &r); /* acc+r < denominator: cannot overflow */
        }
    }
    ibz_set(&small, (int32_t)carry);
    if (!add_nonnegative_checked(quotient, &scaled_q, &small))
        goto cleanup;
    ok = 1;

cleanup:
    ibz_finalize(&q);
    ibz_finalize(&r);
    ibz_finalize(&acc);
    ibz_finalize(&gap);
    ibz_finalize(&scaled_q);
    ibz_finalize(&small);
    return ok;
}

/* Return ceil((17/32)*x).  17/32 is a dyadic, rigorously certified upper
 * bound for every Gram--Schmidt coefficient accepted below. */
static int
ceil_eta_times(ibz_t *out, const ibz_t *x)
{
    int ok = 0;
    int32_t rem_i;
    ibz_t denominator, numerator, q, r, tail;
    ibz_init(&denominator);
    ibz_init(&numerator);
    ibz_init(&q);
    ibz_init(&r);
    ibz_init(&tail);
    ibz_set(&denominator, 32);
    ibz_set(&numerator, 17);

    if (ibz_is_negative(x))
        goto cleanup;
    ibz_div(&q, &r, x, &denominator);
    ibz_mul(out, &q, &numerator);
    rem_i = ibz_get(&r);
    ibz_set(&tail, (17 * rem_i + 31) / 32);
    if (!add_nonnegative_checked(out, out, &tail))
        goto cleanup;
    ok = 1;

cleanup:
    ibz_finalize(&denominator);
    ibz_finalize(&numerator);
    ibz_finalize(&q);
    ibz_finalize(&r);
    ibz_finalize(&tail);
    return ok;
}

/* A tiny outward interval layer for certifying the L2 postcondition.  DPE
 * values alone are approximations and therefore cannot justify a sampling
 * bound.  Each input interval contains the exact integer.  Each arithmetic
 * operation expands both endpoints by 2^-36 of the operation scale; IEEE-754
 * binary64 basic operations incur at most 2^-52 relative error, and DPE drops
 * an addend only below half an ulp.  The 16-bit guard therefore covers every
 * endpoint rounding (and normalization rounding) with a factor greater than
 * 2^15.  Interval dependency can make certification fail, but can never make
 * it accept a false inequality. */
#define DPE_CERT_GUARD_BITS 36

typedef struct
{
    dpe_t lo;
    dpe_t hi;
} dpe_interval_t;

static void
dpe_interval_set(dpe_interval_t *out, dpe_interval_t *in)
{
    dpe_set(out->lo, in->lo);
    dpe_set(out->hi, in->hi);
}

static int
dpe_nonzero_exp(const dpe_t x)
{
    return DPE_MANT(x) == 0.0 ? DPE_EXPMIN : DPE_EXP(x);
}

static int
dpe_interval_scale_exp(dpe_interval_t *x)
{
    int lo = dpe_nonzero_exp(x->lo);
    int hi = dpe_nonzero_exp(x->hi);
    return lo > hi ? lo : hi;
}

static void
dpe_interval_expand(dpe_interval_t *x, int scale_exp)
{
    dpe_t error;
    dpe_init(error);
    if (scale_exp == DPE_EXPMIN)
        return;
    DPE_MANT(error) = 0.5;
    DPE_EXP(error) = scale_exp - DPE_CERT_GUARD_BITS + 1;
    dpe_sub(x->lo, x->lo, error);
    dpe_add(x->hi, x->hi, error);
}

static void
dpe_interval_set_z(dpe_interval_t *out, const ibz_t *value)
{
    dpe_set_z(out->lo, value);
    dpe_set(out->hi, out->lo);
    if (!ibz_is_zero(value) && ibz_bitsize(value) > DPE_BITSIZE)
        /* dpe_set_z truncates by less than 2^(bits-53).  Expand by
         * 2^(bits-48), leaving another 5 bits for the endpoint update. */
        dpe_interval_expand(out, ibz_bitsize(value) - 12);
}

static void
dpe_interval_set_dyadic(dpe_interval_t *out, long numerator, unsigned shift)
{
    dpe_set_si(out->lo, numerator);
    dpe_div_2exp(out->lo, out->lo, shift);
    dpe_set(out->hi, out->lo);
}

static void
dpe_interval_sub(dpe_interval_t *out, dpe_interval_t *a, dpe_interval_t *b)
{
    int ea = dpe_interval_scale_exp(a);
    int eb = dpe_interval_scale_exp(b);
    dpe_sub(out->lo, a->lo, b->hi);
    dpe_sub(out->hi, a->hi, b->lo);
    dpe_interval_expand(out, ea > eb ? ea : eb);
}

static void
dpe_interval_mul(dpe_interval_t *out, dpe_interval_t *a, dpe_interval_t *b)
{
    dpe_t products[4];
    int ea = dpe_interval_scale_exp(a);
    int eb = dpe_interval_scale_exp(b);
    for (int i = 0; i < 4; i++)
        dpe_init(products[i]);
    dpe_mul(products[0], a->lo, b->lo);
    dpe_mul(products[1], a->lo, b->hi);
    dpe_mul(products[2], a->hi, b->lo);
    dpe_mul(products[3], a->hi, b->hi);
    dpe_set(out->lo, products[0]);
    dpe_set(out->hi, products[0]);
    for (int i = 1; i < 4; i++) {
        if (dpe_cmp(products[i], out->lo) < 0)
            dpe_set(out->lo, products[i]);
        if (dpe_cmp(products[i], out->hi) > 0)
            dpe_set(out->hi, products[i]);
    }
    if (ea != DPE_EXPMIN && eb != DPE_EXPMIN)
        dpe_interval_expand(out, ea + eb);
}

static int
dpe_interval_div(dpe_interval_t *out, dpe_interval_t *a, dpe_interval_t *b)
{
    dpe_t quotients[4];
    if (dpe_cmp_d(b->lo, 0.0) <= 0 && dpe_cmp_d(b->hi, 0.0) >= 0)
        return 0;
    for (int i = 0; i < 4; i++)
        dpe_init(quotients[i]);
    dpe_div(quotients[0], a->lo, b->lo);
    dpe_div(quotients[1], a->lo, b->hi);
    dpe_div(quotients[2], a->hi, b->lo);
    dpe_div(quotients[3], a->hi, b->hi);
    dpe_set(out->lo, quotients[0]);
    dpe_set(out->hi, quotients[0]);
    for (int i = 1; i < 4; i++) {
        if (dpe_cmp(quotients[i], out->lo) < 0)
            dpe_set(out->lo, quotients[i]);
        if (dpe_cmp(quotients[i], out->hi) > 0)
            dpe_set(out->hi, quotients[i]);
    }
    dpe_interval_expand(out, dpe_interval_scale_exp(out));
    return 1;
}

/* Certify positive definiteness, and for a reduced Gram matrix also certify
 * |mu_ij| <= 17/32 and r_i >= (19/32) r_(i-1).  These dyadic thresholds are
 * exact and leave a large gap to the Section 4 L2 targets. */
static int
certify_gram(const ibz_mat_4x4_t *G, int reduced, dpe_t shortest_lower)
{
    dpe_interval_t A[4][4], r[4], mu[4][4];
    dpe_interval_t numerator, product, temp, eta, neg_eta, alpha, rhs;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (ibz_cmp(&(*G)[i][j], &(*G)[j][i]) != 0)
                return 0;
            dpe_interval_set_z(&A[i][j], &(*G)[i][j]);
        }
    }
    dpe_interval_set_dyadic(&eta, 17, 5);
    dpe_interval_set_dyadic(&neg_eta, -17, 5);
    dpe_interval_set_dyadic(&alpha, 19, 5);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < i; j++) {
            dpe_interval_set(&numerator, &A[i][j]);
            for (int k = 0; k < j; k++) {
                dpe_interval_mul(&product, &mu[i][k], &mu[j][k]);
                dpe_interval_mul(&product, &product, &r[k]);
                dpe_interval_sub(&numerator, &numerator, &product);
            }
            if (!dpe_interval_div(&mu[i][j], &numerator, &r[j]))
                return 0;
            if (reduced &&
                (dpe_cmp(mu[i][j].lo, neg_eta.lo) < 0 ||
                 dpe_cmp(mu[i][j].hi, eta.hi) > 0))
                return 0;
        }

        dpe_interval_set(&r[i], &A[i][i]);
        for (int k = 0; k < i; k++) {
            dpe_interval_mul(&temp, &mu[i][k], &mu[i][k]);
            dpe_interval_mul(&product, &temp, &r[k]);
            dpe_interval_sub(&r[i], &r[i], &product);
        }
        if (dpe_cmp_d(r[i].lo, 0.0) <= 0)
            return 0;
        if (reduced && i > 0) {
            dpe_interval_mul(&rhs, &alpha, &r[i - 1]);
            if (dpe_cmp(r[i].lo, rhs.hi) < 0)
                return 0;
        }
    }
    if (shortest_lower != NULL) {
        dpe_set(shortest_lower, r[0].lo);
        for (int i = 1; i < 4; i++) {
            if (dpe_cmp(r[i].lo, shortest_lower) < 0)
                dpe_set(shortest_lower, r[i].lo);
        }
    }
    return 1;
}

/* An inconclusive interval is not a proof of non-positivity.  For small
 * entries, certify Sylvester's criterion exactly.  If |G_ij| < 2^b, each
 * determinant partial sum is bounded by 24*2^(4b) < 2^(4b+5); the bit gate
 * reserves the sign bit and also bounds all cofactors in the existing
 * inverse helper.  Large inputs keep the interval-only capacity policy. */
static int
certify_gram_small_exact(const ibz_mat_4x4_t *G)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (ibz_cmp(&(*G)[i][j], &(*G)[j][i]) != 0 ||
                ibz_bitsize(&(*G)[i][j]) > (IBZ_BITS - 6) / 4)
                return 0;
        }
    }

    int ok = 1;
    ibz_mat_4x4_t minor;
    ibz_t det;
    ibz_mat_4x4_init(&minor);
    ibz_init(&det);
    for (int n = 1; n <= 4; n++) {
        ibz_mat_4x4_identity(&minor);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                ibz_copy(&minor[i][j], &(*G)[i][j]);
        ibz_mat_4x4_inv_with_det_as_denom(NULL, &det, &minor);
        if (ibz_cmp(&det, &ibz_const_zero) <= 0) {
            ok = 0;
            break;
        }
    }
    ibz_finalize(&det);
    ibz_mat_4x4_finalize(&minor);
    return ok;
}

static int
bound_parallelogram_reduced(ibz_vec_4_t *box,
                            ibz_mat_4x4_t *U,
                            ibz_mat_4x4_t *reduced_G_out,
                            const ibz_mat_4x4_t *G,
                            const ibz_t *radius)
{
    if (ibz_cmp(radius, &ibz_const_zero) <= 0)
        return 0;

    int ok = 0;
    int nonzero = 0;
    static const unsigned diagonal_factors[4] = { 1, 2, 3, 4 };
    ibz_mat_4x4_t reduced_G, transform;
    ibz_vec_4_t base;
    ibz_t ratio, tail_sum, eta_tail;
    dpe_t shortest_lower;
    dpe_interval_t radius_interval;
    ibz_mat_4x4_init(&reduced_G);
    ibz_mat_4x4_init(&transform);
    ibz_vec_4_init(&base);
    ibz_init(&ratio);
    ibz_init(&tail_sum);
    ibz_init(&eta_tail);
    dpe_init(shortest_lower);

    if (!certify_gram(G, 0, NULL) && !certify_gram_small_exact(G))
        goto cleanup;

    /* Work in a primal L2-reduced basis.  If T is the transformation built
     * by quat_lll_core, then reduced_G=T^t*G*T.  The sampler below evaluates
     * U^t*y, so publish U=T^t after the bounds have been computed. */
    ibz_mat_4x4_copy(&reduced_G, G);
    ibz_mat_4x4_identity(&transform);
    if (!quat_lll_core_checked(&reduced_G, &transform) ||
        !certify_gram(&reduced_G, 1, shortest_lower))
        goto cleanup;

    /* Let r_i be the squared Gram--Schmidt lengths and mu_ij the reduced
     * coefficients.  The outward interval certificate above proves the
     * exact dyadic inequalities |mu_ij|<=17/32 and
     * r_i>=(19/32)*r_(i-1).  Consequently
     *
     *   G_ii = r_i + sum_{j<i} mu_ij^2 r_j <= C_i*r_i,
     *
     * with the integer upper bounds C=(1,2,3,4).  Thus every vector of norm
     * at most radius has |z_i|<=sqrt(C_i*radius/G_ii) in orthogonal
     * coordinates.  This replaces the former rounded inverse/cofactor
     * diagonal, which could underbound even when its DPE result was off by
     * only a few ulps after determinant cancellation. */
    for (int i = 0; i < 4; i++) {
        if (ibz_cmp(&reduced_G[i][i], &ibz_const_zero) <= 0 ||
            !floor_ratio_small(&ratio, radius, &reduced_G[i][i], diagonal_factors[i]))
            goto cleanup;
        ibz_sqrt_floor(&base[i], &ratio);
        nonzero |= !ibz_is_zero(&base[i]);
    }

    /* Back substitution gives
     *   |x_i| <= |z_i| + eta * sum_{j>i}|x_j|.
     * Since x_i is integral, floor(|z_i| bound)+ceil(eta*tail) is an
     * integer-outward bound: floor(t+a)<=floor(t)+ceil(a). */
    /* If every base bound is zero, R<G_ii/C_i<=r_i for every i; the highest
     * nonzero coefficient would therefore already contribute more than R.
     * Independently, R<min_i(r_i) proves the same fact.  In either case the
     * origin is the only integral point, and the sampler deliberately excludes
     * it, so fail immediately instead of performing 5,000,001 futile draws. */
    dpe_interval_set_z(&radius_interval, radius);
    if (!nonzero || dpe_cmp(radius_interval.hi, shortest_lower) < 0)
        goto cleanup;

    ibz_copy(&(*box)[3], &base[3]);
    for (int i = 2; i >= 0; i--) {
        ibz_set(&tail_sum, 0);
        for (int j = i + 1; j < 4; j++) {
            if (!add_nonnegative_checked(&tail_sum, &tail_sum, &(*box)[j]))
                goto cleanup;
        }
        if (!ceil_eta_times(&eta_tail, &tail_sum) ||
            !add_nonnegative_checked(&(*box)[i], &base[i], &eta_tail))
            goto cleanup;
    }

    /* ibz_rand_interval samples an interval of width 2*box[i], which must
     * itself remain a positive signed ibz_t. */
    for (int i = 0; i < 4; i++) {
        if (ibz_bitsize(&(*box)[i]) > IBZ_BITS - 2)
            goto cleanup;
    }

    ibz_mat_4x4_transpose(U, &transform);
    if (reduced_G_out != NULL)
        ibz_mat_4x4_copy(reduced_G_out, &reduced_G);
    ok = 1;

cleanup:
    ibz_finalize(&eta_tail);
    ibz_finalize(&tail_sum);
    ibz_finalize(&ratio);
    ibz_vec_4_finalize(&base);
    ibz_mat_4x4_finalize(&transform);
    ibz_mat_4x4_finalize(&reduced_G);
    return ok;
}

int
quat_lattice_bound_parallelogram(ibz_vec_4_t *box,
                                 ibz_mat_4x4_t *U,
                                 const ibz_mat_4x4_t *G,
                                 const ibz_t *radius)
{
    return bound_parallelogram_reduced(box, U, NULL, G, radius);
}

int
quat_lattice_sample_from_ball(quat_alg_elem_t *res,
                              const quat_lattice_t *lattice,
                              const quat_alg_t *alg,
                              const ibz_t *radius)
{
    if (res == NULL || lattice == NULL || alg == NULL || radius == NULL ||
        ibz_cmp(radius, &ibz_const_zero) <= 0 ||
        ibz_is_zero(&lattice->denom) ||
        ibz_cmp(&alg->p, &ibz_const_zero) <= 0)
        return 0;
    if (!lattice_gram_fits(lattice, alg))
        return 0;

    ibz_vec_4_t box;
    ibz_vec_4_init(&box);
    ibz_vec_4_t original_box, result_box;
    ibz_vec_4_init(&original_box);
    ibz_vec_4_init(&result_box);
    ibz_mat_4x4_t U, G, reduced_G;
    ibz_mat_4x4_init(&U);
    ibz_mat_4x4_init(&G);
    ibz_mat_4x4_init(&reduced_G);
    ibz_vec_4_t x;
    ibz_vec_4_init(&x);
    ibz_t rad, tmp, abs_denom;
    ibz_init(&rad);
    ibz_init(&tmp);
    ibz_init(&abs_denom);

    int ok = 0;

    // Compute the Gram matrix of the lattice (preflighted above).
    quat_lattice_gram(&G, lattice, alg);

    // Correct ball radius by the denominator and by 2 (the Gram matrix is
    // twice the norm), with exact capacity checks before every product.
    if (!abs_checked(&abs_denom, &lattice->denom) ||
        !mul_nonnegative_checked(&rad, radius, &abs_denom) ||
        !mul_nonnegative_checked(&rad, &rad, &abs_denom) ||
        !mul_nonnegative_checked(&rad, &rad, &ibz_const_two))
        goto err;

    // Compute an integer-outward bounding parallelogram in a reduced basis.
    // This avoids both the 4K-bit integer cofactor transient and the unsafe
    // inward rounding of the former floating inverse-diagonal estimate.
    ok = bound_parallelogram_reduced(&box, &U, &reduced_G, &G, &rad);
    if (!ok || !qf_box_fits(&reduced_G, &box) ||
        !matrix_box_bound(&original_box, &U, &box, 1) ||
        !matrix_box_bound(&result_box, &lattice->basis, &original_box, 0)) {
        ok = 0;
        goto err;
    }

    // Rejection sampling from the parallelogram
    int cnt = 0;
    do {
        // Sample vector
        for (int i = 0; i < 4; i++) {
            if (ibz_is_zero(&box[i])) {
                ibz_copy(&x[i], &ibz_const_zero);
            } else {
                ibz_add(&tmp, &box[i], &box[i]);
                ok &= ibz_rand_interval(&x[i], &ibz_const_zero, &tmp);
                ibz_sub(&x[i], &x[i], &box[i]);
                if (!ok)
                    goto err;
            }
        }
        // Evaluate in reduced coordinates.  qf_box_fits certified every
        // product and partial sum over this complete box.
        quat_qf_eval(&tmp, &reduced_G, &x);
        cnt++;
        if (cnt > 5000000) {
            ok = 0;
            goto err;
        }
    } while (ibz_is_zero(&tmp) || (ibz_cmp(&tmp, &rad) > 0));

    // Map the accepted reduced coordinate back only once.
    ibz_mat_4x4_eval_t(&x, &x, &U);

    // Evaluate linear combination
    ibz_mat_4x4_eval(&(res->coord), &(lattice->basis), &x);
    ibz_copy(&(res->denom), &(lattice->denom));
    quat_alg_normalize(res);

#ifndef NDEBUG
    // Check norm is smaller than radius
    quat_alg_norm(&tmp, &rad, res, alg);
    ibz_mul(&rad, &rad, radius);
    assert(ibz_cmp(&tmp, &rad) <= 0);
#endif

err:
    ibz_finalize(&rad);
    ibz_finalize(&tmp);
    ibz_finalize(&abs_denom);
    ibz_vec_4_finalize(&x);
    ibz_mat_4x4_finalize(&reduced_G);
    ibz_mat_4x4_finalize(&U);
    ibz_mat_4x4_finalize(&G);
    ibz_vec_4_finalize(&result_box);
    ibz_vec_4_finalize(&original_box);
    ibz_vec_4_finalize(&box);
    return ok;
}

quat_ball_status_t
quat_lattice_sample_from_ball_direct_with_sampler(quat_alg_elem_t *res,
                                                 const quat_lattice_t *lattice,
                                                 const quat_alg_t *alg,
                                                 const ibz_t *bound,
                                                 unsigned trials,
                                                 int (*sample)(ibz_t *, const ibz_t *, const ibz_t *))
{
    if (res == NULL || lattice == NULL || alg == NULL || bound == NULL || sample == NULL ||
        ibz_cmp(&lattice->denom, &ibz_const_zero) <= 0 ||
        ibz_cmp(&alg->p, &ibz_const_zero) <= 0 || ibz_cmp(bound, &ibz_const_zero) <= 0)
        return QUAT_BALL_INVALID;

    static const int kappa[4] = { 14, 6, 4, 5 };
    ibz_t scaled_bound, limit, weighted_limit, tmp, remainder;
    if (!ibz_mul_checked(&scaled_bound, &lattice->denom, &lattice->denom) ||
        !ibz_mul_checked(&scaled_bound, &scaled_bound, bound))
        return QUAT_BALL_PRECISION;
    ibz_sqrt_floor(&limit, &scaled_bound);
    ibz_div(&tmp, &remainder, &scaled_bound, &alg->p);
    ibz_sqrt_floor(&weighted_limit, &tmp);

    ibz_vec_4_t upper, lower, column;
    for (int i = 0; i < 4; i++) {
        ibz_t norm, factor;
        for (int row = 0; row < 4; row++)
            ibz_copy(&column[row], &lattice->basis[row][i]);
        if (!quat_alg_coord_norm(&norm, &column, alg))
            return QUAT_BALL_PRECISION;
        if (ibz_is_zero(&norm))
            return QUAT_BALL_INVALID;
        ibz_set(&factor, kappa[i]);
        if (!ibz_mul_checked(&tmp, &factor, &scaled_bound))
            return QUAT_BALL_PRECISION;
        ibz_div(&tmp, &remainder, &tmp, &norm);
        ibz_sqrt_floor(&upper[i], &tmp);
        ibz_neg(&lower[i], &upper[i]);
    }

    quat_alg_elem_t candidate;
    ibz_copy(&candidate.denom, &lattice->denom);
    for (unsigned trial = 0; trial < trials; trial++) {
        ibz_vec_4_t coefficients;
        for (int i = 0; i < 4; i++) {
            if (ibz_is_zero(&upper[i])) {
                ibz_set(&coefficients[i], 0);
            } else if (!sample(&coefficients[i], &lower[i], &upper[i])) {
                return QUAT_BALL_RNG_FAILURE;
            }
        }
        if (!ibz_mat_4x4_eval_checked(&candidate.coord, &lattice->basis, &coefficients))
            return QUAT_BALL_PRECISION;

        /* Rejected proposals can have a norm above the allocation. Screen
         * coordinates before forming even the first square. */
        int inside = 1;
        for (int i = 0; i < 4; i++) {
            ibz_abs(&tmp, &candidate.coord[i]);
            if (ibz_is_negative(&tmp) || ibz_cmp(&tmp, i < 2 ? &limit : &weighted_limit) > 0) {
                inside = 0;
                break;
            }
        }
        if (!inside)
            continue;

        /* Each surviving term is at most S. If their sum exceeds signed
         * capacity it is above S as well, so it is a rejected proposal. */
        if (quat_alg_coord_norm(&tmp, &candidate.coord, alg) && ibz_cmp(&tmp, &scaled_bound) <= 0) {
            quat_alg_elem_copy(res, &candidate);
            return QUAT_BALL_SUCCESS;
        }
    }
    return QUAT_BALL_EXHAUSTED;
}

quat_ball_status_t
quat_lattice_sample_from_ball_direct(quat_alg_elem_t *res,
                                    const quat_lattice_t *lattice,
                                    const quat_alg_t *alg,
                                    const ibz_t *bound)
{
    return quat_lattice_sample_from_ball_direct_with_sampler(res, lattice, alg, bound, 65536, ibz_rand_interval);
}
