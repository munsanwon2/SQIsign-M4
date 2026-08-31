#include <quaternion.h>
#include "internal.h"
#include <stdint.h>
// get shortest equivalent ideal

#define DEBUG_PRINTS 0

/* Keep phase-local fixed-precision temporaries out of their callers' stack
 * frames.  This matters on Cortex-M4, where every ibz_t has a compile-time
 * maximum size and a 4x4 matrix is consequently several kilobytes at level
 * V. */
#if defined(__GNUC__) || defined(__clang__)
#define QLAPOTI_NOINLINE __attribute__((noinline))
#else
#define QLAPOTI_NOINLINE
#endif

/* All probabilistic searches have a build-time ceiling in addition to their
 * caller-supplied budget.  The defaults are deliberately much larger than
 * the expected number of trials, yet finite on a microcontroller with a bad
 * entropy source or adversarial input. */
#ifndef SQISIGN_QLAPOTI_GENERATOR_MAX_ATTEMPTS
#define SQISIGN_QLAPOTI_GENERATOR_MAX_ATTEMPTS 4096U
#endif

#ifndef SQISIGN_QLAPOTI_LAMBDA_MAX_ATTEMPTS
#define SQISIGN_QLAPOTI_LAMBDA_MAX_ATTEMPTS 4096U
#endif

#ifndef SQISIGN_QLAPOTI_MAX_COUNTER_ALPHA
#define SQISIGN_QLAPOTI_MAX_COUNTER_ALPHA 65536U
#endif

#ifndef SQISIGN_QLAPOTI_COUNTER_MULTIPLIER
#define SQISIGN_QLAPOTI_COUNTER_MULTIPLIER 20U
#endif

#ifndef SQISIGN_QLAPOTI_MAX_PRIMALITY_ITERATIONS
#define SQISIGN_QLAPOTI_MAX_PRIMALITY_ITERATIONS 256
#endif

#if SQISIGN_QLAPOTI_GENERATOR_MAX_ATTEMPTS == 0
#error "SQISIGN_QLAPOTI_GENERATOR_MAX_ATTEMPTS must be positive"
#endif

#if SQISIGN_QLAPOTI_LAMBDA_MAX_ATTEMPTS == 0
#error "SQISIGN_QLAPOTI_LAMBDA_MAX_ATTEMPTS must be positive"
#endif

#if SQISIGN_QLAPOTI_MAX_COUNTER_ALPHA == 0
#error "SQISIGN_QLAPOTI_MAX_COUNTER_ALPHA must be positive"
#endif

#if SQISIGN_QLAPOTI_COUNTER_MULTIPLIER == 0
#error "SQISIGN_QLAPOTI_COUNTER_MULTIPLIER must be positive"
#endif

#if SQISIGN_QLAPOTI_MAX_PRIMALITY_ITERATIONS <= 0
#error "SQISIGN_QLAPOTI_MAX_PRIMALITY_ITERATIONS must be positive"
#endif

_Static_assert(SQISIGN_QLAPOTI_MAX_COUNTER_ALPHA <=
                   UINT32_MAX / SQISIGN_QLAPOTI_COUNTER_MULTIPLIER,
               "Qlapoti total attempt cap does not fit uint32_t");

#define QLAPOTI_MAX_CORNACCHIA_PRIMES 101U

/* Compute the short generator while the two 4x4 reduction matrices are the
 * only large Qlapoti-owned objects live.  The caller's element is private, so
 * writing it before success does not expose a partially computed public
 * output. */
static QLAPOTI_NOINLINE int
quat_lideal_shortest_generator(quat_alg_elem_t *new_alpha,
                               const quat_left_ideal_t *lideal,
                               const quat_alg_t *alg)
{
    ibz_mat_4x4_t gram, red;
    ibz_t n, remainder;
    int status = 0;

    ibz_mat_4x4_init(&gram);
    ibz_mat_4x4_init(&red);
    ibz_init(&n);
    ibz_init(&remainder);

    if (!quat_lideal_reduce_basis(&red, &gram, lideal, alg))
        goto cleanup;

    ibz_set(&new_alpha->coord[0], 1);
    ibz_set(&new_alpha->coord[1], 0);
    ibz_set(&new_alpha->coord[2], 0);
    ibz_set(&new_alpha->coord[3], 0);
    ibz_mat_4x4_eval(&new_alpha->coord, &red, &new_alpha->coord);
    ibz_copy(&new_alpha->denom, &lideal->lattice.denom);
    if (!quat_lattice_contains(NULL, &lideal->lattice, new_alpha) ||
        !quat_alg_norm(&n, &remainder, new_alpha, alg) ||
        !ibz_is_one(&remainder))
        goto cleanup;
    ibz_div(&n, &remainder, &n, &lideal->norm);
    if (!ibz_is_zero(&remainder))
        goto cleanup;

    status = 1;

cleanup:
    ibz_finalize(&remainder);
    ibz_finalize(&n);
    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&red);
    return status;
}

int
quat_lideal_shortest_equivalent(quat_left_ideal_t *equiv,
                                 quat_alg_elem_t *elem,
                                 const quat_left_ideal_t *lideal,
                                 const quat_alg_t *alg)
{
    quat_alg_elem_t new_alpha;
    ibz_t original_denom;
    int status = 0;

    if (equiv == NULL || elem == NULL || lideal == NULL || alg == NULL ||
        ibz_cmp(&lideal->norm, &ibz_const_zero) <= 0 ||
        ibz_is_zero(&lideal->lattice.denom))
        return 0;

    quat_alg_elem_init(&new_alpha);
    ibz_init(&original_denom);

    if (!quat_lideal_shortest_generator(&new_alpha, lideal, alg))
        goto cleanup;

    ibz_copy(&original_denom, &new_alpha.denom);
    ibz_mul(&new_alpha.denom, &new_alpha.denom, &lideal->norm);
    ibz_neg(&new_alpha.coord[0], &new_alpha.coord[0]);
    /* quat_lideal_mul is itself transactional and alias-safe.  A second
     * quat_left_ideal_t candidate here merely kept about 7 KiB live while
     * its own candidate was being constructed at level V. */
    if (!quat_lideal_mul(equiv, lideal, &new_alpha, alg))
        goto cleanup;

    ibz_neg(&new_alpha.coord[0], &new_alpha.coord[0]);
    ibz_copy(&new_alpha.denom, &original_denom);
    quat_alg_elem_copy(elem, &new_alpha);
    status = 1;

cleanup:
    ibz_finalize(&original_denom);
    quat_alg_elem_finalize(&new_alpha);
    return status;
}

int
quat_lideal_generator_small_coprime(quat_alg_elem_t *gen,
                                     const quat_left_ideal_t *lideal,
                                     const quat_alg_t *alg,
                                     int sampling_bound_bits,
                                     uint32_t max_attempts)
{
    int found = 0;
    int status = QUAT_QLAPOTI_RETRY;
    ibz_t n, d, n2, gcd;
    ibz_vec_4_t coeffs;
    quat_alg_elem_t candidate;

    if (gen == NULL || lideal == NULL || alg == NULL ||
        sampling_bound_bits <= 0 || sampling_bound_bits >= IBZ_BITS - 1 ||
        max_attempts == 0 || ibz_cmp(&lideal->norm, &ibz_const_zero) <= 0 ||
        ibz_is_zero(&lideal->lattice.denom))
        return QUAT_QLAPOTI_FATAL;

    uint32_t attempt_limit = max_attempts;
    if (attempt_limit > SQISIGN_QLAPOTI_GENERATOR_MAX_ATTEMPTS)
        attempt_limit = SQISIGN_QLAPOTI_GENERATOR_MAX_ATTEMPTS;

    ibz_init(&n);
    ibz_init(&d);
    ibz_init(&n2);
    ibz_init(&gcd);
    ibz_vec_4_init(&coeffs);
    quat_alg_elem_init(&candidate);
    ibz_copy(&candidate.denom, &lideal->lattice.denom);
    if (!ibz_is_one(&candidate.denom) &&
        ibz_cmp(&candidate.denom, &ibz_const_two) != 0) {
        status = QUAT_QLAPOTI_FATAL;
        goto cleanup;
    }
    ibz_mul(&n2, &lideal->norm, &lideal->norm);
    for (uint32_t attempt = 0; attempt < attempt_limit && !found; ++attempt) {
        for (int i = 0; i < 4; i++) {
            if (!ibz_rand_interval_bits(&(coeffs[i]), (uint32_t)sampling_bound_bits)) {
                status = QUAT_QLAPOTI_FATAL;
                goto cleanup;
            }
        }
        ibz_mat_4x4_eval(&(candidate.coord), &(lideal->lattice.basis), &coeffs);

        // check a_alpha invertible
        if (ibz_is_one(&(candidate.denom))) {
            ibz_mul(&gcd, &(candidate.coord[0]), &ibz_const_two);
            ibz_gcd(&gcd, &lideal->norm, &gcd);
            found = ibz_is_one(&gcd);
        } else if (0 == ibz_cmp(&ibz_const_two, &(candidate.denom))) {
            ibz_gcd(&gcd, &lideal->norm, &(candidate.coord[0]));
            found = ibz_is_one(&gcd);
        }
        // check generator
        if (found) {
            if (!quat_alg_norm(&n, &d, &candidate, alg) || !ibz_is_one(&d)) {
                status = QUAT_QLAPOTI_FATAL;
                goto cleanup;
            }
            ibz_gcd(&gcd, &n2, &n);
            found = (ibz_cmp(&gcd, &lideal->norm) == 0);
        }
    }

    if (found) {
        quat_alg_elem_copy(gen, &candidate);
        status = QUAT_QLAPOTI_SUCCESS;
    }

cleanup:
    quat_alg_elem_finalize(&candidate);
    ibz_vec_4_finalize(&coeffs);
    ibz_finalize(&n);
    ibz_finalize(&d);
    ibz_finalize(&n2);
    ibz_finalize(&gcd);
    return status;
}

// enumerate in dim 2

// helper for cvp
void
ibz_rounded_div(ibz_t *q, const ibz_t *a, const ibz_t *b)
{
    ibz_t r, abs_b;
    ibz_init(&r);
    ibz_init(&abs_b);

    // assumed to round towards 0
    const int quotient_negative =
        ibz_is_negative(a) != ibz_is_negative(b);
    ibz_abs(&abs_b, b);
    ibz_div(q, &r, a, b);
    ibz_abs(&r, &r);
    ibz_add(&r, &r, &r);
    if (ibz_cmp(&r, &abs_b) > 0) {
        /* The correction has sign(a/b).  Derive it from operand signs
         * instead of multiplying two full-width integers only to inspect the
         * sign of the product. */
        if (quotient_negative)
            ibz_sub(q, q, &ibz_const_one);
        else
            ibz_add(q, q, &ibz_const_one);
    }
    ibz_finalize(&r);
    ibz_finalize(&abs_b);
}

int
quat_dim2_lattice_contains(const ibz_mat_2x2_t *basis, const ibz_t *coord1, const ibz_t *coord2)
{
    int res = 1;
    ibz_t prod, sum, det, r;
    ibz_init(&det);
    ibz_init(&r);
    ibz_init(&sum);
    ibz_init(&prod);
    // compute det, then both coordinates (inverse*det)*vec, where vec is (coord1, coord2) and check wthether det
    // divides both results
    ibz_mat_2x2_det_from_ibz(&det, &((*basis)[0][0]), &((*basis)[0][1]), &((*basis)[1][0]), &((*basis)[1][1]));
    ibz_mul(&sum, coord1, &((*basis)[1][1]));
    ibz_mul(&prod, coord2, &((*basis)[0][1]));
    ibz_sub(&sum, &sum, &prod);
    ibz_div(&prod, &r, &sum, &det);
    res = res & ibz_is_zero(&r);
    ibz_mul(&sum, coord2, &((*basis)[0][0]));
    ibz_mul(&prod, coord1, &((*basis)[1][0]));
    ibz_sub(&sum, &sum, &prod);
    ibz_div(&prod, &r, &sum, &det);
    res = res & ibz_is_zero(&r);
    ibz_finalize(&det);
    ibz_finalize(&r);
    ibz_finalize(&sum);
    ibz_finalize(&prod);
    return (res);
}

void
quat_dim2_lattice_norm(ibz_t *norm, const ibz_t *coord1, const ibz_t *coord2, const ibz_t *norm_q)
{
    ibz_t prod;
    ibz_init(&prod);
    ibz_mul(&prod, coord2, coord2);
    if (!ibz_is_one(norm_q))
        ibz_mul(&prod, &prod, norm_q);
    /* Every potentially aliasing input has now been consumed. */
    ibz_mul(norm, coord1, coord1);
    ibz_add(norm, norm, &prod);
    ibz_finalize(&prod);
}

void
quat_dim2_lattice_bilinear(ibz_t *res,
                           const ibz_t *v11,
                           const ibz_t *v12,
                           const ibz_t *v21,
                           const ibz_t *v22,
                           const ibz_t *norm_q)
{
    ibz_t prod;
    ibz_init(&prod);
    ibz_mul(&prod, v12, v22);
    if (!ibz_is_one(norm_q))
        ibz_mul(&prod, &prod, norm_q);
    /* Every potentially aliasing second-term input is consumed first. */
    ibz_mul(res, v11, v21);
    ibz_add(res, res, &prod);
    ibz_finalize(&prod);
}

// algo 3.1.14 Cohen (exact solution for shortest vector in dimension 2, than take a second, orthogonal vector)
void
quat_dim2_lattice_short_basis(ibz_mat_2x2_t *reduced, const ibz_mat_2x2_t *basis, const ibz_t *norm_q)
{
    ibz_vec_2_t a, b;
    ibz_t prod, norm_a, norm_b, r, norm_t, n;
    ibz_vec_2_init(&a);
    ibz_vec_2_init(&b);
    ibz_init(&prod);
    ibz_init(&r);
    ibz_init(&n);
    ibz_init(&norm_t);
    ibz_init(&norm_a);
    ibz_init(&norm_b);
    // init a,b
    ibz_copy(&(a[0]), &((*basis)[0][0]));
    ibz_copy(&(a[1]), &((*basis)[1][0]));
    ibz_copy(&(b[0]), &((*basis)[0][1]));
    ibz_copy(&(b[1]), &((*basis)[1][1]));
    // compute initial norms
    quat_dim2_lattice_norm(&norm_a, &(a[0]), &(a[1]), norm_q);
    quat_dim2_lattice_norm(&norm_b, &(b[0]), &(b[1]), norm_q);
    // exchange if needed
    if (ibz_cmp(&norm_a, &norm_b) < 0) {
        ibz_swap(&a[0], &b[0]);
        ibz_swap(&a[1], &b[1]);
        ibz_swap(&norm_a, &norm_b);
    }
    int test = 1;
    while (test) {
        // compute n
        quat_dim2_lattice_bilinear(&n, &(a[0]), &(a[1]), &(b[0]), &(b[1]), norm_q);
        // set r
        ibz_rounded_div(&r, &n, &norm_b);
        // compute t_norm
        ibz_add(&prod, &n, &n);
        ibz_mul(&prod, &prod, &r);
        ibz_sub(&norm_t, &norm_a, &prod);
        ibz_mul(&prod, &r, &r);
        ibz_mul(&prod, &prod, &norm_b);
        ibz_add(&norm_t, &norm_t, &prod);
        // test:
        if (ibz_cmp(&norm_b, &norm_t) > 0) {
            // compute t, a, b
            ibz_copy(&norm_a, &norm_b);
            ibz_copy(&norm_b, &norm_t);
            /* norm_t and n are dead until the next iteration; use them for
             * t = a-rb while rotating (a,b) without a third vector. */
            ibz_mul(&prod, &r, &(b[0]));
            ibz_sub(&norm_t, &(a[0]), &prod);
            ibz_mul(&prod, &r, &(b[1]));
            ibz_sub(&n, &(a[1]), &prod);
            ibz_copy(&(a[0]), &(b[0]));
            ibz_copy(&(a[1]), &(b[1]));
            ibz_copy(&(b[0]), &norm_t);
            ibz_copy(&(b[1]), &n);
        } else {
            test = 0;
        }
    }
    // output : now b is short: need to get 2nd short vector: idea: take shortest among t and a
    if (ibz_cmp(&norm_t, &norm_a) < 0) {
        ibz_mul(&prod, &r, &(b[0]));
        ibz_sub(&(a[0]), &(a[0]), &prod);
        ibz_mul(&prod, &r, &(b[1]));
        ibz_sub(&(a[1]), &(a[1]), &prod);
    }
    ibz_copy(&((*reduced)[0][0]), &(b[0]));
    ibz_copy(&((*reduced)[1][0]), &(b[1]));
    ibz_copy(&((*reduced)[0][1]), &(a[0]));
    ibz_copy(&((*reduced)[1][1]), &(a[1]));

    ibz_finalize(&prod);
    ibz_finalize(&norm_a);
    ibz_finalize(&norm_b);
    ibz_finalize(&norm_t);
    ibz_vec_2_finalize(&a);
    ibz_vec_2_finalize(&b);
    ibz_finalize(&r);
    ibz_finalize(&n);
}

// qlapoti
int
quat_qlapoti_check_mod_condition(const ibz_t *m, const ibz_t *a, const ibz_t *b)
{
    int ok = 1;
    int A = ibz_get(a) & 1;
    int B = ibz_get(b) & 1;
    int M8 = ibz_get(m) & 7;
    int M = M8 & 3;
#ifndef NDEBUG
    if (DEBUG_PRINTS > 3)
        ibz_printf("m %d a %d b %d\n", M, A, B);
#endif
    if (M8 == 0)
        ok = 0;
    if (ok) {
        if ((A == B) && (A == 0)) {
            if (M != 0) {
                ok = 0;
#ifndef NDEBUG
                if (DEBUG_PRINTS > 2)
                    ibz_printf("case 0\n");
#endif
            }
        } else {
            if (((A == B) && (A == 1))) {
                if (M != 2) {
                    ok = 0;
#ifndef NDEBUG
                    if (DEBUG_PRINTS > 2)
                        ibz_printf("case 1\n");
#endif
                }
            } else {
                if (M != 1) {
#ifndef NDEBUG
                    if (DEBUG_PRINTS > 2)
                        ibz_printf("mod 4\n");
#endif
                    ok = 0;
                }
            }
        }
    }

    return (ok);
}

int
quat_dim2_lattice_qlapoti_cvp_condition(quat_alg_elem_t *elem, const ibz_vec_2_t *vec, const void *params)
{
    // remember a_alpha and b_alpha are integers
    int found = 1;
    const qlapoti_enumeration_parameters_t *q_params =
        (const qlapoti_enumeration_parameters_t *)params;
    if (elem == NULL || vec == NULL || q_params == NULL ||
        q_params->m == NULL || q_params->n == NULL ||
        q_params->a_alpha == NULL || q_params->b_alpha == NULL ||
        q_params->cornacchia_params == NULL ||
        ibz_cmp(q_params->n, &ibz_const_zero) <= 0 ||
        q_params->cornacchia_params->q <= 0 ||
        q_params->cornacchia_params->primality_test_iterations <= 0 ||
        q_params->cornacchia_params->primality_test_iterations >
            SQISIGN_QLAPOTI_MAX_PRIMALITY_ITERATIONS ||
        q_params->cornacchia_params->prime_list_length >
            QLAPOTI_MAX_CORNACCHIA_PRIMES ||
        (q_params->cornacchia_params->prime_list_length != 0 &&
         q_params->cornacchia_params->prime_list == NULL))
        return QUAT_QLAPOTI_FATAL;
    ibz_t m2, tmp, A, B;
    ibz_vec_2_t a, b;
    ibz_init(&m2);
    ibz_init(&tmp);
    ibz_init(&A);
    ibz_init(&B);
    ibz_vec_2_init(&a);
    ibz_vec_2_init(&b);

    // compute A,B
    // Assumes vec is target-close (up to one sign)
    ibz_neg(&A, &((*vec)[0]));
    ibz_neg(&B, &((*vec)[1]));
#ifndef NDEBUG
    if (DEBUG_PRINTS > 3)
        ibz_printf("A %Zd\nB %Zd\n", &A, &B);
    // assert (2*(a_alpha*A + b_alpha*B)) % N == M % N
    ibz_mul(&tmp, &A, q_params->a_alpha);
    ibz_mul(&m2, &B, q_params->b_alpha);
    ibz_add(&tmp, &tmp, &m2);
    ibz_mod(&tmp, &tmp, q_params->n);
    ibz_mod(&m2, q_params->m, q_params->n);
    assert(0 == ibz_cmp(&m2, &tmp));
#endif
    // compute m2
    // M2 = M - 2 *a_alpha *A - 2 *b_alpha *B; M2 = ZZ(M2 / N)
    ibz_mul(&tmp, &A, q_params->a_alpha);
    ibz_mul(&m2, &B, q_params->b_alpha);
    ibz_add(&tmp, &tmp, &m2);
    ibz_sub(&m2, q_params->m, &tmp);
    ibz_div(&m2, &tmp, &m2, q_params->n);
    if (!ibz_is_zero(&tmp)) {
        found = QUAT_QLAPOTI_FATAL;
        goto cleanup;
    }
    // Complete the square
    // Comp: M3 = M2 - A ** 2 - B ** 2
    ibz_mul(&tmp, &A, &A);
    ibz_sub(&m2, &m2, &tmp);
    ibz_mul(&tmp, &B, &B);
    ibz_sub(&m2, &m2, &tmp);

    // Comp; M4 = 2 * M3 + A * *2 + B * *2
    ibz_add(&m2, &m2, &m2);
    ibz_add(&m2, &m2, &tmp);
    ibz_mul(&tmp, &A, &A);
    ibz_add(&m2, &m2, &tmp);

#ifndef NDEBUG
    if ((DEBUG_PRINTS > 2))
        ibz_printf("M4 %Zd\n", &m2);
    // If the first one is too small, there is no point in trying others...
    // if M4 < 0:if first_vec : break else : continue
    // Must be communicated to exterior loop (enforce change of alpha)? Might be fine already since better enum (if
    // correct bound)
    // Test for unsolvable cases
    if ((DEBUG_PRINTS > 2) && found && !(ibz_cmp(&ibz_const_zero, &m2) < 0))
        ibz_printf("size\n");
#endif
    if (ibz_cmp(&m2, &ibz_const_zero) <= 0)
        found = 0;
    if (found > 0)
        found = quat_qlapoti_check_mod_condition(&m2, &A, &B);
    if (found > 0) {
        // cornacchia
        found = ibz_cornacchia_extended(&(a[0]), &(b[0]), &m2, q_params->cornacchia_params);
#ifndef NDEBUG
        if ((DEBUG_PRINTS > 2) && !found)
            ibz_printf("cor \n");
#endif
    }
    if (found) {
// treat output
#ifndef NDEBUG
        if (DEBUG_PRINTS > 3)
            ibz_printf("a_0 %Zd, b_0 %Zd\n", &(a[0]), &(b[0]));
#endif
        if (!((ibz_get(&A) & 1) == (ibz_get(&(a[0])) & 1)))
            ibz_swap(&(a[0]), &(b[0]));
        if (((ibz_get(&A) & 1) != (ibz_get(&(a[0])) & 1)) ||
            ((ibz_get(&B) & 1) != (ibz_get(&(b[0])) & 1))) {
            found = QUAT_QLAPOTI_FATAL;
            goto cleanup;
        }
        // a1 = ZZ((ad1 + A)/2);b1 = ZZ((bd1 + B)/2)
        ibz_add(&(a[0]), &(a[0]), &A);
        ibz_div(&(a[0]), &tmp, &(a[0]), &ibz_const_two);
        if (!ibz_is_zero(&tmp)) {
            found = QUAT_QLAPOTI_FATAL;
            goto cleanup;
        }
        ibz_add(&(b[0]), &(b[0]), &B);
        ibz_div(&(b[0]), &tmp, &(b[0]), &ibz_const_two);
        if (!ibz_is_zero(&tmp)) {
            found = QUAT_QLAPOTI_FATAL;
            goto cleanup;
        }
        // a2 = A - a1;b2 = B - b1
        ibz_sub(&(a[1]), &A, &(a[0]));
        ibz_sub(&(b[1]), &B, &(b[0]));
        ibz_copy(&(elem->coord[0]), &(a[0]));
        ibz_copy(&(elem->coord[1]), &(a[1]));
        ibz_copy(&(elem->coord[2]), &(b[0]));
        ibz_copy(&(elem->coord[3]), &(b[1]));
    }
cleanup:
    ibz_finalize(&m2);
    ibz_finalize(&tmp);
    ibz_finalize(&A);
    ibz_finalize(&B);
    ibz_vec_2_finalize(&a);
    ibz_vec_2_finalize(&b);
    return (found);
}

int
quat_elem_is_odd_norm(const quat_alg_elem_t *elem)
{
    int found = 0;
    if ((ibz_get(&elem->coord[0]) & 1) == (ibz_get(&elem->coord[2]) & 1)) {
        return found;
    }
    for (int i = 0; i < 4; i++) {
        if (found) {
            if ((ibz_get(&elem->coord[i]) & 3) == 2) {
                return 0;
            }
        } else {
            if ((ibz_get(&elem->coord[i]) & 3) == 2) {
                found = 1;
            }
        }
    }
    return found;
}

int
get_endtype(const quat_alg_elem_t *elem)
{
    int t1, t2, t3, t4;
    // transformlist = {1: [[2, 0, 1, 2], [2, 2, 3, 0], [2, 2, 1, 0], [2, 0, 3, 2]], 2: [[0, 2, 2, 1], [0, 2, 2, 3], [2,
    // 2, 0, 1], [2, 2, 0, 3]]}
    t1 = ibz_get(&elem->coord[0]) & 3;
    t2 = ibz_get(&elem->coord[1]) & 3;
    // NB! These are swapped because I computed this using j and k swapped
    t3 = ibz_get(&elem->coord[3]) & 3;
    t4 = ibz_get(&elem->coord[2]) & 3;
    // End NB!
    if (t1 == 2) {
        if (t2 == 2) {
            if (t3 == 1 && t4 == 0) {
                return 1;
            } else if (t3 == 3 && t4 == 2) {
                return 1;
            } else if (t3 == 0) {
                if (t4 == 1) {
                    return 2;
                } else if (t4 == 3) {
                    return 2;
                } else {
                    return 0;
                }
            } else {
                return 0;
            }
        } else if (t2 == 0) {
            if (t3 == 1 && t4 == 2) {
                return 1;
            } else if (t3 == 3 && t4 == 2) {
                return 1;
            } else {
                return 0;
            }
        } else {
            return 0;
        }
    } else if (t1 == 0) {
        if (t2 == 2 && t3 == 2) {
            if (t4 == 1) {
                return 2;
            } else if (t4 == 3) {
                return 2;
            } else {
                return 0;
            }
        } else {
            return 0;
        }
    } else {
        return 0;
    }
}
/* Reduce the equivalent ideal in a separate phase so its 4x4 Gram matrix is
 * not retained by the bounded Qlapoti search. */
static QLAPOTI_NOINLINE int
quat_qlapoti_reduce_small_basis(quat_left_ideal_t *small,
                                const quat_alg_t *alg)
{
    ibz_mat_4x4_t gram;
    int status;

    ibz_mat_4x4_init(&gram);
    status = quat_lideal_reduce_basis(
        &small->lattice.basis, &gram, small, alg);
    ibz_mat_4x4_finalize(&gram);
    return status;
}

/* Build and reduce the two-dimensional congruence lattice.  Only its
 * adjugate (called L_inv in the original implementation) survives this
 * phase.  quat_dim2_lattice_short_basis first copies both input columns, so
 * using L as both input and output is safe. */
static QLAPOTI_NOINLINE int
quat_qlapoti_prepare_2d_lattice(ibz_mat_2x2_t *L_inv,
                                int *keep_alpha,
                                const ibz_t *n,
                                const ibz_t *b_alpha,
                                const ibz_t *a_alpha,
                                const ibz_t *a_alpha_inv,
                                const ibz_t *psqrt)
{
    ibz_mat_2x2_t L;
    int status = 0;

    ibz_mat_2x2_init(&L);

    /* x = b_alpha/a_alpha (mod n), followed by
     * L = [[n-x, n], [1, 0]]. */
    ibz_mul(&L[0][0], b_alpha, a_alpha_inv);
    ibz_mod(&L[0][0], &L[0][0], n);
    ibz_sub(&L[0][0], n, &L[0][0]);
    ibz_copy(&L[0][1], n);
    ibz_set(&L[1][0], 1);
    ibz_set(&L[1][1], 0);
    quat_dim2_lattice_short_basis(&L, &L, &ibz_const_one);

    if (!ibz_mat_2x2_inv_with_det_as_denom(L_inv, NULL, &L))
        goto cleanup;

    /* L is dead after its adjugate has been published.  Reuse its entries
     * for the suitability test instead of adding scalar temporaries. */
    ibz_gcd(&L[0][0], a_alpha, b_alpha);
    *keep_alpha = ibz_is_one(&L[0][0]);
    ibz_mul(&L[0][1], &L[0][1], &L[0][1]);
    ibz_mul(&L[1][1], &L[1][1], &L[1][1]);
    ibz_add(&L[0][1], &L[0][1], &L[1][1]);
    *keep_alpha = *keep_alpha && (ibz_cmp(&L[0][1], psqrt) < 0);
    status = 1;

cleanup:
    ibz_mat_2x2_finalize(&L);
    return status;
}

/* L_inv is adj(L_red).  Evaluate L_red * vec directly from the adjugate so
 * the search does not need to retain both 2x2 matrices.  The private result
 * also makes res == vec safe. */
static QLAPOTI_NOINLINE void
quat_qlapoti_red_eval_from_inverse(ibz_vec_2_t *res,
                                    const ibz_mat_2x2_t *L_inv,
                                    const ibz_vec_2_t *vec)
{
    ibz_vec_2_t product;
    ibz_t tmp;

    ibz_vec_2_init(&product);
    ibz_init(&tmp);

    ibz_mul(&product[0], &(*L_inv)[1][1], &(*vec)[0]);
    ibz_mul(&tmp, &(*L_inv)[0][1], &(*vec)[1]);
    ibz_sub(&product[0], &product[0], &tmp);
    ibz_mul(&product[1], &(*L_inv)[0][0], &(*vec)[1]);
    ibz_mul(&tmp, &(*L_inv)[1][0], &(*vec)[0]);
    ibz_sub(&product[1], &product[1], &tmp);

    ibz_copy(&(*res)[0], &product[0]);
    ibz_copy(&(*res)[1], &product[1]);
    ibz_finalize(&tmp);
    ibz_vec_2_finalize(&product);
}

// document all functions
static QLAPOTI_NOINLINE int
quat_qlapoti_search(quat_alg_elem_t *mu1,
                    quat_alg_elem_t *mu2,
                    quat_alg_elem_t *theta,
                    const quat_left_ideal_t *small,
                    const quat_alg_t *alg,
                    int max_counter_alpha,
                    int gen_sampling_bound_bits,
                    int two_power,
                    const ibz_cornacchia_extended_params_t *cornacchia_params)
{
    int found = 0;
    int transformed = 0;
    int keep_alpha = 0;
    int result = QUAT_QLAPOTI_RETRY;
    quat_alg_elem_t alpha, encoding, alpha_0;
    ibz_mat_2x2_t L_inv;
    ibz_vec_2_t v_target, v_close;
    ibz_t n, m, a_alpha, b_alpha, alpha0norm, psqrt;
    ibz_t tmp, lam, two_e;
    qlapoti_enumeration_parameters_t params;
    quat_alg_elem_t gamma1, gamma2;
    quat_alg_elem_init(&gamma1);
    quat_alg_elem_init(&gamma2);
    ibz_vec_2_init(&v_target);
    ibz_vec_2_init(&v_close);
    ibz_mat_2x2_init(&L_inv);
    ibz_init(&tmp);
    ibz_init(&n);
    ibz_init(&m);
    ibz_init(&two_e);
    ibz_init(&lam);
    ibz_init(&a_alpha);
    ibz_init(&b_alpha);
    ibz_init(&alpha0norm);
    ibz_init(&psqrt);
    quat_alg_elem_init(&alpha);
    quat_alg_elem_init(&alpha_0);
    quat_alg_elem_init(&encoding);
    ibz_sqrt_floor(&psqrt, &alg->p);
    params.cornacchia_params = cornacchia_params;
    ibz_copy(&n, &small->norm);
    if (ibz_cmp(&n, &ibz_const_zero) <= 0) {
        result = QUAT_QLAPOTI_FATAL;
        goto cleanup;
    }
    ibz_pow(&two_e, &ibz_const_two, two_power);

    // int num_alphas_tried = 0;
    // int num_loops_total = 0;
    const uint32_t counter_alpha = (uint32_t)max_counter_alpha;
    const uint32_t total_attempts =
        counter_alpha * SQISIGN_QLAPOTI_COUNTER_MULTIPLIER;
    for (uint32_t counter = 0; counter < total_attempts; counter++) {
        /* Success is established only after CVP and all output invariants.
         * Keeping a stale true value across a bounded-search `continue` can
         * otherwise turn a final generator/lambda/inverse miss into success. */
        found = 0;

        // compute alpha, lam, alpha_0, alpha0norm
        {
            if (counter % counter_alpha == 0) {
                keep_alpha = 0;
                ibz_set(&lam, 1);
            }
            if (!keep_alpha) {
                // num_alphas_tried++;
                int generator_status = quat_lideal_generator_small_coprime(
                    &alpha_0,
                    small,
                    alg,
                    gen_sampling_bound_bits,
                    counter_alpha);
                if (generator_status == QUAT_QLAPOTI_FATAL) {
                    result = QUAT_QLAPOTI_FATAL;
                    goto cleanup;
                }
                if (generator_status == QUAT_QLAPOTI_RETRY) {
                    keep_alpha = 0;
                    continue;
                }
                quat_alg_elem_copy(&alpha, &alpha_0);
                if (!quat_alg_norm(&alpha0norm, &b_alpha, &alpha_0, alg) ||
                    !ibz_is_one(&b_alpha)) {
                    result = QUAT_QLAPOTI_FATAL;
                    goto cleanup;
                }
                ibz_set(&lam, 1);
            } else {
                uint32_t lambda_attempts = 0;
                do {
                    ibz_add(&lam, &lam, &ibz_const_one);
                    /* alpha and alpha_0 retain the same positive denominator
                     * (generator invariant: 1 or 2).  Accumulate their raw
                     * numerators and defer normalization until CVP succeeds;
                     * generic rational addition/normalization here dominated
                     * every lambda retry while producing the same element. */
                    ibz_vec_4_add(&alpha.coord,
                                  &alpha.coord,
                                  &alpha_0.coord);
                    ibz_gcd(&b_alpha, &small->norm, &lam);
                    lambda_attempts++;
                } while (!ibz_is_one(&b_alpha) &&
                         lambda_attempts <
                             SQISIGN_QLAPOTI_LAMBDA_MAX_ATTEMPTS);
                if (!ibz_is_one(&b_alpha)) {
                    keep_alpha = 0;
                    continue;
                }
            }
            if (!quat_lattice_contains(NULL, &small->lattice, &alpha)) {
                result = QUAT_QLAPOTI_FATAL;
                goto cleanup;
            }
        }

        // compute m
        // 2^e-(2norm(alpha0)lam^2)/n
        {
            ibz_add(&m, &alpha0norm, &alpha0norm);
            ibz_mul(&tmp, &lam, &lam);
            ibz_mul(&m, &m, &tmp);
            ibz_div_floor(&m, &tmp, &m, &n);
            if (!ibz_is_zero(&tmp)) {
                result = QUAT_QLAPOTI_FATAL;
                goto cleanup;
            }
            ibz_sub(&m, &two_e, &m);
            // Discard bad alphas directly: this case should never happen
            if (ibz_cmp(&m, &ibz_const_zero) <= 0) {
                found = 0;
                continue;
            }
        }
        // set a_alpha, b_alpha (to double of what they are in sage)
        {
            ibz_copy(&a_alpha, &(alpha.coord[0]));
            ibz_copy(&b_alpha, &(alpha.coord[1]));
            if (ibz_is_one(&alpha.denom)) {
                ibz_add(&a_alpha, &a_alpha, &a_alpha);
                ibz_add(&b_alpha, &b_alpha, &b_alpha);
            }
            else if (ibz_cmp(&alpha.denom, &ibz_const_two) != 0) {
                result = QUAT_QLAPOTI_FATAL;
                goto cleanup;
            }
        }
        // Prepare target vector
        {
            // #A + (2*b_alpha/2*a_alpha)B = M/(2*a_alpha) (mod N)
            // x = ZZ(Z_N(2*b_alpha)*(Z_N(2*a_alpha)**-1))
            if (!ibz_invmod(&tmp, &a_alpha, &n)) {
                keep_alpha = 0;
                continue;
            }
            // T = ZZ(Z_N(M)*Z_N(2*a_alpha)**-1)
            ibz_mul(&v_target[0], &tmp, &m);
            ibz_mod(&v_target[0], &v_target[0], &n);
            // v = vector(ZZ, [-T, 0])
            ibz_set(&(v_target[1]), 0);
            ibz_neg(&(v_target[0]), &(v_target[0]));
        }
        // Prepare the reduced 2d lattice adjugate. Set keep_alpha accordingly.
        // tmp must contain invmod(a_alpha,n) here
        if (!keep_alpha) {
            if (!quat_qlapoti_prepare_2d_lattice(
                    &L_inv,
                    &keep_alpha,
                    &n,
                    &b_alpha,
                    &a_alpha,
                    &tmp,
                    &psqrt)) {
                result = QUAT_QLAPOTI_FATAL;
                goto cleanup;
            }
        }
        // find close vector
        {
            /* v_target[1] is set to zero above and is unchanged here.  Avoid
             * the two zero products and the generic matrix-evaluation frame. */
            ibz_mul(&v_close[0], &L_inv[0][0], &v_target[0]);
            ibz_mul(&v_close[1], &L_inv[1][0], &v_target[0]);
            /* det(adj(L_red)) = det(L_red) in dimension two.  tmp is dead
             * after lattice preparation and can hold this iteration's det. */
            ibz_mat_2x2_det_from_ibz(&tmp,
                                     &L_inv[0][0],
                                     &L_inv[0][1],
                                     &L_inv[1][0],
                                     &L_inv[1][1]);
            if (ibz_is_zero(&tmp)) {
                result = QUAT_QLAPOTI_FATAL;
                goto cleanup;
            }
            ibz_rounded_div(&(v_close[0]), &(v_close[0]), &tmp);
            ibz_rounded_div(&(v_close[1]), &(v_close[1]), &tmp);
            quat_qlapoti_red_eval_from_inverse(&v_close, &L_inv, &v_close);
            ibz_sub(&(v_target[0]), &(v_target[0]), &(v_close[0]));
            ibz_sub(&(v_target[1]), &(v_target[1]), &(v_close[1]));
        }
        // set parameters and call condition
        params.a_alpha = &a_alpha;
        params.b_alpha = &b_alpha;
        params.m = &m;
        params.n = &n;
        int cvp_status =
            quat_dim2_lattice_qlapoti_cvp_condition(&encoding, &v_target, &params);
        if (cvp_status == QUAT_QLAPOTI_FATAL) {
            result = QUAT_QLAPOTI_FATAL;
            goto cleanup;
        }
        found = cvp_status;

        if (!found) {
            continue;
        } else {
            /* Restore the canonical representation once on the successful
             * path.  This preserves the former downstream gamma encodings
             * while avoiding normalization on every rejected attempt. */
            quat_alg_normalize(&alpha);
            // extract ideals from encoding of output
            ibz_mul(&(gamma1.coord[0]), &(encoding.coord[0]), &n);
            ibz_mul(&(gamma2.coord[0]), &(encoding.coord[1]), &n);
            ibz_mul(&(gamma1.coord[1]), &(encoding.coord[2]), &n);
            ibz_mul(&(gamma2.coord[1]), &(encoding.coord[3]), &n);
            quat_alg_add(&gamma1, &gamma1, &alpha);
            quat_alg_add(&gamma2, &gamma2, &alpha);
#ifndef NDEBUG
            ibz_t n1, n2, d;
            ibz_init(&n1);
            ibz_init(&n2);
            ibz_init(&d);
            assert(quat_lattice_contains(NULL, &small->lattice, &gamma1));
            assert(quat_lattice_contains(NULL, &small->lattice, &gamma2));
            assert(quat_alg_norm(&n1, &d, &gamma1, alg));
            assert(ibz_is_one(&d));
            assert(quat_alg_norm(&n2, &d, &gamma2, alg));
            assert(ibz_is_one(&d));
            ibz_div(&n1, &d, &n1, &small->norm);
            assert(ibz_is_zero(&d));
            ibz_div(&n2, &d, &n2, &small->norm);
            assert(ibz_is_zero(&d));
            ibz_add(&n2, &n1, &n2);
            assert(0 == ibz_cmp(&two_e, &n2));
            ibz_finalize(&n1);
            ibz_finalize(&n2);
            ibz_finalize(&d);
#endif
            /* `encoding` is no longer needed after extraction, so reuse it as
             * the transactional theta candidate.  alpha_0 is also dead once
             * CVP succeeds; reuse it as conjugate scratch instead of keeping
             * an additional search element live.  Every rejection below resets
             * keep_alpha, so the next iteration resamples alpha_0. */
            quat_alg_conj(&alpha_0, &gamma1);
            ibz_mul(&alpha_0.denom, &alpha_0.denom, &n);
            quat_alg_mul(&encoding, &gamma2, &alpha_0, alg);
#ifndef NDEBUG
            quat_left_ideal_t I1, I2;
            quat_left_ideal_init(&I1);
            quat_left_ideal_init(&I2);
            int ideals_ok = quat_lideal_mul(&I1, small, &alpha_0, alg);
            quat_alg_elem_copy(&alpha_0, &gamma2);
            ibz_neg(&(alpha_0.coord[0]), &(alpha_0.coord[0]));
            ibz_mul(&alpha_0.denom, &alpha_0.denom, &n);
            ideals_ok = ideals_ok && quat_lideal_mul(&I2, small, &alpha_0, alg);
            if (ideals_ok) {
                ibz_add(&tmp, &I1.norm, &I2.norm);
                assert(ibz_cmp(&tmp, &two_e) == 0);
            }
            assert(ideals_ok);
            quat_left_ideal_finalize(&I1);
            quat_left_ideal_finalize(&I2);
#endif

            // enforces suitable output
            quat_alg_normalize(&encoding);
            if (ibz_is_one(&encoding.denom)) {
                int endtype = get_endtype(&encoding);
                if (endtype == 1) {
                    quat_alg_add(&alpha_0, &gamma1, &gamma2);
                    ibz_add(&alpha_0.denom, &alpha_0.denom, &alpha_0.denom);
                    quat_alg_sub(&gamma1, &gamma1, &gamma2);
                    ibz_add(&gamma1.denom, &gamma1.denom, &gamma1.denom);
                    quat_alg_elem_copy(&gamma2, &alpha_0);
                } else if (endtype == 2) {
                    quat_alg_elem_set(&alpha_0, 1, 0, 1, 0, 0);
                    quat_alg_mul(&gamma2, &alpha_0, &gamma2, alg);
                    quat_alg_add(&alpha_0, &gamma1, &gamma2);
                    ibz_add(&alpha_0.denom, &alpha_0.denom, &alpha_0.denom);
                    quat_alg_sub(&gamma1, &gamma1, &gamma2);
                    ibz_add(&gamma1.denom, &gamma1.denom, &gamma1.denom);
                    quat_alg_elem_copy(&gamma2, &alpha_0);
                } else {
                    keep_alpha = 0;
                    ibz_set(&lam, 1);
                    ibz_vec_4_set(&gamma1.coord, 0, 0, 0, 0);
                    ibz_vec_4_set(&gamma2.coord, 0, 0, 0, 0);
                    ibz_set(&gamma1.denom, 1);
                    ibz_set(&gamma2.denom, 1);
                    found = 0;
                    transformed = 0;
                    continue;
                }

                if (!quat_lattice_contains(NULL, &small->lattice, &gamma1) ||
                    !quat_lattice_contains(NULL, &small->lattice, &gamma2)) {
                    result = QUAT_QLAPOTI_FATAL;
                    goto cleanup;
                }
                quat_alg_conj(&alpha_0, &gamma1);
                ibz_mul(&alpha_0.denom, &alpha_0.denom, &n);
                quat_alg_mul(&encoding, &gamma2, &alpha_0, alg);
                quat_alg_normalize(&encoding);
                transformed = 1;
            } else {
                if (!quat_elem_is_odd_norm(&encoding)) {
                    keep_alpha = 0;
                    ibz_set(&lam, 1);
                    ibz_vec_4_set(&gamma1.coord, 0, 0, 0, 0);
                    ibz_vec_4_set(&gamma2.coord, 0, 0, 0, 0);
                    ibz_set(&gamma1.denom, 1);
                    ibz_set(&gamma2.denom, 1);
                    found = 0;
                    transformed = 0;
                    continue;
                }
            }
            found = 1;
            // num_loops_total = counter;
            break;
        }
    }

    if (found > 0) {
        quat_alg_elem_copy(mu1, &gamma1);
        quat_alg_elem_copy(mu2, &gamma2);
        quat_alg_elem_copy(theta, &encoding);
        result = transformed ? QUAT_QLAPOTI_SUCCESS_TRANSFORMED
                             : QUAT_QLAPOTI_SUCCESS;
    }

cleanup:
    quat_alg_elem_finalize(&gamma1);
    quat_alg_elem_finalize(&gamma2);
    ibz_finalize(&n);
    ibz_finalize(&m);
    ibz_finalize(&tmp);
    ibz_finalize(&lam);
    ibz_finalize(&two_e);
    ibz_finalize(&a_alpha);
    ibz_finalize(&b_alpha);
    ibz_finalize(&alpha0norm);
    ibz_finalize(&psqrt);
    quat_alg_elem_finalize(&alpha);
    quat_alg_elem_finalize(&alpha_0);
    quat_alg_elem_finalize(&encoding);
    ibz_vec_2_finalize(&v_target);
    ibz_vec_2_finalize(&v_close);
    ibz_mat_2x2_finalize(&L_inv);

    return result;
}

int
quat_qlapoti(quat_alg_elem_t *mu1,
             quat_alg_elem_t *mu2,
             quat_alg_elem_t *theta,
             quat_alg_elem_t *smallest,
             const quat_left_ideal_t *lideal,
             const quat_alg_t *alg,
             int max_counter_alpha,
             int gen_sampling_bound_bits,
             int two_power,
             const ibz_cornacchia_extended_params_t *cornacchia_params)
{
    quat_left_ideal_t small;
    quat_alg_elem_t smallest_candidate;
    int result = QUAT_QLAPOTI_FATAL;

    if (mu1 == NULL || mu2 == NULL || theta == NULL || smallest == NULL ||
        lideal == NULL || alg == NULL || cornacchia_params == NULL ||
        max_counter_alpha <= 0 ||
        (uint32_t)max_counter_alpha > SQISIGN_QLAPOTI_MAX_COUNTER_ALPHA ||
        gen_sampling_bound_bits <= 0 ||
        gen_sampling_bound_bits >= IBZ_BITS - 1 ||
        two_power <= 0 || two_power >= IBZ_BITS - 1 ||
        cornacchia_params->q <= 0 ||
        cornacchia_params->primality_test_iterations <= 0 ||
        cornacchia_params->primality_test_iterations >
            SQISIGN_QLAPOTI_MAX_PRIMALITY_ITERATIONS ||
        cornacchia_params->prime_list_length > QLAPOTI_MAX_CORNACCHIA_PRIMES ||
        (cornacchia_params->prime_list_length != 0 &&
         cornacchia_params->prime_list == NULL) ||
        ibz_cmp(&lideal->norm, &ibz_const_zero) <= 0 ||
        ibz_is_zero(&lideal->lattice.denom))
        return QUAT_QLAPOTI_FATAL;
    for (unsigned i = 0; i < cornacchia_params->prime_list_length; ++i) {
        if (cornacchia_params->prime_list[i] <= 1)
            return QUAT_QLAPOTI_FATAL;
    }

    quat_left_ideal_init(&small);
    quat_alg_elem_init(&smallest_candidate);

    /* The shortest-equivalent reduction, reduced-basis preparation, and
     * probabilistic search have disjoint large scratch lifetimes. */
    if (!quat_lideal_shortest_equivalent(
            &small, &smallest_candidate, lideal, alg))
        goto cleanup;
    if (!quat_qlapoti_reduce_small_basis(&small, alg))
        goto cleanup;

    result = quat_qlapoti_search(mu1,
                                 mu2,
                                 theta,
                                 &small,
                                 alg,
                                 max_counter_alpha,
                                 gen_sampling_bound_bits,
                                 two_power,
                                 cornacchia_params);
    if (result == QUAT_QLAPOTI_SUCCESS ||
        result == QUAT_QLAPOTI_SUCCESS_TRANSFORMED)
        quat_alg_elem_copy(smallest, &smallest_candidate);

cleanup:
    quat_alg_elem_finalize(&smallest_candidate);
    quat_left_ideal_finalize(&small);
    return result;
}
