#include <quaternion.h>

#include "internal.h"
#include <stdint.h>

#if defined(__GNUC__) || defined(__clang__)
#define NORMEQ_NOINLINE __attribute__((noinline))
#else
#define NORMEQ_NOINLINE
#endif

/* These loops succeed with constant (gamma) or overwhelming (beta)
 * probability for the prime norms used by SQIsign.  The generous limits are
 * defense-in-depth against a broken entropy source or malformed inputs, and
 * may be lowered by an embedded build after measuring its parameter sets. */
#ifndef SQISIGN_RANDOM_IDEAL_GAMMA_MAX_ATTEMPTS
#define SQISIGN_RANDOM_IDEAL_GAMMA_MAX_ATTEMPTS 4096U
#endif

#ifndef SQISIGN_RANDOM_IDEAL_BETA_MAX_ATTEMPTS
#define SQISIGN_RANDOM_IDEAL_BETA_MAX_ATTEMPTS 4096U
#endif

#ifndef SQISIGN_REPRESENT_INTEGER_MAX_ATTEMPTS
#define SQISIGN_REPRESENT_INTEGER_MAX_ATTEMPTS 1048576U
#endif

#ifndef SQISIGN_NORMEQ_MAX_PRIMALITY_ITERATIONS
#define SQISIGN_NORMEQ_MAX_PRIMALITY_ITERATIONS 256
#endif

#if SQISIGN_RANDOM_IDEAL_GAMMA_MAX_ATTEMPTS == 0
#error "SQISIGN_RANDOM_IDEAL_GAMMA_MAX_ATTEMPTS must be positive"
#endif

#if SQISIGN_RANDOM_IDEAL_BETA_MAX_ATTEMPTS == 0
#error "SQISIGN_RANDOM_IDEAL_BETA_MAX_ATTEMPTS must be positive"
#endif

#if SQISIGN_REPRESENT_INTEGER_MAX_ATTEMPTS == 0
#error "SQISIGN_REPRESENT_INTEGER_MAX_ATTEMPTS must be positive"
#endif

#if SQISIGN_NORMEQ_MAX_PRIMALITY_ITERATIONS <= 0
#error "SQISIGN_NORMEQ_MAX_PRIMALITY_ITERATIONS must be positive"
#endif

_Static_assert(SQISIGN_RANDOM_IDEAL_GAMMA_MAX_ATTEMPTS <= UINT32_MAX &&
                   SQISIGN_RANDOM_IDEAL_BETA_MAX_ATTEMPTS <= UINT32_MAX &&
                   SQISIGN_REPRESENT_INTEGER_MAX_ATTEMPTS <= UINT32_MAX,
               "norm-equation attempt caps must fit uint32_t");

/** @file
 *
 * @authors Antonin Leroux
 *
 * @brief Functions related to norm equation solving or special extremal orders
 */

/* paper Issue 14: check N^2 | nrd(x), with precondition N | nrd(x).
 * Returns 1 iff N^2 | nrd(x), 0 otherwise.
 * Implementation: nrd = N*Q + R where for each i,
 *   x[i]^2 = q[i]*N + r[i], Q = q[0]+q[1]+p*(q[2]+q[3]),
 *   R = r[0]+r[1]+p*(r[2]+r[3]). Since N | nrd and N | N*Q, also N | R.
 * Then nrd/N = Q + R/N, and (nrd/N) mod N == 0 iff N^2 | nrd.
 * Max ibz_mul transient: x[i]^2 = 2*log2(N) bit. Under paper Lemma 4N^2. */
static int
quat_alg_nrd_N2_divides(const quat_alg_elem_t *x, const ibz_t *N, const quat_alg_t *alg)
{
    int result = 0;
    ibz_t q[4], r[4], sq, acc, tmp, R_div_N, rem;
    for (int i = 0; i < 4; i++) {
        ibz_init(&q[i]);
        ibz_init(&r[i]);
    }
    ibz_init(&sq);
    ibz_init(&acc);
    ibz_init(&tmp);
    ibz_init(&R_div_N);
    ibz_init(&rem);

    if (ibz_cmp(N, &ibz_const_zero) <= 0) {
        result = -1;
        goto cleanup;
    }

    for (int i = 0; i < 4; i++) {
        ibz_mul(&sq, &(x->coord[i]), &(x->coord[i]));   /* 2*log2(N) bit transient */
        ibz_div(&q[i], &r[i], &sq, N);                  /* q[i], r[i] < N */
    }

    /* Q = q[0] + q[1] + p*(q[2] + q[3]) */
    ibz_add(&acc, &q[0], &q[1]);
    ibz_add(&tmp, &q[2], &q[3]);
    ibz_mul(&tmp, &tmp, &(alg->p));                     /* p_bits + log2(N) bit */
    ibz_add(&acc, &acc, &tmp);                          /* acc = Q */

    /* R = (r[0]+r[1]) + p*(r[2]+r[3]) */
    ibz_add(&tmp, &r[0], &r[1]);
    /* reuse R_div_N as R scratch */
    ibz_copy(&R_div_N, &tmp);
    ibz_add(&tmp, &r[2], &r[3]);
    ibz_mul(&tmp, &tmp, &(alg->p));
    ibz_add(&R_div_N, &R_div_N, &tmp);                  /* R_div_N = R */
    ibz_div(&R_div_N, &rem, &R_div_N, N);               /* R / N (rem should be 0) */
    if (!ibz_is_zero(&rem)) {
        result = -1;
        goto cleanup;
    }

    /* nrd/N = Q + R/N; check (nrd/N) mod N == 0 */
    ibz_add(&acc, &acc, &R_div_N);
    ibz_mod(&acc, &acc, N);
    result = ibz_is_zero(&acc);

cleanup:
    for (int i = 0; i < 4; i++) {
        ibz_finalize(&q[i]);
        ibz_finalize(&r[i]);
    }
    ibz_finalize(&sq);
    ibz_finalize(&acc);
    ibz_finalize(&tmp);
    ibz_finalize(&R_div_N);
    ibz_finalize(&rem);
    return result;
}

void
quat_lattice_O0_set(quat_lattice_t *O0)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(O0->basis[i][j]), 0);
        }
    }
    ibz_set(&(O0->denom), 2);
    ibz_set(&(O0->basis[0][0]), 2);
    ibz_set(&(O0->basis[1][1]), 2);
    ibz_set(&(O0->basis[2][2]), 1);
    ibz_set(&(O0->basis[1][2]), 1);
    ibz_set(&(O0->basis[3][3]), 1);
    ibz_set(&(O0->basis[0][3]), 1);
}

void
quat_lattice_O0_set_extremal(quat_p_extremal_maximal_order_t *O0)
{
    ibz_set(&O0->z.coord[1], 1);
    ibz_set(&O0->t.coord[2], 1);
    ibz_set(&O0->z.denom, 1);
    ibz_set(&O0->t.denom, 1);
    O0->q = 1;
    quat_lattice_O0_set(&(O0->order));
}

void
quat_order_elem_create(quat_alg_elem_t *elem,
                       const quat_p_extremal_maximal_order_t *order,
                       const ibz_vec_4_t *coeffs,
                       const quat_alg_t *Bpoo)
{

    // var dec
    quat_alg_elem_t quat_temp;

    // var init
    quat_alg_elem_init(&quat_temp);

    // elem = x
    quat_alg_scalar(elem, &(*coeffs)[0], &ibz_const_one);

    // quat_temp = i*y
    quat_alg_scalar(&quat_temp, &((*coeffs)[1]), &ibz_const_one);
    quat_alg_mul(&quat_temp, &order->z, &quat_temp, Bpoo);

    // elem = x + i*y
    quat_alg_add(elem, elem, &quat_temp);

    // quat_temp = z * j
    quat_alg_scalar(&quat_temp, &(*coeffs)[2], &ibz_const_one);
    quat_alg_mul(&quat_temp, &order->t, &quat_temp, Bpoo);

    // elem = x + i* + z*j
    quat_alg_add(elem, elem, &quat_temp);

    // quat_temp = t * j * i
    quat_alg_scalar(&quat_temp, &(*coeffs)[3], &ibz_const_one);
    quat_alg_mul(&quat_temp, &order->t, &quat_temp, Bpoo);
    quat_alg_mul(&quat_temp, &quat_temp, &order->z, Bpoo);

    // elem =  x + i*y + j*z + j*i*t
    quat_alg_add(elem, elem, &quat_temp);

    quat_alg_elem_finalize(&quat_temp);
}

/* The q == 1 shortcut below is valid for the canonical O0 representation,
 * whose basis columns are 1, i, (i+j)/2, (1+ij)/2.  Keep a cheap, allocation-
 * free structural guard so a caller supplying a non-canonical q == 1 order
 * retains the generic behavior instead of silently taking the shortcut. */
static NORMEQ_NOINLINE int
quat_represent_integer_has_canonical_O0(
    const quat_p_extremal_maximal_order_t *order)
{
    if (order->q != 1 ||
        ibz_is_zero(&order->z.denom) ||
        ibz_cmp_int32(&order->z.coord[0], 0) != 0 ||
        ibz_cmp(&order->z.coord[1], &order->z.denom) != 0 ||
        ibz_cmp_int32(&order->z.coord[2], 0) != 0 ||
        ibz_cmp_int32(&order->z.coord[3], 0) != 0 ||
        ibz_is_zero(&order->t.denom) ||
        ibz_cmp_int32(&order->t.coord[0], 0) != 0 ||
        ibz_cmp_int32(&order->t.coord[1], 0) != 0 ||
        ibz_cmp(&order->t.coord[2], &order->t.denom) != 0 ||
        ibz_cmp_int32(&order->t.coord[3], 0) != 0 ||
        ibz_cmp_int32(&order->order.denom, 2) != 0)
        return 0;

    return
        ibz_cmp_int32(&order->order.basis[0][0], 2) == 0 &&
        ibz_cmp_int32(&order->order.basis[0][1], 0) == 0 &&
        ibz_cmp_int32(&order->order.basis[0][2], 0) == 0 &&
        ibz_cmp_int32(&order->order.basis[0][3], 1) == 0 &&
        ibz_cmp_int32(&order->order.basis[1][0], 0) == 0 &&
        ibz_cmp_int32(&order->order.basis[1][1], 2) == 0 &&
        ibz_cmp_int32(&order->order.basis[1][2], 1) == 0 &&
        ibz_cmp_int32(&order->order.basis[1][3], 0) == 0 &&
        ibz_cmp_int32(&order->order.basis[2][0], 0) == 0 &&
        ibz_cmp_int32(&order->order.basis[2][1], 0) == 0 &&
        ibz_cmp_int32(&order->order.basis[2][2], 1) == 0 &&
        ibz_cmp_int32(&order->order.basis[2][3], 0) == 0 &&
        ibz_cmp_int32(&order->order.basis[3][0], 0) == 0 &&
        ibz_cmp_int32(&order->order.basis[3][1], 0) == 0 &&
        ibz_cmp_int32(&order->order.basis[3][2], 0) == 0 &&
        ibz_cmp_int32(&order->order.basis[3][3], 1) == 0;
}

/* Isolate the large quaternion temporary and the generic lattice-membership
 * conversion from the hot represent-integer frame. */
static NORMEQ_NOINLINE int
quat_represent_integer_make_primitive_generic(
    ibz_vec_4_t *coeffs,
    ibz_t *content,
    const quat_represent_integer_params_t *params,
    const ibz_t *expected_norm)
{
    quat_alg_elem_t elem;
    quat_alg_elem_init(&elem);
    quat_order_elem_create(&elem, params->order, coeffs, params->algebra);

#ifndef NDEBUG
    ibz_t norm_num, norm_denom;
    ibz_init(&norm_num);
    ibz_init(&norm_denom);
    assert(quat_alg_norm(&norm_num, &norm_denom, &elem, params->algebra));
    assert(ibz_is_one(&norm_denom));
    assert(ibz_cmp(&norm_num, expected_norm) == 0);
    assert(quat_lattice_contains(NULL, &params->order->order, &elem));
    ibz_finalize(&norm_denom);
    ibz_finalize(&norm_num);
#else
    (void)expected_norm;
#endif

    const int ok = quat_alg_make_primitive(
        coeffs, content, &elem, &params->order->order);
    quat_alg_elem_finalize(&elem);
    return ok;
}

typedef enum {
    QUAT_REPRESENT_INTEGER_FATAL = -1,
    QUAT_REPRESENT_INTEGER_RETRY = 0,
    QUAT_REPRESENT_INTEGER_SUCCESS = 1,
} quat_represent_integer_status_t;

/* Internal tri-state form.  The public API remains boolean, while callers
 * that themselves expose fatal/retry can distinguish a direct entropy-source
 * failure from an ordinary bounded-search miss. */
static quat_represent_integer_status_t
quat_represent_integer_with_status(
    quat_alg_elem_t *gamma,
    const ibz_t *n_gamma,
    int non_diag,
    const quat_represent_integer_params_t *params)
{

    if (gamma == NULL || n_gamma == NULL || params == NULL ||
        params->order == NULL || params->algebra == NULL ||
        params->order->q <= 0 || params->primality_test_iterations <= 0 ||
        params->primality_test_iterations >
            SQISIGN_NORMEQ_MAX_PRIMALITY_ITERATIONS ||
        ibz_cmp(&params->algebra->p, &ibz_const_zero) <= 0 ||
        ibz_cmp(n_gamma, &ibz_const_zero) <= 0 || ibz_is_even(n_gamma)) {
        return QUAT_REPRESENT_INTEGER_FATAL;
    }
    // var dec
    int found;
    quat_represent_integer_status_t status = QUAT_REPRESENT_INTEGER_RETRY;
    ibz_t cornacchia_target;
    ibz_t adjusted_n_gamma, q;
    ibz_t bound, sq_bound, temp;
#ifndef NDEBUG
    ibz_t test;
#endif
    ibz_t counter;
    ibz_vec_4_t coeffs; // coeffs = [x,y,z,t]

    if (non_diag && params->order->q % 4 != 1)
        return QUAT_REPRESENT_INTEGER_FATAL;

    // var init
    found = 0;
    ibz_init(&bound);
#ifndef NDEBUG
    ibz_init(&test);
#endif
    ibz_init(&temp);
    ibz_init(&q);
    ibz_init(&sq_bound);
    ibz_vec_4_init(&coeffs);
    ibz_init(&adjusted_n_gamma);
    ibz_init(&cornacchia_target);
    ibz_init(&counter);

    ibz_set(&q, params->order->q);

    // this could be removed in the current state
    int standard_order = (params->order->q == 1);
    const int canonical_O0 = standard_order &&
                             quat_represent_integer_has_canonical_O0(
                                 params->order);

    // adjusting the norm of gamma (multiplying by 4 to find a solution in an order of odd level)
    if (non_diag || standard_order) {
        ibz_mul(&adjusted_n_gamma, n_gamma, &ibz_const_two);
        ibz_mul(&adjusted_n_gamma, &adjusted_n_gamma, &ibz_const_two);
    } else {
        ibz_copy(&adjusted_n_gamma, n_gamma);
    }
    // computation of the first bound = sqrt (adjust_n_gamma / p - q)
    ibz_div(&sq_bound, &bound, &adjusted_n_gamma, &((params->algebra)->p));
    ibz_set(&temp, params->order->q);
    ibz_sub(&sq_bound, &sq_bound, &temp);
    if (ibz_cmp(&sq_bound, &ibz_const_zero) < 0)
        goto cleanup;
    ibz_sqrt_floor(&bound, &sq_bound);
    if (ibz_cmp(&bound, &ibz_const_one) < 0)
        goto cleanup;

    // the size of the search space is roughly n_gamma / (p√q)
    ibz_mul(&temp, &temp, &((params->algebra)->p));
    ibz_mul(&temp, &temp, &((params->algebra)->p));
    ibz_sqrt_floor(&temp, &temp);
    ibz_div(&counter, &temp, &adjusted_n_gamma, &temp);

    // entering the main loop
    for (uint32_t attempt = 0;
         !found && attempt < SQISIGN_REPRESENT_INTEGER_MAX_ATTEMPTS &&
         ibz_cmp(&counter, &ibz_const_zero) != 0;
         ++attempt) {
        // decreasing the counter
        ibz_sub(&counter, &counter, &ibz_const_one);

        // we start by sampling the first coordinate
        if (!ibz_rand_interval(&coeffs[2], &ibz_const_one, &bound)) {
            status = QUAT_REPRESENT_INTEGER_FATAL;
            goto cleanup;
        }

        // then, we sample the second coordinate
        // computing the second bound in temp as sqrt( (adjust_n_gamma - p*coeffs[2]²)/qp )
        ibz_mul(&cornacchia_target, &coeffs[2], &coeffs[2]);
        ibz_mul(&temp, &cornacchia_target, &(params->algebra->p));
        ibz_sub(&temp, &adjusted_n_gamma, &temp);
        ibz_mul(&sq_bound, &q, &(params->algebra->p));
        if (ibz_cmp(&sq_bound, &ibz_const_zero) <= 0)
            goto cleanup;
        ibz_div(&temp, &sq_bound, &temp, &sq_bound);
        if (ibz_cmp(&temp, &ibz_const_zero) < 0)
            continue;
        ibz_sqrt_floor(&temp, &temp);

        if (ibz_cmp(&temp, &ibz_const_zero) == 0) {
            continue;
        }
        // sampling the second value
        if (!ibz_rand_interval(&coeffs[3], &ibz_const_one, &temp)) {
            status = QUAT_REPRESENT_INTEGER_FATAL;
            goto cleanup;
        }

        // compute cornacchia_target = n_gamma - p * (z² + q*t²)
        ibz_mul(&temp, &coeffs[3], &coeffs[3]);
        ibz_mul(&temp, &q, &temp);
        ibz_add(&cornacchia_target, &cornacchia_target, &temp);
        ibz_mul(&cornacchia_target, &cornacchia_target, &((params->algebra)->p));
        ibz_sub(&cornacchia_target, &adjusted_n_gamma, &cornacchia_target);
        if (ibz_cmp(&cornacchia_target, &ibz_const_zero) <= 0)
            continue;

        // applying cornacchia
        if (ibz_probab_prime(&cornacchia_target, params->primality_test_iterations))
            found = ibz_cornacchia_prime(&(coeffs[0]), &(coeffs[1]), &q, &cornacchia_target);
        else
            found = 0;

        if (found && non_diag && standard_order) {
            // check that we can divide by two at least once
            // the treatmeat depends if the basis contains (1+j)/2 or (1+k)/2
            // we must have x = t mod 2 and y = z mod 2
            // if q=1 we can simply swap x and y
            if (ibz_is_odd(&coeffs[0]) != ibz_is_odd(&coeffs[3])) {
                ibz_swap(&coeffs[1], &coeffs[0]);
            }
            // we further check that (x-t)/2 = 1 mod 2 and (y-z)/2 = 1 mod 2 to ensure that the
            // resulting endomorphism will behave well for dim 2 computations
            found = found && ((ibz_get(&coeffs[0]) - ibz_get(&coeffs[3])) % 4 == 2) &&
                    ((ibz_get(&coeffs[1]) - ibz_get(&coeffs[2])) % 4 == 2);
        }
        if (found) {

#ifndef NDEBUG
            ibz_set(&temp, (params->order->q));
            ibz_mul(&temp, &temp, &(coeffs[1]));
            ibz_mul(&temp, &temp, &(coeffs[1]));
            ibz_mul(&test, &(coeffs[0]), &(coeffs[0]));
            ibz_add(&temp, &temp, &test);
            assert(0 == ibz_cmp(&temp, &cornacchia_target));

            ibz_mul(&cornacchia_target, &(coeffs[3]), &(coeffs[3]));
            ibz_mul(&cornacchia_target, &cornacchia_target, &(params->algebra->p));
            ibz_mul(&temp, &(coeffs[1]), &(coeffs[1]));
            ibz_add(&cornacchia_target, &cornacchia_target, &temp);
            ibz_set(&temp, (params->order->q));
            ibz_mul(&cornacchia_target, &cornacchia_target, &temp);
            ibz_mul(&temp, &(coeffs[0]), &coeffs[0]);
            ibz_add(&cornacchia_target, &cornacchia_target, &temp);
            ibz_mul(&temp, &(coeffs[2]), &coeffs[2]);
            ibz_mul(&temp, &temp, &(params->algebra->p));
            ibz_add(&cornacchia_target, &cornacchia_target, &temp);
            assert(0 == ibz_cmp(&cornacchia_target, &adjusted_n_gamma));
#endif
            int primitive_ok;
            if (canonical_O0) {
#ifndef NDEBUG
                /* Differentially retain the former generic conversion in
                 * debug builds.  None of this storage exists in release. */
                ibz_vec_4_t generic_coeffs;
                ibz_t generic_content;
                ibz_vec_4_init(&generic_coeffs);
                ibz_init(&generic_content);
                ibz_vec_4_copy(&generic_coeffs, &coeffs);
                const int generic_ok =
                    quat_represent_integer_make_primitive_generic(
                        &generic_coeffs,
                        &generic_content,
                        params,
                        &adjusted_n_gamma);
#endif

                /* quat_order_elem_create gives [x,y,z,-t] in the algebra
                 * basis.  In the canonical O0 basis
                 *   1, i, (i+j)/2, (1+ij)/2,
                 * its integral coordinates are
                 *   [x+t, y-z, 2z, -2t].
                 * Preserve the original z,t until the first two updates. */
                ibz_add(&coeffs[0], &coeffs[0], &coeffs[3]);
                ibz_sub(&coeffs[1], &coeffs[1], &coeffs[2]);
                ibz_add(&coeffs[2], &coeffs[2], &coeffs[2]);
                ibz_add(&coeffs[3], &coeffs[3], &coeffs[3]);
                ibz_neg(&coeffs[3], &coeffs[3]);

                ibz_vec_4_content(&temp, &coeffs);
                primitive_ok = !ibz_is_zero(&temp) &&
                               ibz_vec_4_scalar_div(
                                   &coeffs, &temp, &coeffs);

#ifndef NDEBUG
                assert(generic_ok == primitive_ok);
                assert(!primitive_ok ||
                       ibz_cmp(&generic_content, &temp) == 0);
                if (primitive_ok) {
                    for (int i = 0; i < 4; ++i)
                        assert(ibz_cmp(&generic_coeffs[i], &coeffs[i]) == 0);
                }
                ibz_finalize(&generic_content);
                ibz_vec_4_finalize(&generic_coeffs);
#endif
            } else {
                /* The generic q != 1 path (and defensive fallback for a
                 * non-canonical q == 1 order) stays behavior-identical. */
                primitive_ok =
                    quat_represent_integer_make_primitive_generic(
                        &coeffs, &temp, params, &adjusted_n_gamma);
            }

            if (!primitive_ok) {
                found = 0;
                continue;
            }

            if (non_diag || standard_order)
                found = (ibz_cmp(&temp, &ibz_const_two) == 0);
            else
                found = (ibz_cmp(&temp, &ibz_const_one) == 0);
        }
    }

    if (found) {
        /* Publish only after every stochastic and arithmetic condition has
         * succeeded.  Failed attempts therefore leave gamma unchanged. */
        ibz_mat_4x4_eval(
            &gamma->coord, &params->order->order.basis, &coeffs);
        ibz_copy(&gamma->denom, &params->order->order.denom);
        status = QUAT_REPRESENT_INTEGER_SUCCESS;
    }
cleanup:
    // var finalize
    ibz_finalize(&counter);
    ibz_finalize(&bound);
    ibz_finalize(&temp);
    ibz_finalize(&sq_bound);
    ibz_vec_4_finalize(&coeffs);
    ibz_finalize(&adjusted_n_gamma);
    ibz_finalize(&cornacchia_target);
    ibz_finalize(&q);
#ifndef NDEBUG
    ibz_finalize(&test);
#endif

    return status;
}

int
quat_represent_integer(quat_alg_elem_t *gamma,
                       const ibz_t *n_gamma,
                       int non_diag,
                       const quat_represent_integer_params_t *params)
{
    return quat_represent_integer_with_status(
               gamma, n_gamma, non_diag, params) ==
           QUAT_REPRESENT_INTEGER_SUCCESS;
}

quat_random_ideal_status_t
quat_sampling_random_ideal_O0_given_norm(quat_left_ideal_t *lideal,
                                         const ibz_t *norm,
                                         int is_prime,
                                         const quat_represent_integer_params_t *params,
                                         const ibz_t *prime_cofactor)
{

    ibz_t n_temp;
#ifndef NDEBUG
    ibz_t norm_d;
#endif
    ibz_t disc;
    quat_alg_elem_t gen, gen_rerand;
    int found = 0;
    ibz_init(&n_temp);
#ifndef NDEBUG
    ibz_init(&norm_d);
#endif
    ibz_init(&disc);
    quat_alg_elem_init(&gen);
    quat_alg_elem_init(&gen_rerand);
    quat_random_ideal_status_t status = QUAT_RANDOM_IDEAL_RETRY;

    if (lideal == NULL || norm == NULL || params == NULL ||
        params->order == NULL || params->algebra == NULL ||
        params->primality_test_iterations <= 0 ||
        params->primality_test_iterations >
            SQISIGN_NORMEQ_MAX_PRIMALITY_ITERATIONS ||
        ibz_cmp(norm, &ibz_const_one) <= 0) {
        status = QUAT_RANDOM_IDEAL_FATAL;
        goto cleanup;
    }

    // when the norm is prime we can be quite efficient
    // by avoiding to run represent integer
    // the first step is to generate one ideal of the correct norm
    if (is_prime) {

        // paper Algorithm RandomIdealGivenPrimeNorm (04Sampling.tex:38)
        // step 1: sample gamma = g1*i + g2*j + g3*k with sqrt(-nrd) mod N adjusted
        for (uint32_t attempt = 0;
             attempt < SQISIGN_RANDOM_IDEAL_GAMMA_MAX_ATTEMPTS && !found;
             ++attempt) {
            // generating a trace-zero element at random
            ibz_set(&gen.coord[0], 0);
            ibz_sub(&n_temp, norm, &ibz_const_one);
            for (int i = 1; i < 4; i++) {
                if (!ibz_rand_interval(&gen.coord[i], &ibz_const_zero, &n_temp)) {
                    status = QUAT_RANDOM_IDEAL_FATAL;
                    goto cleanup;
                }
            }

            // paper Issue 14: compute (-nrd(γ)) mod N modularly so max
            // ibz transient stays <= 2*log2(N) (vs full nrd ~ p*N^2).
            if (!quat_alg_norm_mod(&n_temp, &gen, norm, (params->algebra))) {
                status = QUAT_RANDOM_IDEAL_FATAL;
                goto cleanup;
            }

            ibz_neg(&disc, &n_temp);
            ibz_mod(&disc, &disc, norm);
            // check (-nrd:N) = 1 and compute sqrt if so
            found = ibz_sqrt_mod_p(&gen.coord[0], &disc, norm);
            found = found && !quat_alg_elem_is_zero(&gen);

            // paper while condition: gcd(nrd(gamma), N^2) = N
            // (reject if N^2 | nrd, since then gamma is in N*O_0)
            // paper Issue 14: avoid full nrd transient (~p*N^2 bit). Helper
            // decomposes nrd = N*Q + R and checks (Q + R/N) mod N == 0,
            // keeping max transient at 2*log2(N) (paper Lemma 4N^2).
            if (found) {
                int n2_divides = quat_alg_nrd_N2_divides(
                    &gen, norm, (params->algebra));
                if (n2_divides < 0) {
                    status = QUAT_RANDOM_IDEAL_FATAL;
                    goto cleanup;
                }
                if (n2_divides)
                    found = 0;
            }
        }
        if (!found)
            goto cleanup;

        // paper Algorithm steps 2-3: sample beta and set g <- gamma*beta mod N
        //   Issue 5 (author reply): mod N is REQUIRED; otherwise gamma*beta
        //   blows up modified LLL (intermediate values explode).
        {
            ibz_t beta_nrd;
            int beta_ok = 0;
            ibz_init(&beta_nrd);

            // beta = a + b*i + c*j + d*k, with a^2+b^2+p(c^2+d^2) != 0 (mod N)
            // paper Issue 14: norm-mod kept under 2*log2(N) instead of p*N^2
            for (uint32_t attempt = 0;
                 attempt < SQISIGN_RANDOM_IDEAL_BETA_MAX_ATTEMPTS && !beta_ok;
                 ++attempt) {
                ibz_sub(&n_temp, norm, &ibz_const_one);
                for (int i = 0; i < 4; i++) {
                    if (!ibz_rand_interval(
                            &gen_rerand.coord[i], &ibz_const_zero, &n_temp)) {
                        status = QUAT_RANDOM_IDEAL_FATAL;
                        ibz_finalize(&beta_nrd);
                        goto cleanup;
                    }
                }
                if (!quat_alg_norm_mod(
                        &beta_nrd, &gen_rerand, norm, (params->algebra))) {
                    status = QUAT_RANDOM_IDEAL_FATAL;
                    ibz_finalize(&beta_nrd);
                    goto cleanup;
                }
                beta_ok = (ibz_cmp(&beta_nrd, &ibz_const_zero) != 0);
            }

            if (!beta_ok) {
                ibz_finalize(&beta_nrd);
                goto cleanup;
            }

            // paper Issue 14: g = gamma*beta computed mod N*O_0 directly
            // (paper Algorithm RandomIdealGivenPrimeNorm line 59), so coord
            // intermediates stay within paper Lemma 4N^2 bound.
            if (!quat_alg_mul_mod(
                    &gen, &gen, &gen_rerand, norm, (params->algebra))) {
                status = QUAT_RANDOM_IDEAL_FATAL;
                ibz_finalize(&beta_nrd);
                goto cleanup;
            }

            ibz_finalize(&beta_nrd);
        }
    } else {
        if (prime_cofactor == NULL || ibz_is_zero(norm)) {
            status = QUAT_RANDOM_IDEAL_FATAL;
            goto cleanup;
        }
        // if it is not prime or we don't know if it is prime, we may just use represent integer
        // and use a precomputed prime as cofactor
        ibz_mul(&n_temp, prime_cofactor, norm);
        quat_represent_integer_status_t represent_status =
            quat_represent_integer_with_status(&gen, &n_temp, 0, params);
        if (represent_status == QUAT_REPRESENT_INTEGER_FATAL) {
            status = QUAT_RANDOM_IDEAL_FATAL;
            goto cleanup;
        }
        found = represent_status == QUAT_REPRESENT_INTEGER_SUCCESS;
        found = found && !quat_alg_elem_is_zero(&gen);
    }
#ifndef NDEBUG
    if (found) {
        // first, we compute the norm of the gen
        quat_alg_norm(&n_temp, &norm_d, &gen, (params->algebra));
        assert(ibz_is_one(&norm_d));
        ibz_mod(&n_temp, &n_temp, norm);
        assert(ibz_cmp(&n_temp, &ibz_const_zero) == 0);
    }
#endif

    // now we just have to rerandomize the class of the ideal generated by gen
    /*found = 0;
    while (!found) {
        for (int i = 0; i < 4; i++) {
            ibz_rand_interval(&gen_rerand.coord[i], &ibz_const_one, norm);
        }
        quat_alg_norm(&n_temp, &norm_d, &gen_rerand, (params->algebra));
        assert(ibz_is_one(&norm_d));
        ibz_gcd(&disc, &n_temp, norm);
        found = ibz_is_one(&disc);
        found = found && !quat_alg_elem_is_zero(&gen_rerand);
    }

    quat_alg_mul(&gen, &gen, &gen_rerand, (params->algebra));*/
    // in both cases, whether norm is prime or not prime,
    // gen is not divisible by any integer factor of the target norm
    // therefore the call below will yield an ideal of the correct norm
    // paper Issue 9: norm = caller-supplied N (RandomIdealGivenPrimeNorm
    // guarantees nrd(ideal) = N via Lemma nrd-mod). Skip quat_alg_norm(gen)
    // which produces a ~2*N + p_bits transient (~1275 bit at lvl1 commit).
    if (found) {
        status = quat_lideal_create_with_norm(
                     lideal, &gen, norm, &((params->order)->order), (params->algebra))
                     ? QUAT_RANDOM_IDEAL_SUCCESS
                     : QUAT_RANDOM_IDEAL_FATAL;
    }
    if (status == QUAT_RANDOM_IDEAL_SUCCESS) {
        const int norm_matches = ibz_cmp(norm, &lideal->norm) == 0;
#ifndef NDEBUG
        assert(norm_matches);
#endif
        if (!norm_matches)
            status = QUAT_RANDOM_IDEAL_FATAL;
    }

cleanup:
    ibz_finalize(&n_temp);
    quat_alg_elem_finalize(&gen);
    quat_alg_elem_finalize(&gen_rerand);
#ifndef NDEBUG
    ibz_finalize(&norm_d);
#endif
    ibz_finalize(&disc);
    return status;
}

int
quat_change_to_O0_basis(ibz_vec_4_t *vec, const quat_alg_elem_t *el)
{
    if (ibz_is_zero(&el->denom))
        return 0;
    ibz_copy(&(*vec)[2], &el->coord[2]);
    ibz_add(&(*vec)[2], &(*vec)[2], &(*vec)[2]); // double (not optimal if el->denom is even...)
    ibz_copy(&(*vec)[3], &el->coord[3]);         // double (not optimal if el->denom is even...)
    ibz_add(&(*vec)[3], &(*vec)[3], &(*vec)[3]);
    ibz_sub(&(*vec)[0], &el->coord[0], &el->coord[3]);
    ibz_sub(&(*vec)[1], &el->coord[1], &el->coord[2]);

    return ibz_vec_4_scalar_div(vec, &el->denom, vec);
}
