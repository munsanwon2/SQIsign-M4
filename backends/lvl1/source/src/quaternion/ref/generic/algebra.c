#include <quaternion.h>
#include "internal.h"

#if defined(__GNUC__) || defined(__clang__)
#define QUAT_ALG_NOINLINE __attribute__((noinline))
#else
#define QUAT_ALG_NOINLINE
#endif

// Internal helper functions

void
quat_alg_init_set_ui(quat_alg_t *alg, unsigned int p)
{
    ibz_t bp;
    ibz_init(&bp);
    ibz_set(&bp, p);
    quat_alg_init_set(alg, &bp);
    ibz_finalize(&bp);
}

void
quat_alg_coord_mul(ibz_vec_4_t *res, const ibz_vec_4_t *a, const ibz_vec_4_t *b, const quat_alg_t *alg)
{
    ibz_t prod;
    ibz_vec_4_t sum;
    ibz_init(&prod);
    ibz_vec_4_init(&sum);

    ibz_set(&(sum[0]), 0);
    ibz_set(&(sum[1]), 0);
    ibz_set(&(sum[2]), 0);
    ibz_set(&(sum[3]), 0);

    // compute 1 coordinate
    ibz_mul(&prod, &((*a)[2]), &((*b)[2]));
    ibz_sub(&(sum[0]), &(sum[0]), &prod);
    ibz_mul(&prod, &((*a)[3]), &((*b)[3]));
    ibz_sub(&(sum[0]), &(sum[0]), &prod);
    ibz_mul(&(sum[0]), &(sum[0]), &(alg->p));
    ibz_mul(&prod, &((*a)[0]), &((*b)[0]));
    ibz_add(&(sum[0]), &(sum[0]), &prod);
    ibz_mul(&prod, &((*a)[1]), &((*b)[1]));
    ibz_sub(&(sum[0]), &(sum[0]), &prod);
    // compute i coordiante
    ibz_mul(&prod, &((*a)[2]), &((*b)[3]));
    ibz_add(&(sum[1]), &(sum[1]), &prod);
    ibz_mul(&prod, &((*a)[3]), &((*b)[2]));
    ibz_sub(&(sum[1]), &(sum[1]), &prod);
    ibz_mul(&(sum[1]), &(sum[1]), &(alg->p));
    ibz_mul(&prod, &((*a)[0]), &((*b)[1]));
    ibz_add(&(sum[1]), &(sum[1]), &prod);
    ibz_mul(&prod, &((*a)[1]), &((*b)[0]));
    ibz_add(&(sum[1]), &(sum[1]), &prod);
    // compute j coordiante
    ibz_mul(&prod, &((*a)[0]), &((*b)[2]));
    ibz_add(&(sum[2]), &(sum[2]), &prod);
    ibz_mul(&prod, &((*a)[2]), &((*b)[0]));
    ibz_add(&(sum[2]), &(sum[2]), &prod);
    ibz_mul(&prod, &((*a)[1]), &((*b)[3]));
    ibz_sub(&(sum[2]), &(sum[2]), &prod);
    ibz_mul(&prod, &((*a)[3]), &((*b)[1]));
    ibz_add(&(sum[2]), &(sum[2]), &prod);
    // compute ij coordiante
    ibz_mul(&prod, &((*a)[0]), &((*b)[3]));
    ibz_add(&(sum[3]), &(sum[3]), &prod);
    ibz_mul(&prod, &((*a)[3]), &((*b)[0]));
    ibz_add(&(sum[3]), &(sum[3]), &prod);
    ibz_mul(&prod, &((*a)[2]), &((*b)[1]));
    ibz_sub(&(sum[3]), &(sum[3]), &prod);
    ibz_mul(&prod, &((*a)[1]), &((*b)[2]));
    ibz_add(&(sum[3]), &(sum[3]), &prod);

    ibz_copy(&((*res)[0]), &(sum[0]));
    ibz_copy(&((*res)[1]), &(sum[1]));
    ibz_copy(&((*res)[2]), &(sum[2]));
    ibz_copy(&((*res)[3]), &(sum[3]));

    ibz_finalize(&prod);
    ibz_vec_4_finalize(&sum);
}

void
quat_alg_equal_denom(quat_alg_elem_t *res_a, quat_alg_elem_t *res_b, const quat_alg_elem_t *a, const quat_alg_elem_t *b)
{
    ibz_t gcd, r;
    ibz_init(&gcd);
    ibz_init(&r);
    ibz_gcd(&gcd, &(a->denom), &(b->denom));
    // temporarily set res_a.denom to a.denom/gcd, and res_b.denom to b.denom/gcd
    ibz_div(&(res_a->denom), &r, &(a->denom), &gcd);
    ibz_div(&(res_b->denom), &r, &(b->denom), &gcd);
    for (int i = 0; i < 4; i++) {
        // multiply coordiates by reduced denominators from the other element
        ibz_mul(&(res_a->coord[i]), &(a->coord[i]), &(res_b->denom));
        ibz_mul(&(res_b->coord[i]), &(b->coord[i]), &(res_a->denom));
    }
    // multiply both reduced denominators
    ibz_mul(&(res_a->denom), &(res_a->denom), &(res_b->denom));
    // multiply them by the gcd to get the new common denominator
    ibz_mul(&(res_b->denom), &(res_a->denom), &gcd);
    ibz_mul(&(res_a->denom), &(res_a->denom), &gcd);
    ibz_finalize(&gcd);
    ibz_finalize(&r);
}

// Public Functions

static QUAT_ALG_NOINLINE void
quat_alg_add_different_denom(quat_alg_elem_t *res,
                             const quat_alg_elem_t *a,
                             const quat_alg_elem_t *b)
{
    quat_alg_elem_t res_a, res_b;
    quat_alg_elem_init(&res_a);
    quat_alg_elem_init(&res_b);
    quat_alg_equal_denom(&res_a, &res_b, a, b);
    ibz_copy(&res->denom, &res_a.denom);
    ibz_vec_4_add(&res->coord, &res_a.coord, &res_b.coord);
    quat_alg_elem_finalize(&res_a);
    quat_alg_elem_finalize(&res_b);
}

static QUAT_ALG_NOINLINE void
quat_alg_sub_different_denom(quat_alg_elem_t *res,
                             const quat_alg_elem_t *a,
                             const quat_alg_elem_t *b)
{
    quat_alg_elem_t res_a, res_b;
    quat_alg_elem_init(&res_a);
    quat_alg_elem_init(&res_b);
    quat_alg_equal_denom(&res_a, &res_b, a, b);
    ibz_copy(&res->denom, &res_a.denom);
    ibz_vec_4_sub(&res->coord, &res_a.coord, &res_b.coord);
    quat_alg_elem_finalize(&res_a);
    quat_alg_elem_finalize(&res_b);
}

void
quat_alg_add(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_elem_t *b)
{
    if (ibz_cmp(&a->denom, &b->denom) == 0) {
        ibz_copy(&res->denom, &a->denom);
        ibz_vec_4_add(&res->coord, &a->coord, &b->coord);
        return;
    }
    quat_alg_add_different_denom(res, a, b);
}

void
quat_alg_sub(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_elem_t *b)
{
    if (ibz_cmp(&a->denom, &b->denom) == 0) {
        ibz_copy(&res->denom, &a->denom);
        ibz_vec_4_sub(&res->coord, &a->coord, &b->coord);
        return;
    }
    quat_alg_sub_different_denom(res, a, b);
}

void
quat_alg_mul(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_elem_t *b, const quat_alg_t *alg)
{
    // denominator: product of denominators
    ibz_mul(&(res->denom), &(a->denom), &(b->denom));
    quat_alg_coord_mul(&(res->coord), &(a->coord), &(b->coord), alg);
}

int
quat_alg_norm_mod(ibz_t *res, const quat_alg_elem_t *x, const ibz_t *N, const quat_alg_t *alg)
{
    /* paper Issue 14: nrd(x) mod N = (x0^2 + x1^2 + p*x2^2 + p*x3^2) mod N
     * with each ibz_mul reduced before the next operation. Max transient
     * is max(2*log2(N), p_bits + log2(N)), matching paper Lemma 4N^2 when N >= p. */
    if (!ibz_is_one(&(x->denom)) || ibz_cmp(N, &ibz_const_zero) <= 0) {
        ibz_set(res, 0);
        return 0;
    }
    ibz_t tmp;
    ibz_init(&tmp);

    /* res = x0^2 mod N */
    ibz_mul(res, &(x->coord[0]), &(x->coord[0]));
    ibz_mod(res, res, N);

    /* res += x1^2 mod N */
    ibz_mul(&tmp, &(x->coord[1]), &(x->coord[1]));
    ibz_mod(&tmp, &tmp, N);
    ibz_add(res, res, &tmp);

    /* res += p * (x2^2 mod N) mod N */
    ibz_mul(&tmp, &(x->coord[2]), &(x->coord[2]));
    ibz_mod(&tmp, &tmp, N);
    ibz_mul(&tmp, &tmp, &(alg->p));
    ibz_mod(&tmp, &tmp, N);
    ibz_add(res, res, &tmp);

    /* res += p * (x3^2 mod N) mod N */
    ibz_mul(&tmp, &(x->coord[3]), &(x->coord[3]));
    ibz_mod(&tmp, &tmp, N);
    ibz_mul(&tmp, &tmp, &(alg->p));
    ibz_mod(&tmp, &tmp, N);
    ibz_add(res, res, &tmp);

    ibz_mod(res, res, N);

    ibz_finalize(&tmp);
    return 1;
}

int
quat_alg_mul_mod(quat_alg_elem_t *res, const quat_alg_elem_t *a, const quat_alg_elem_t *b, const ibz_t *N, const quat_alg_t *alg)
{
    /* paper Issue 14: (a*b) mod N*O_0 — coords computed with each ibz_mul
     * followed by ibz_mod, max transient = max(2*log2(N), p_bits + log2(N)).
     * Use a local sum[4] so aliasing res==a or res==b is safe. */
    if (!ibz_is_one(&(a->denom)) || !ibz_is_one(&(b->denom)) ||
        ibz_cmp(N, &ibz_const_zero) <= 0)
        return 0;
    ibz_t prod;
    ibz_vec_4_t sum;
    ibz_init(&prod);
    ibz_vec_4_init(&sum);

    /* sum[0] = a[0]*b[0] - a[1]*b[1] - p*(a[2]*b[2] + a[3]*b[3]) */
    ibz_mul(&prod, &(a->coord[2]), &(b->coord[2]));
    ibz_mod(&prod, &prod, N);
    ibz_sub(&sum[0], &sum[0], &prod);
    ibz_mul(&prod, &(a->coord[3]), &(b->coord[3]));
    ibz_mod(&prod, &prod, N);
    ibz_sub(&sum[0], &sum[0], &prod);
    ibz_mod(&sum[0], &sum[0], N);
    ibz_mul(&sum[0], &sum[0], &(alg->p));
    ibz_mod(&sum[0], &sum[0], N);
    ibz_mul(&prod, &(a->coord[0]), &(b->coord[0]));
    ibz_mod(&prod, &prod, N);
    ibz_add(&sum[0], &sum[0], &prod);
    ibz_mul(&prod, &(a->coord[1]), &(b->coord[1]));
    ibz_mod(&prod, &prod, N);
    ibz_sub(&sum[0], &sum[0], &prod);
    ibz_mod(&sum[0], &sum[0], N);

    /* sum[1] = a[0]*b[1] + a[1]*b[0] + p*(a[2]*b[3] - a[3]*b[2]) */
    ibz_mul(&prod, &(a->coord[2]), &(b->coord[3]));
    ibz_mod(&prod, &prod, N);
    ibz_add(&sum[1], &sum[1], &prod);
    ibz_mul(&prod, &(a->coord[3]), &(b->coord[2]));
    ibz_mod(&prod, &prod, N);
    ibz_sub(&sum[1], &sum[1], &prod);
    ibz_mod(&sum[1], &sum[1], N);
    ibz_mul(&sum[1], &sum[1], &(alg->p));
    ibz_mod(&sum[1], &sum[1], N);
    ibz_mul(&prod, &(a->coord[0]), &(b->coord[1]));
    ibz_mod(&prod, &prod, N);
    ibz_add(&sum[1], &sum[1], &prod);
    ibz_mul(&prod, &(a->coord[1]), &(b->coord[0]));
    ibz_mod(&prod, &prod, N);
    ibz_add(&sum[1], &sum[1], &prod);
    ibz_mod(&sum[1], &sum[1], N);

    /* sum[2] = a[0]*b[2] + a[2]*b[0] - a[1]*b[3] + a[3]*b[1] */
    ibz_mul(&prod, &(a->coord[0]), &(b->coord[2]));
    ibz_mod(&prod, &prod, N);
    ibz_add(&sum[2], &sum[2], &prod);
    ibz_mul(&prod, &(a->coord[2]), &(b->coord[0]));
    ibz_mod(&prod, &prod, N);
    ibz_add(&sum[2], &sum[2], &prod);
    ibz_mul(&prod, &(a->coord[1]), &(b->coord[3]));
    ibz_mod(&prod, &prod, N);
    ibz_sub(&sum[2], &sum[2], &prod);
    ibz_mul(&prod, &(a->coord[3]), &(b->coord[1]));
    ibz_mod(&prod, &prod, N);
    ibz_add(&sum[2], &sum[2], &prod);
    ibz_mod(&sum[2], &sum[2], N);

    /* sum[3] = a[0]*b[3] + a[3]*b[0] - a[2]*b[1] + a[1]*b[2] */
    ibz_mul(&prod, &(a->coord[0]), &(b->coord[3]));
    ibz_mod(&prod, &prod, N);
    ibz_add(&sum[3], &sum[3], &prod);
    ibz_mul(&prod, &(a->coord[3]), &(b->coord[0]));
    ibz_mod(&prod, &prod, N);
    ibz_add(&sum[3], &sum[3], &prod);
    ibz_mul(&prod, &(a->coord[2]), &(b->coord[1]));
    ibz_mod(&prod, &prod, N);
    ibz_sub(&sum[3], &sum[3], &prod);
    ibz_mul(&prod, &(a->coord[1]), &(b->coord[2]));
    ibz_mod(&prod, &prod, N);
    ibz_add(&sum[3], &sum[3], &prod);
    ibz_mod(&sum[3], &sum[3], N);

    for (int i = 0; i < 4; i++)
        ibz_copy(&(res->coord[i]), &sum[i]);
    ibz_set(&(res->denom), 1);

    ibz_finalize(&prod);
    ibz_vec_4_finalize(&sum);
    return 1;
}

int
quat_alg_norm(ibz_t *res_num, ibz_t *res_denom, const quat_alg_elem_t *a, const quat_alg_t *alg)
{
    int ok = 0;
    ibz_t r, g, num, denom;
    ibz_init(&r);
    ibz_init(&g);
    ibz_init(&num);
    ibz_init(&denom);

    if (res_num == NULL || res_denom == NULL || a == NULL || alg == NULL ||
        ibz_is_zero(&a->denom) || ibz_cmp(&alg->p, &ibz_const_zero) <= 0)
        goto cleanup;

    /* nrd(x0+x1*i+x2*j+x3*ij) =
     * x0^2+x1^2+p*(x2^2+x3^2).  This is exactly the scalar coordinate of
     * a*conj(a), without constructing a conjugate or the three identically
     * zero imaginary coordinates.  Keep the former scalar-coordinate order
     * so fixed-precision intermediate bounds are unchanged. */
    ibz_mul(&num, &a->coord[2], &a->coord[2]);
    ibz_mul(&r, &a->coord[3], &a->coord[3]);
    ibz_add(&num, &num, &r);
    ibz_mul(&num, &num, &alg->p);
    ibz_mul(&r, &a->coord[0], &a->coord[0]);
    ibz_add(&num, &num, &r);
    ibz_mul(&r, &a->coord[1], &a->coord[1]);
    ibz_add(&num, &num, &r);
    ibz_mul(&denom, &a->denom, &a->denom);

    ibz_gcd(&g, &num, &denom);
    if (ibz_is_zero(&g))
        goto cleanup;
    ibz_div(&num, &r, &num, &g);
    if (!ibz_is_zero(&r))
        goto cleanup;
    ibz_div(&denom, &r, &denom, &g);
    if (!ibz_is_zero(&r) || ibz_is_zero(&denom))
        goto cleanup;
    ibz_abs(&denom, &denom);
    ibz_abs(&num, &num);
    ibz_copy(res_num, &num);
    ibz_copy(res_denom, &denom);
    ok = 1;

cleanup:
    if (!ok) {
        if (res_num != NULL)
            ibz_set(res_num, 0);
        if (res_denom != NULL)
            ibz_set(res_denom, 0);
    }
    ibz_finalize(&denom);
    ibz_finalize(&num);
    ibz_finalize(&r);
    ibz_finalize(&g);
    return ok;
}

void
quat_alg_scalar(quat_alg_elem_t *elem, const ibz_t *numerator, const ibz_t *denominator)
{
    ibz_copy(&(elem->denom), denominator);
    ibz_copy(&(elem->coord[0]), numerator);
    ibz_set(&(elem->coord[1]), 0);
    ibz_set(&(elem->coord[2]), 0);
    ibz_set(&(elem->coord[3]), 0);
}

void
quat_alg_conj(quat_alg_elem_t *conj, const quat_alg_elem_t *x)
{
    ibz_copy(&(conj->denom), &(x->denom));
    ibz_copy(&(conj->coord[0]), &(x->coord[0]));
    ibz_neg(&(conj->coord[1]), &(x->coord[1]));
    ibz_neg(&(conj->coord[2]), &(x->coord[2]));
    ibz_neg(&(conj->coord[3]), &(x->coord[3]));
}

int
quat_alg_make_primitive(ibz_vec_4_t *primitive_x, ibz_t *content, const quat_alg_elem_t *x, const quat_lattice_t *order)
{
    if (!quat_lattice_contains(primitive_x, order, x)) {
        ibz_set(content, 0);
        return 0;
    }
    ibz_vec_4_content(content, primitive_x);
    if (ibz_is_zero(content) ||
        !ibz_vec_4_scalar_div(primitive_x, content, primitive_x))
        return 0;
    return 1;
}

void
quat_alg_normalize(quat_alg_elem_t *x)
{
    ibz_t gcd;
    ibz_init(&gcd);
    ibz_vec_4_content(&gcd, &(x->coord));
    ibz_gcd(&gcd, &gcd, &(x->denom));
    if (!ibz_is_one(&gcd)) {
        ibz_div(&(x->denom), NULL, &(x->denom), &gcd);
        (void)ibz_vec_4_scalar_div(&(x->coord), &gcd, &(x->coord));
    }
    if (ibz_is_negative(&x->denom)) {
        ibz_neg(&x->denom, &x->denom);
        ibz_vec_4_negate(&x->coord, &x->coord);
    }
    ibz_finalize(&gcd);
}

int
quat_alg_elem_equal(const quat_alg_elem_t *a, const quat_alg_elem_t *b)
{
    quat_alg_elem_t diff;
    quat_alg_elem_init(&diff);
    quat_alg_sub(&diff, a, b);
    int res = quat_alg_elem_is_zero(&diff);
    quat_alg_elem_finalize(&diff);
    return (res);
}

int
quat_alg_elem_is_zero(const quat_alg_elem_t *x)
{
    int res = ibz_vec_4_is_zero(&(x->coord));
    return (res);
}

void
quat_alg_elem_set(quat_alg_elem_t *elem, int32_t denom, int32_t coord0, int32_t coord1, int32_t coord2, int32_t coord3)
{
    ibz_set(&(elem->coord[0]), coord0);
    ibz_set(&(elem->coord[1]), coord1);
    ibz_set(&(elem->coord[2]), coord2);
    ibz_set(&(elem->coord[3]), coord3);

    ibz_set(&(elem->denom), denom);
}

void
quat_alg_elem_copy(quat_alg_elem_t *copy, const quat_alg_elem_t *copied)
{
    ibz_copy(&copy->denom, &copied->denom);
    ibz_copy(&copy->coord[0], &copied->coord[0]);
    ibz_copy(&copy->coord[1], &copied->coord[1]);
    ibz_copy(&copy->coord[2], &copied->coord[2]);
    ibz_copy(&copy->coord[3], &copied->coord[3]);
}

// helper functions for lattices
void
quat_alg_elem_copy_ibz(quat_alg_elem_t *elem,
                       const ibz_t *denom,
                       const ibz_t *coord0,
                       const ibz_t *coord1,
                       const ibz_t *coord2,
                       const ibz_t *coord3)
{
    ibz_copy(&(elem->coord[0]), coord0);
    ibz_copy(&(elem->coord[1]), coord1);
    ibz_copy(&(elem->coord[2]), coord2);
    ibz_copy(&(elem->coord[3]), coord3);

    ibz_copy(&(elem->denom), denom);
}

void
quat_alg_elem_mul_by_scalar(quat_alg_elem_t *res, const ibz_t *scalar, const quat_alg_elem_t *elem)
{
    for (int i = 0; i < 4; i++) {
        ibz_mul(&(res->coord[i]), &(elem->coord[i]), scalar);
    }
    ibz_copy(&(res->denom), &(elem->denom));
}
