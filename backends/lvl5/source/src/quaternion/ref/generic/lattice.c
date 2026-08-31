#include <quaternion.h>
#include <rng.h>
#include "internal.h"
#include "lll_internals.h" /* for quat_ml2 (paper Issue 8 wire) */

static int
quat_lattice_add_products_fit(const quat_lattice_t *lat1,
                              const quat_lattice_t *lat2)
{
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

/* ML2 exposes basis vectors as rows, while lattice matrices store them as
 * columns.  Six limb-array swaps avoid ibz_mat_4x4_transpose's full matrix
 * temporary and are safe for an in-place candidate. */
static void
quat_lattice_basis_transpose_in_place(ibz_mat_4x4_t *basis)
{
    for (int row = 0; row < 4; row++)
        for (int column = row + 1; column < 4; column++)
            ibz_swap(&(*basis)[row][column], &(*basis)[column][row]);
}

// helper functions
int
quat_lattice_equal(const quat_lattice_t *lat1, const quat_lattice_t *lat2)
{
    int equal = 1;
    quat_lattice_t a, b;
    quat_lattice_init(&a);
    quat_lattice_init(&b);
    if (!quat_lattice_reduce_denom(&a, lat1) ||
        !quat_lattice_reduce_denom(&b, lat2)) {
        equal = 0;
        goto cleanup;
    }
    ibz_abs(&(a.denom), &(a.denom));
    ibz_abs(&(b.denom), &(b.denom));
    if (!quat_lattice_hnf(&a) || !quat_lattice_hnf(&b)) {
        equal = 0;
        goto cleanup;
    }
    equal = equal && (ibz_cmp(&(a.denom), &(b.denom)) == 0);
    equal = equal && ibz_mat_4x4_equal(&(a.basis), &(b.basis));
cleanup:
    quat_lattice_finalize(&a);
    quat_lattice_finalize(&b);
    return (equal);
}

// sublattice test
int
quat_lattice_inclusion(const quat_lattice_t *sublat, const quat_lattice_t *overlat)
{
    int res;
    quat_lattice_t sum;
    quat_lattice_init(&sum);
    res = quat_lattice_add(&sum, overlat, sublat) && quat_lattice_equal(&sum, overlat);
    quat_lattice_finalize(&sum);
    return (res);
}

int
quat_lattice_reduce_denom(quat_lattice_t *reduced, const quat_lattice_t *lat)
{
    int ok = 0;
    ibz_t gcd, remainder, candidate_denom;
    ibz_mat_4x4_t candidate_basis;
    ibz_init(&gcd);
    ibz_init(&remainder);
    ibz_init(&candidate_denom);
    ibz_mat_4x4_init(&candidate_basis);
    ibz_mat_4x4_gcd(&gcd, &(lat->basis));
    ibz_gcd(&gcd, &gcd, &(lat->denom));
    if (ibz_is_zero(&gcd) ||
        !ibz_mat_4x4_scalar_div(&candidate_basis, &gcd, &(lat->basis)))
        goto cleanup;
    ibz_div(&candidate_denom, &remainder, &(lat->denom), &gcd);
    if (!ibz_is_zero(&remainder) || ibz_is_zero(&candidate_denom))
        goto cleanup;
    ibz_abs(&candidate_denom, &candidate_denom);
    ibz_mat_4x4_copy(&reduced->basis, &candidate_basis);
    ibz_copy(&reduced->denom, &candidate_denom);
    ok = 1;

cleanup:
    ibz_mat_4x4_finalize(&candidate_basis);
    ibz_finalize(&candidate_denom);
    ibz_finalize(&remainder);
    ibz_finalize(&gcd);
    return ok;
}

void
quat_lattice_conjugate_without_hnf(quat_lattice_t *conj, const quat_lattice_t *lat)
{
    ibz_mat_4x4_copy(&(conj->basis), &(lat->basis));
    ibz_copy(&(conj->denom), &(lat->denom));

    for (int row = 1; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            ibz_neg(&(conj->basis[row][col]), &(conj->basis[row][col]));
        }
    }
}

// Method described in https://cseweb.ucsd.edu/classes/sp14/cse206A-a/lec4.pdf consulted on 19 of
// May 2023, 12h40 CEST
int
quat_lattice_dual_without_hnf(quat_lattice_t *dual, const quat_lattice_t *lat)
{
    /* paper Algorithm LatticeDual (91proof.tex:435): (A, Delta) = adj+det,
     * M^# = r * A^T, r^# = Delta, then GCD-normalize (Lines 7-8). The GCD
     * step is essential for the paper bound lem:dual-bound; without it,
     * dual.basis can grow to 4*|entry| bits (cofactor expansion). */
    if (dual == NULL || lat == NULL || ibz_is_zero(&lat->denom))
        return 0;

    int ok = 0;
    ibz_mat_4x4_t inv;
    ibz_t det;
    quat_lattice_t candidate;
    ibz_init(&det);
    ibz_mat_4x4_init(&inv);
    quat_lattice_init(&candidate);
    if (!ibz_mat_4x4_inv_with_det_as_denom(&inv, &det, &(lat->basis))) {
        goto cleanup;
    }
    ibz_mat_4x4_transpose(&inv, &inv);
    // dual_denom = det/lat_denom
    ibz_mat_4x4_scalar_mul(&candidate.basis, &(lat->denom), &inv);
    ibz_copy(&candidate.denom, &det);

    /* paper Line 7-8: GCD normalization */
    if (!quat_lattice_reduce_denom(&candidate, &candidate))
        goto cleanup;

    ibz_mat_4x4_copy(&dual->basis, &candidate.basis);
    ibz_copy(&dual->denom, &candidate.denom);
    ok = 1;

cleanup:
    quat_lattice_finalize(&candidate);
    ibz_finalize(&det);
    ibz_mat_4x4_finalize(&inv);
    return ok;
}

int
quat_lattice_add(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t *lat2)
{
    ibz_vec_4_t generators[8];
    quat_lattice_t candidate;
    int ok = 0;
    if (!quat_lattice_add_products_fit(lat1, lat2))
        return 0;
    for (int i = 0; i < 8; i++)
        ibz_vec_4_init(&(generators[i]));
    quat_lattice_init(&candidate);

    /* The candidate is unpublished until the complete operation succeeds, so
     * its basis can also serve as the short-lived scaled-matrix scratch. */
    ibz_mat_4x4_scalar_mul(
        &candidate.basis, &(lat1->denom), &(lat2->basis));
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_copy(&(generators[j][i]), &(candidate.basis[i][j]));
        }
    }
    ibz_mat_4x4_scalar_mul(
        &candidate.basis, &(lat2->denom), &(lat1->basis));
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_copy(&(generators[4 + j][i]), &(candidate.basis[i][j]));
        }
    }
    /* Compact arithmetic uses ML2 directly on the eight generators.  The old
     * HNF fallback required eager 4x4 determinants; those cofactors can exceed
     * the fixed-width budget even when ML2's Gram matrix fits. */
    {
        int rho = quat_ml2_retry(
            &candidate.basis[0], 4, generators, 8, NULL);
        if (rho == 4) {
            quat_lattice_basis_transpose_in_place(&candidate.basis);
            ibz_mul(&candidate.denom, &(lat1->denom), &(lat2->denom));
            if (quat_lattice_reduce_denom(&candidate, &candidate)) {
                ibz_mat_4x4_copy(&res->basis, &candidate.basis);
                ibz_copy(&res->denom, &candidate.denom);
                ok = 1;
            }
        }
    }
    quat_lattice_finalize(&candidate);
    for (int i = 0; i < 8; i++)
        ibz_vec_4_finalize(&(generators[i]));
    return ok;
}

// method described in https://cseweb.ucsd.edu/classes/sp14/cse206A-a/lec4.pdf consulted on 19 of
// May 2023, 12h40 CEST
int
quat_lattice_intersect(quat_lattice_t *res, const quat_lattice_t *lat1, const quat_lattice_t *lat2)
{
    int ok = 0;
    quat_lattice_t dual1, dual2, dual_res;
    quat_lattice_init(&dual1);
    quat_lattice_init(&dual2);
    quat_lattice_init(&dual_res);
    if (!quat_lattice_dual_without_hnf(&dual1, lat1) ||
        !quat_lattice_dual_without_hnf(&dual2, lat2))
        goto cleanup;
    if (!quat_lattice_add(&dual_res, &dual1, &dual2))
        goto cleanup;
    if (!quat_lattice_dual_without_hnf(res, &dual_res))
        goto cleanup;
    /* paper Issue 8 (2026-05-17): removed the redundant final HNF.
     * The lattice basis is already a valid (non-HNF) basis; consumers via
     * sample_response only need a basis matrix, not a unique form. Original
     * comment in upstream: "could be removed if we do not expect HNF any more". */
    /* quat_lattice_hnf(res); */
    ok = 1;
cleanup:
    quat_lattice_finalize(&dual1);
    quat_lattice_finalize(&dual2);
    quat_lattice_finalize(&dual_res);
    return ok;
}

void
quat_lattice_mat_alg_coord_mul_without_hnf(ibz_mat_4x4_t *prod,
                                           const ibz_mat_4x4_t *lat,
                                           const ibz_vec_4_t *coord,
                                           const quat_alg_t *alg)
{
    ibz_vec_4_t p, a;
    ibz_vec_4_init(&p);
    ibz_vec_4_init(&a);
    for (int i = 0; i < 4; i++) {
        ibz_vec_4_copy_ibz(&a, &((*lat)[0][i]), &((*lat)[1][i]), &((*lat)[2][i]), &((*lat)[3][i]));
        quat_alg_coord_mul(&p, &a, coord, alg);
        ibz_copy(&((*prod)[0][i]), &(p[0]));
        ibz_copy(&((*prod)[1][i]), &(p[1]));
        ibz_copy(&((*prod)[2][i]), &(p[2]));
        ibz_copy(&((*prod)[3][i]), &(p[3]));
    }
    ibz_vec_4_finalize(&p);
    ibz_vec_4_finalize(&a);
}

int
quat_lattice_alg_elem_mul(quat_lattice_t *prod,
                          const quat_lattice_t *lat,
                          const quat_alg_elem_t *elem,
                          const quat_alg_t *alg)
{
    quat_lattice_t candidate;
    int ok = 0;
    quat_lattice_init(&candidate);
    quat_lattice_mat_alg_coord_mul_without_hnf(
        &candidate.basis, &lat->basis, &elem->coord, alg);
    ibz_mul(&candidate.denom, &lat->denom, &elem->denom);
    /* paper Issue 8 extension (replaces quat_lattice_hnf): use ML2 on the
     * 4 column generators to obtain a LLL-reduced basis of the lattice
     * (without cofactor blow-up). Transpose the private candidate so its rows
     * are the four generators, then use ML2's documented input/output overlap:
     * the complete input is copied before rank-four publication. */
    {
        quat_lattice_basis_transpose_in_place(&candidate.basis);
        int rho = quat_ml2_retry(&candidate.basis[0],
                                 4,
                                 (const ibz_vec_4_t *)&candidate.basis[0],
                                 4,
                                 alg);
        if (rho != 4)
            goto cleanup;
        quat_lattice_basis_transpose_in_place(&candidate.basis);
    }
    if (!quat_lattice_reduce_denom(&candidate, &candidate))
        goto cleanup;
    ibz_mat_4x4_copy(&prod->basis, &candidate.basis);
    ibz_copy(&prod->denom, &candidate.denom);
    ok = 1;

cleanup:
    quat_lattice_finalize(&candidate);
    return ok;
}

int
quat_lattice_mul(quat_lattice_t *res,
                 const quat_lattice_t *lat1,
                 const quat_lattice_t *lat2,
                 const quat_alg_t *alg,
                 const ibz_t *norm1,
                 const ibz_t *norm2)
{
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
    int p_bits = ibz_bitsize(&alg->p);
    int generator_product_bound = basis1_bits + basis2_bits + p_bits + 3;
    int denom1_bits = ibz_bitsize(&lat1->denom);
    int denom2_bits = ibz_bitsize(&lat2->denom);
    int hnfmod_bound = 2 * ibz_bitsize(norm1) + 2 * ibz_bitsize(norm2) +
                       4 * denom1_bits + 4 * denom2_bits;
    int safety_margin = 64;
    int use_ml2 = hnfmod_bound > IBZ_HNF_ROUTE_BITS - 1 ||
                  2 * hnfmod_bound + safety_margin > IBZ_HNF_ROUTE_BITS;

    if (ibz_cmp(norm1, &ibz_const_zero) <= 0 ||
        ibz_cmp(norm2, &ibz_const_zero) <= 0 ||
        generator_product_bound > IBZ_BITS - 1 ||
        denom1_bits + denom2_bits > IBZ_BITS - 1)
        return 0;

    ibz_vec_4_t elem1, elem2, elem_res;
    ibz_vec_4_t generators[16];
    // ibz_mat_4x4_t detmat;
    // ibz_t det;
    ibz_t hnfmod, r1, r2;
    quat_lattice_t lat_res;
    int ok = 0;
    // ibz_init(&det);
    // ibz_mat_4x4_init(&detmat);
    ibz_init(&hnfmod);
    ibz_init(&r1);
    ibz_init(&r2);
    quat_lattice_init(&lat_res);
    ibz_vec_4_init(&elem1);
    ibz_vec_4_init(&elem2);
    ibz_vec_4_init(&elem_res);
    for (int i = 0; i < 16; i++)
        ibz_vec_4_init(&(generators[i]));
    for (int k = 0; k < 4; k++) {
        ibz_vec_4_copy_ibz(
            &elem1, &(lat1->basis[0][k]), &(lat1->basis[1][k]), &(lat1->basis[2][k]), &(lat1->basis[3][k]));
        for (int i = 0; i < 4; i++) {
            ibz_vec_4_copy_ibz(
                &elem2, &(lat2->basis[0][i]), &(lat2->basis[1][i]), &(lat2->basis[2][i]), &(lat2->basis[3][i]));
            quat_alg_coord_mul(&elem_res, &elem1, &elem2, alg);
            for (int j = 0; j < 4; j++) {
                /*if (k == 0)
                    ibz_copy(&(detmat[i][j]), &(elem_res[j]));*/
                ibz_copy(&(generators[4 * k + i][j]), &(elem_res[j]));
            }
        }
    }
    /* Compute the HNF modulus only when its complete multiplication chain is
     * known to fit.  The ML2 branch does not need this value. */
    if (!use_ml2) {
        ibz_mul(&hnfmod, norm1, norm1);
        ibz_mul(&hnfmod, &hnfmod, norm2);
        ibz_mul(&hnfmod, &hnfmod, norm2);
        ibz_mul(&r1, &(lat1->denom), &(lat1->denom));
        ibz_mul(&r1, &r1, &r1);
        ibz_mul(&r2, &(lat2->denom), &(lat2->denom));
        ibz_mul(&r2, &r2, &r2);
        ibz_mul(&hnfmod, &hnfmod, &r1);
        ibz_mul(&hnfmod, &hnfmod, &r2);
    }
    // ibz_t t1, t2;
    // ibz_init(&t1);
    // ibz_init(&t2);

    // // hnfmod = norm1^2 * norm2^2
    // ibz_mul(&t1, norm1, norm1);   // t1 = norm1^2
    // ibz_mul(&t2, norm2, norm2);   // t2 = norm2^2
    // ibz_mul(&hnfmod, &t1, &t2);   // hnfmod = norm1^2 * norm2^2

    // // verify hnfmod > 0 (debug, may remove later)
    // assert(ibz_cmp(&hnfmod, &ibz_const_zero) > 0);

    // ibz_finalize(&t1);
    // ibz_finalize(&t2);

    /* paper Issue 11/12 spec: when the validated compact-width routing says
     * that HNF is unsafe, fall through to the spec-faithful "form
     * bar(J_t)·I product, run LLL directly on its 16 column generators" path
     * via ML2(d=16). This matches paper §SuitableIdeals "form ideal product,
     * compute Gram, run LLL" without the canonical-HNF cofactor expansion.
     *
     * The worst-case branch widens ibz_t but deliberately retains the original
     * compact routing threshold: extra storage bits alone do not bound HNF's
     * extended-GCD coefficients and must not silently select a previously
     * untested HNF path. */
    {
        if (use_ml2) {
            int rho = quat_ml2_retry(
                &lat_res.basis[0], 4, generators, 16, alg);
            if (rho == 4) {
                quat_lattice_basis_transpose_in_place(&lat_res.basis);
                ok = 1;
            }
        } else {
            ibz_mat_4xn_hnf_mod_core(&(lat_res.basis), 16, generators, &hnfmod);
            ok = 1;
        }
    }
    if (!ok)
        goto cleanup;
    ibz_mul(&lat_res.denom, &lat1->denom, &lat2->denom);
    if (!quat_lattice_reduce_denom(&lat_res, &lat_res)) {
        ok = 0;
        goto cleanup;
    }
    ibz_mat_4x4_copy(&res->basis, &lat_res.basis);
    ibz_copy(&res->denom, &lat_res.denom);
cleanup:
    ibz_vec_4_finalize(&elem1);
    ibz_vec_4_finalize(&elem2);
    ibz_vec_4_finalize(&elem_res);
    quat_lattice_finalize(&lat_res);
    // ibz_finalize(&det);
    // ibz_mat_4x4_finalize(&(detmat));
    ibz_finalize(&hnfmod);
    ibz_finalize(&r1);
    ibz_finalize(&r2);
    for (int i = 0; i < 16; i++)
        ibz_vec_4_finalize(&(generators[i]));
    return ok;
}

// lattice assumed of full rank
int
quat_lattice_contains(ibz_vec_4_t *coord, const quat_lattice_t *lat, const quat_alg_elem_t *x)
{
    int divisible = 0;
    ibz_vec_4_t work_coord;
    ibz_mat_4x4_t inv;
    ibz_t det, prod;
    ibz_init(&prod);
    ibz_init(&det);
    ibz_vec_4_init(&work_coord);
    ibz_mat_4x4_init(&inv);
    if (!ibz_mat_4x4_inv_with_det_as_denom(&inv, &det, &(lat->basis)))
        goto cleanup;
    ibz_mat_4x4_eval(&work_coord, &inv, &(x->coord));
    ibz_vec_4_scalar_mul(&(work_coord), &(lat->denom), &work_coord);
    ibz_mul(&prod, &(x->denom), &det);
    if (!ibz_is_zero(&prod))
        divisible = ibz_vec_4_scalar_div(&work_coord, &prod, &work_coord);
    // copy result
    if (divisible && (coord != NULL)) {
        for (int i = 0; i < 4; i++) {
            ibz_copy(&((*coord)[i]), &(work_coord[i]));
        }
    }
cleanup:
    ibz_finalize(&prod);
    ibz_finalize(&det);
    ibz_mat_4x4_finalize(&inv);
    ibz_vec_4_finalize(&work_coord);
    return (divisible);
}

int
quat_lattice_index(ibz_t *index, const quat_lattice_t *sublat, const quat_lattice_t *overlat)
{
    int ok = 0;
    ibz_t tmp, det, candidate;
    if (index == NULL || sublat == NULL || overlat == NULL ||
        ibz_is_zero(&sublat->denom) || ibz_is_zero(&overlat->denom))
        return 0;
    ibz_init(&tmp);
    ibz_init(&det);
    ibz_init(&candidate);

    // det = det(sublat->basis)
    if (!ibz_mat_4x4_inv_with_det_as_denom(NULL, &det, &sublat->basis))
        goto cleanup;
    // tmp = (overlat->denom)⁴
    ibz_mul(&tmp, &overlat->denom, &overlat->denom);
    ibz_mul(&tmp, &tmp, &tmp);
    // index = (overlat->denom)⁴ · det(sublat->basis)
    ibz_mul(&candidate, &det, &tmp);
    // tmp = (sublat->denom)⁴
    ibz_mul(&tmp, &sublat->denom, &sublat->denom);
    ibz_mul(&tmp, &tmp, &tmp);
    // det = det(overlat->basis)
    if (!ibz_mat_4x4_inv_with_det_as_denom(NULL, &det, &overlat->basis))
        goto cleanup;
    // tmp = (sublat->denom)⁴ · det(overlat->basis)
    ibz_mul(&tmp, &tmp, &det);
    // index = index / tmp
    if (ibz_is_zero(&tmp))
        goto cleanup;
    ibz_div(&candidate, &tmp, &candidate, &tmp);
    if (!ibz_is_zero(&tmp))
        goto cleanup;
    // index = |index|
    ibz_abs(&candidate, &candidate);
    ibz_copy(index, &candidate);
    ok = 1;

cleanup:
    ibz_finalize(&candidate);
    ibz_finalize(&tmp);
    ibz_finalize(&det);
    return ok;
}

int
quat_lattice_hnf(quat_lattice_t *lat)
{
    ibz_t mod;
    ibz_vec_4_t generators[4];
    quat_lattice_t candidate;
    int ok = 0;
    ibz_init(&mod);
    quat_lattice_init(&candidate);
    ibz_mat_4x4_copy(&candidate.basis, &lat->basis);
    ibz_copy(&candidate.denom, &lat->denom);
    if (!ibz_mat_4x4_inv_with_det_as_denom(NULL, &mod, &candidate.basis))
        goto cleanup;
    ibz_abs(&mod, &mod);
    if (ibz_is_zero(&mod))
        goto cleanup;
    for (int i = 0; i < 4; i++)
        ibz_vec_4_init(&(generators[i]));
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_copy(&(generators[j][i]), &(candidate.basis[i][j]));
        }
    }
    ibz_mat_4xn_hnf_mod_core(&candidate.basis, 4, generators, &mod);
    if (quat_lattice_reduce_denom(&candidate, &candidate)) {
        ibz_mat_4x4_copy(&lat->basis, &candidate.basis);
        ibz_copy(&lat->denom, &candidate.denom);
        ok = 1;
    }
    for (int i = 0; i < 4; i++)
        ibz_vec_4_finalize(&(generators[i]));

cleanup:
    quat_lattice_finalize(&candidate);
    ibz_finalize(&mod);
    return ok;
}

void
quat_lattice_gram(ibz_mat_4x4_t *G, const quat_lattice_t *lattice, const quat_alg_t *alg)
{
    ibz_t tmp;
    ibz_init(&tmp);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j <= i; j++) {
            ibz_set(&(*G)[i][j], 0);
            for (int k = 0; k < 4; k++) {
                ibz_mul(&tmp, &(lattice->basis)[k][i], &(lattice->basis)[k][j]);
                if (k >= 2)
                    ibz_mul(&tmp, &tmp, &alg->p);
                ibz_add(&(*G)[i][j], &(*G)[i][j], &tmp);
            }
            ibz_mul(&(*G)[i][j], &(*G)[i][j], &ibz_const_two);
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int j = i + 1; j < 4; j++) {
            ibz_copy(&(*G)[i][j], &(*G)[j][i]);
        }
    }
    ibz_finalize(&tmp);
}
