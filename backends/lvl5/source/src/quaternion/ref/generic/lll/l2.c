#include <quaternion.h>
#include "lll_internals.h"
#include "internal.h"

#include "dpe.h"

// Access entry of symmetric matrix
#define SYM(M, i, j) (i < j ? &M[j][i] : &M[i][j])

static int
lll_mul_checked(ibz_t *product, const ibz_t *a, const ibz_t *b)
{
    int ok = 0;
    ibz_t abs_a, abs_b, limit, quotient, remainder;
    ibz_init(&abs_a); ibz_init(&abs_b); ibz_init(&limit);
    ibz_init(&quotient); ibz_init(&remainder);
    ibz_abs(&abs_a, a);
    ibz_abs(&abs_b, b);
    if (ibz_is_negative(&abs_a) || ibz_is_negative(&abs_b))
        goto cleanup;
    if (ibz_is_zero(&abs_a) || ibz_is_zero(&abs_b)) {
        ibz_set(product, 0);
        ok = 1;
        goto cleanup;
    }
    for (int i = 0; i < IBZ_LIMBS; i++)
        limit[i] = UINT64_MAX;
    limit[IBZ_LIMBS - 1] >>= 1;
    ibz_div(&quotient, &remainder, &limit, &abs_b);
    if (ibz_cmp(&abs_a, &quotient) > 0)
        goto cleanup;
    ibz_mul(product, a, b);
    ok = 1;
cleanup:
    ibz_finalize(&abs_a); ibz_finalize(&abs_b); ibz_finalize(&limit);
    ibz_finalize(&quotient); ibz_finalize(&remainder);
    return ok;
}

/* Compute a-b in a temporary two's-complement array and inspect the signs
 * before publishing it.  This remains a status path even when the global
 * SQISIGN_INTBIG_OVERFLOW_CHECK build would otherwise abort. */
static int
lll_sub_checked(ibz_t *difference, const ibz_t *a, const ibz_t *b)
{
    ibz_t work;
    uint64_t borrow = 0;
    int a_negative = ibz_is_negative(a);
    int b_negative = ibz_is_negative(b);
    ibz_init(&work);
    for (int i = 0; i < IBZ_LIMBS; i++) {
        uint64_t first = (*a)[i] - borrow;
        borrow = first > (*a)[i];
        uint64_t second = first - (*b)[i];
        borrow += second > first;
        work[i] = second;
    }
    if (a_negative != b_negative && ibz_is_negative(&work) != a_negative) {
        ibz_finalize(&work);
        return 0;
    }
    ibz_copy(difference, &work);
    ibz_finalize(&work);
    return 1;
}

static int
lll_gram_fits(const quat_lattice_t *lattice, const quat_alg_t *alg)
{
    const int p_bits = ibz_bitsize(&alg->p);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            int max_product_bits = 0;
            for (int k = 0; k < 4; k++) {
                if (ibz_is_zero(&lattice->basis[k][i]) ||
                    ibz_is_zero(&lattice->basis[k][j]))
                    continue;
                int product_bits = ibz_bitsize(&lattice->basis[k][i]) +
                                   ibz_bitsize(&lattice->basis[k][j]);
                if (k >= 2)
                    product_bits += p_bits;
                if (product_bits > max_product_bits)
                    max_product_bits = product_bits;
            }
            /* Four terms and the final factor two require at most three
             * additional magnitude bits. */
            if (max_product_bits > IBZ_BITS - 4)
                return 0;
        }
    }
    return 1;
}

int
quat_lll_core_checked(ibz_mat_4x4_t *G, ibz_mat_4x4_t *basis)
{
    int status = 0;
    uint64_t iterations = 0;
    if (G == NULL || basis == NULL || ibz_cmp(&(*G)[0][0], &ibz_const_zero) <= 0)
        return 0;
    dpe_t dpe_const_one, dpe_const_DELTABAR;

    dpe_init(dpe_const_one);
    dpe_set_ui(dpe_const_one, 1);

    dpe_init(dpe_const_DELTABAR);
    dpe_set_d(dpe_const_DELTABAR, DELTABAR);

    // fp variables for Gram-Schmidt orthogonalization and Lovasz' conditions
    dpe_t r[4][4], u[4][4], lovasz[4];
    for (int i = 0; i < 4; i++) {
        dpe_init(lovasz[i]);
        for (int j = 0; j <= i; j++) {
            dpe_init(r[i][j]);
            dpe_init(u[i][j]);
        }
    }

    // threshold for swaps
    dpe_t delta_bar;
    dpe_init(delta_bar);
    dpe_set_d(delta_bar, DELTABAR);

    // Other work variables
    dpe_t Xf, tmpF;
    dpe_init(Xf);
    dpe_init(tmpF);
    ibz_t X, tmpI;
    ibz_init(&X);
    ibz_init(&tmpI);

    // Main L² loop
    // dpe_set_z(r[0][0], (*G)[0][0]);
    dpe_set_z(r[0][0], &(*G)[0][0]); 

    int kappa = 1;
    while (kappa < 4) {
        if (++iterations > UINT64_C(1000000))
            goto cleanup;
        // size reduce b_κ
        int done = 0;
        while (!done) {
            // Recompute the κ-th row of the Choleski Factorisation
            // Loop invariant:
            //     r[κ][j] ≈ u[κ][j] ‖b_j*‖² ≈ 〈b_κ, b_j*〉
            for (int j = 0; j <= kappa; j++) {
                // dpe_set_z(r[kappa][j], (*G)[kappa][j]);
                dpe_set_z(r[kappa][j], &(*G)[kappa][j]);
                for (int k = 0; k < j; k++) {
                    dpe_mul(tmpF, r[kappa][k], u[j][k]);
                    dpe_sub(r[kappa][j], r[kappa][j], tmpF);
                }
                if (j < kappa) {
                    if (dpe_cmp_d(r[j][j], 0.0) <= 0)
                        goto cleanup;
                    dpe_div(u[kappa][j], r[kappa][j], r[j][j]);
                }
            }

            done = 1;
            // size reduce
            for (int i = kappa - 1; i >= 0; i--) {
                if (dpe_cmp_d(u[kappa][i], ETABAR) > 0 || dpe_cmp_d(u[kappa][i], -ETABAR) < 0) {
                    done = 0;
                    dpe_set(Xf, u[kappa][i]);
                    dpe_round(Xf, Xf);
                    if (DPE_EXP(Xf) >= IBZ_BITS - 1)
                        goto cleanup;
                    dpe_get_z(&X, Xf);
                    if (ibz_is_zero(&X))
                        goto cleanup;

                    // Update basis: b_κ ← b_κ - X·b_i
                    for (int j = 0; j < 4; j++) {
                        if (!lll_mul_checked(&tmpI, &X, &(*basis)[j][i]) ||
                            !lll_sub_checked(&(*basis)[j][kappa], &(*basis)[j][kappa], &tmpI))
                            goto cleanup;
                    }
                    // Update lower half of the Gram matrix
                    // <b_κ - X·b_i, b_κ - X·b_i> = <b_κ, b_κ> - 2X<b_κ, b_i> + X²<b_i, b_i> =
                    // <b_κ,b_κ> - X<b_κ,b_i> - X(<b_κ,b_i> - X·<b_i, b_i>)
                    //// 〈b_κ, b_κ〉 ← 〈b_κ, b_κ〉 - X·〈b_κ, b_i〉
                    if (!lll_mul_checked(&tmpI, &X, &(*G)[kappa][i]) ||
                        !lll_sub_checked(&(*G)[kappa][kappa], &(*G)[kappa][kappa], &tmpI))
                        goto cleanup;
                    for (int j = 0; j < 4; j++) { // works because i < κ
                        // 〈b_κ, b_j〉 ← 〈b_κ, b_j〉 - X·〈b_i, b_j〉
                        if (!lll_mul_checked(&tmpI, &X, SYM((*G), i, j)) ||
                            !lll_sub_checked(SYM((*G), kappa, j), SYM((*G), kappa, j), &tmpI))
                            goto cleanup;
                    }
                    // After the loop:
                    //// 〈b_κ,b_κ〉 ← 〈b_κ,b_κ〉 - X·〈b_κ,b_i〉 - X·(〈b_κ,b_i〉 - X·〈b_i,
                    /// b_i〉) = 〈b_κ - X·b_i, b_κ - X·b_i〉
                    //
                    // Update u[kappa][j]
                    for (int j = 0; j < i; j++) {
                        dpe_mul(tmpF, Xf, u[i][j]);
                        dpe_sub(u[kappa][j], u[kappa][j], tmpF);
                    }
                }
            }
        }

        // Check Lovasz' conditions
        // lovasz[0] = ‖b_κ‖²
        // dpe_set_z(lovasz[0], (*G)[kappa][kappa]);
        dpe_set_z(lovasz[0], &(*G)[kappa][kappa]);
        // lovasz[i] = lovasz[i-1] - u[κ][i-1]·r[κ][i-1]
        for (int i = 1; i < kappa; i++) {
            dpe_mul(tmpF, u[kappa][i - 1], r[kappa][i - 1]);
            dpe_sub(lovasz[i], lovasz[i - 1], tmpF);
        }
        int swap;
        for (swap = kappa; swap > 0; swap--) {
            dpe_mul(tmpF, delta_bar, r[swap - 1][swap - 1]);
            if (dpe_cmp(tmpF, lovasz[swap - 1]) < 0)
                break;
        }

        // Insert b_κ before b_swap
        if (kappa != swap) {
            // Insert b_κ before b_swap in the basis and in the lower half Gram matrix
            for (int j = kappa; j > swap; j--) {
                for (int i = 0; i < 4; i++) {
                    ibz_swap(&(*basis)[i][j], &(*basis)[i][j - 1]);
                    if (i == j - 1)
                        ibz_swap(&(*G)[i][i], &(*G)[j][j]);
                    else if (i != j)
                        ibz_swap(SYM((*G), i, j), SYM((*G), i, j - 1));
                }
            }
            // Copy row u[κ] and r[κ] in swap position, ignore what follows
            for (int i = 0; i < swap; i++) {
                dpe_set(u[swap][i], u[kappa][i]);
                dpe_set(r[swap][i], r[kappa][i]);
            }
            dpe_set(r[swap][swap], lovasz[swap]);
            // swap complete
            kappa = swap;
        }

        kappa += 1;
    }

#ifndef NDEBUG
    // Check size-reducedness
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < i; j++) {
            dpe_abs(u[i][j], u[i][j]);
            assert(dpe_cmp_d(u[i][j], ETABAR) <= 0);
        }
    // Check Lovasz' conditions
    for (int i = 1; i < 4; i++) {
        dpe_mul(tmpF, u[i][i - 1], u[i][i - 1]);
        dpe_sub(tmpF, dpe_const_DELTABAR, tmpF);
        dpe_mul(tmpF, tmpF, r[i - 1][i - 1]);
        assert(dpe_cmp(tmpF, r[i][i]) <= 0);
    }
#endif

    // Fill in the upper half of the Gram matrix
    for (int i = 0; i < 4; i++) {
        for (int j = i + 1; j < 4; j++)
            ibz_copy(&(*G)[i][j], &(*G)[j][i]);
    }

    status = 1;

    // Clearinghouse
cleanup:
    ibz_finalize(&X);
    ibz_finalize(&tmpI);
    dpe_clear(dpe_const_one);
    dpe_clear(dpe_const_DELTABAR);
    dpe_clear(Xf);
    dpe_clear(tmpF);
    dpe_clear(delta_bar);
    for (int i = 0; i < 4; i++) {
        dpe_clear(lovasz[i]);
        for (int j = 0; j <= i; j++) {
            dpe_clear(r[i][j]);
            dpe_clear(u[i][j]);
        }
    }
    return status;
}

void
quat_lll_core(ibz_mat_4x4_t *G, ibz_mat_4x4_t *basis)
{
    (void)quat_lll_core_checked(G, basis);
}

int
quat_lattice_lll(ibz_mat_4x4_t *red, const quat_lattice_t *lattice, const quat_alg_t *alg)
{
    if (red == NULL || lattice == NULL || alg == NULL ||
        ibz_is_zero(&lattice->denom) ||
        ibz_cmp(&alg->p, &ibz_const_zero) <= 0 ||
        !lll_gram_fits(lattice, alg))
        return 0;

    ibz_mat_4x4_t G; // Gram Matrix
    ibz_mat_4x4_t candidate;
    ibz_mat_4x4_init(&G);
    ibz_mat_4x4_init(&candidate);
    quat_lattice_gram(&G, lattice, alg);
    if (ibz_cmp(&G[0][0], &ibz_const_zero) <= 0) {
        ibz_mat_4x4_finalize(&candidate);
        ibz_mat_4x4_finalize(&G);
        return 0;
    }
    ibz_mat_4x4_copy(&candidate, &lattice->basis);
    if (!quat_lll_core_checked(&G, &candidate)) {
        ibz_mat_4x4_finalize(&candidate);
        ibz_mat_4x4_finalize(&G);
        return 0;
    }
    ibz_mat_4x4_copy(red, &candidate);
    ibz_mat_4x4_finalize(&candidate);
    ibz_mat_4x4_finalize(&G);
    return 1;
}
