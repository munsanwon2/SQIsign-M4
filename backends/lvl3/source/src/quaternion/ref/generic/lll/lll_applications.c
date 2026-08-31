#include <quaternion.h>
#include <internal.h>
#include <stdint.h>
#include "lll_internals.h"

int
quat_lideal_reduce_basis(ibz_mat_4x4_t *reduced,
                         ibz_mat_4x4_t *gram,
                         const quat_left_ideal_t *lideal,
                         const quat_alg_t *alg)
{
    if (!quat_order_is_maximal((lideal->parent_order), alg))
        return 0;
    ibz_t gram_corrector;
    ibz_init(&gram_corrector);
    ibz_mul(&gram_corrector, &(lideal->lattice.denom), &(lideal->lattice.denom));
    if (!quat_lideal_class_gram(gram, lideal, alg)) {
        ibz_finalize(&gram_corrector);
        return 0;
    }
    ibz_mat_4x4_copy(reduced, &(lideal->lattice.basis));
    if (!quat_lll_core_checked(gram, reduced)) {
        ibz_finalize(&gram_corrector);
        return 0;
    }
    ibz_mat_4x4_scalar_mul(gram, &gram_corrector, gram);
    for (int i = 0; i < 4; i++) {
        ibz_div_2exp(&((*gram)[i][i]), &((*gram)[i][i]), 1);
        for (int j = i + 1; j < 4; j++) {
            ibz_set(&((*gram)[i][j]), 0);
        }
    }
    ibz_finalize(&gram_corrector);
    return 1;
}

int
quat_lideal_lideal_mul_reduced(quat_left_ideal_t *prod,
                               ibz_mat_4x4_t *gram,
                               const quat_left_ideal_t *lideal1,
                               const quat_left_ideal_t *lideal2,
                               const quat_alg_t *alg)
{
    if (prod == NULL || gram == NULL || lideal1 == NULL || lideal2 == NULL ||
        alg == NULL)
        return 0;

    int ok = 0;
    ibz_mat_4x4_t red, candidate_gram;
    quat_left_ideal_t candidate_prod;
    ibz_mat_4x4_init(&red);
    ibz_mat_4x4_init(&candidate_gram);
    quat_left_ideal_init(&candidate_prod);

    if (!quat_lattice_mul(&candidate_prod.lattice,
                          &(lideal1->lattice),
                          &(lideal2->lattice),
                          alg,
                          (&lideal1->norm),
                          &(lideal2->norm))) {
        goto cleanup;
    }

    candidate_prod.parent_order = lideal1->parent_order;
    if (!quat_lideal_norm(&candidate_prod) ||
        !quat_lideal_reduce_basis(&red, &candidate_gram, &candidate_prod, alg))
        goto cleanup;
    ibz_mat_4x4_copy(&candidate_prod.lattice.basis, &red);

    quat_lideal_copy(prod, &candidate_prod);
    ibz_mat_4x4_copy(gram, &candidate_gram);
    ok = 1;

cleanup:
    quat_left_ideal_finalize(&candidate_prod);
    ibz_mat_4x4_finalize(&candidate_gram);
    ibz_mat_4x4_finalize(&red);
    return ok;
}

int
quat_lideal_prime_norm_reduced_equivalent(quat_left_ideal_t *lideal,
                                          const quat_alg_t *alg,
                                          const int primality_num_iter,
                                          const int equiv_bound_coeff)
{
    ibz_mat_4x4_t gram, red;
    ibz_mat_4x4_init(&gram);
    ibz_mat_4x4_init(&red);

    int found = 0;

    if (primality_num_iter <= 0 || equiv_bound_coeff < 0)
        goto cleanup_matrices;

    // computing the reduced basis
    if (!quat_lideal_reduce_basis(&red, &gram, lideal, alg))
        goto cleanup_matrices;

    quat_alg_elem_t new_alpha;
    quat_left_ideal_t candidate;
    quat_alg_elem_init(&new_alpha);
    quat_left_ideal_init(&candidate);
    ibz_t tmp, remainder, adjusted_norm;
    ibz_init(&tmp);
    ibz_init(&remainder);
    ibz_init(&adjusted_norm);

    ibz_mul(&adjusted_norm, &lideal->lattice.denom, &lideal->lattice.denom);

    uint64_t ctr = 0;

    // equiv_num_iter = (2 * equiv_bound_coeff + 1)^4
    uint64_t side = 2 * (uint64_t)equiv_bound_coeff + 1;
    if (side > UINT16_MAX)
        goto cleanup;
    uint64_t equiv_num_iter = side * side;
    equiv_num_iter *= equiv_num_iter;

    while (!found && ctr < equiv_num_iter) {
        ctr++;
        // we select our linear combination at random
        if (!ibz_rand_interval_minm_m(&new_alpha.coord[0], equiv_bound_coeff) ||
            !ibz_rand_interval_minm_m(&new_alpha.coord[1], equiv_bound_coeff) ||
            !ibz_rand_interval_minm_m(&new_alpha.coord[2], equiv_bound_coeff) ||
            !ibz_rand_interval_minm_m(&new_alpha.coord[3], equiv_bound_coeff))
            goto cleanup;

        // computation of the norm of the vector sampled
        quat_qf_eval(&tmp, &gram, &new_alpha.coord);

        // compute the norm of the equivalent ideal
        // can be improved by removing the power of two first and the odd part only if the trial
        // division failed (this should always be called on an ideal of norm 2^x * N for some
        // big prime N )
        if (ibz_is_zero(&adjusted_norm))
            goto cleanup;
        ibz_div(&tmp, &remainder, &tmp, &adjusted_norm);
        if (!ibz_is_zero(&remainder))
            goto cleanup;

        // pseudo-primality test
        if (ibz_probab_prime(&tmp, primality_num_iter)) {

            // computes the generator using a matrix multiplication
            ibz_mat_4x4_eval(&new_alpha.coord, &red, &new_alpha.coord);
            ibz_copy(&new_alpha.denom, &lideal->lattice.denom);
            if (!quat_lattice_contains(NULL, &lideal->lattice, &new_alpha))
                goto cleanup;

            quat_alg_conj(&new_alpha, &new_alpha);
            ibz_mul(&new_alpha.denom, &new_alpha.denom, &lideal->norm);
            if (!quat_lideal_mul(&candidate, lideal, &new_alpha, alg))
                goto cleanup;
            /* tmp was accepted as prime above and becomes the norm of the
             * equivalent ideal.  Re-running the randomized primality test in
             * an assert consumed the application DRBG only in Debug builds,
             * making KAT output depend on NDEBUG. */
            if (ibz_cmp(&candidate.norm, &tmp) != 0)
                goto cleanup;

            quat_lideal_copy(lideal, &candidate);

            found = 1;
            break;
        }
    }
cleanup:
    ibz_finalize(&tmp);
    ibz_finalize(&remainder);
    ibz_finalize(&adjusted_norm);
    quat_left_ideal_finalize(&candidate);
    quat_alg_elem_finalize(&new_alpha);

cleanup_matrices:
    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&red);

    return found;
}
