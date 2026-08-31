#include <signature.h>
#include <quaternion_constants.h>
#include <quaternion_data.h>
#include <id2iso.h>
#include <torsion_constants.h>

#ifndef SQISIGN_KEYGEN_MAX_ATTEMPTS
#define SQISIGN_KEYGEN_MAX_ATTEMPTS 1024u
#endif

void
secret_key_init(secret_key_t *sk)
{
    quat_left_ideal_init(&(sk->secret_ideal));
    ibz_mat_2x2_init(&(sk->mat_BAcan_to_BA0_two));
    ec_curve_init(&sk->curve);
}

void
secret_key_finalize(secret_key_t *sk)
{
    quat_left_ideal_finalize(&(sk->secret_ideal));
    ibz_mat_2x2_finalize(&(sk->mat_BAcan_to_BA0_two));
}

int
protocols_keygen(public_key_t *pk, secret_key_t *sk)
{
    int found = 0;
    ec_basis_t B_0_two;

    // Probabilistic failures retry, but entropy or arithmetic faults must not
    // leave an embedded device in an unbounded loop.
    for (unsigned attempt = 0;
         !found && attempt < SQISIGN_KEYGEN_MAX_ATTEMPTS;
         ++attempt) {
        quat_random_ideal_status_t sampling_status =
            quat_sampling_random_ideal_O0_given_norm(
            &sk->secret_ideal, &SEC_DEGREE, 1, &QUAT_represent_integer_params, NULL);
        if (sampling_status == QUAT_RANDOM_IDEAL_FATAL)
            return 0;
        found = (sampling_status == QUAT_RANDOM_IDEAL_SUCCESS);

        // replacing the secret key ideal by a shorter equivalent one for efficiency
        found = found && quat_lideal_prime_norm_reduced_equivalent(
                             &sk->secret_ideal, &QUATALG_PINFTY, QUAT_primality_num_iter, QUAT_equiv_bound_coeff);

        // ideal to isogeny clapotis
        if (found) {
            int isogeny_status = dim2id2iso_arbitrary_isogeny_evaluation(
                &B_0_two, &sk->curve, &sk->secret_ideal);
            if (isogeny_status == ID2ISO_STATUS_FATAL)
                return 0;
            found = (isogeny_status > 0);
        }
    }

    if (!found)
        return 0;

    // Assert the isogeny was found and images have the correct order
    assert(test_basis_order_twof(&B_0_two, &sk->curve, TORSION_EVEN_POWER));

    // Compute a deterministic basis with a hint to speed up verification
    if (!ec_curve_to_basis_2f_to_hint(
            &sk->canonical_basis, &sk->curve, TORSION_EVEN_POWER, &pk->hint_pk)) {
        return 0;
    }

    // Assert the deterministic basis we computed has the correct order
    assert(test_basis_order_twof(&sk->canonical_basis, &sk->curve, TORSION_EVEN_POWER));

    // Compute the 2x2 matrix basis change from the canonical basis to the evaluation of our secret
    // isogeny
    change_of_basis_matrix_tate(
        &sk->mat_BAcan_to_BA0_two, &sk->canonical_basis, &B_0_two, &sk->curve, TORSION_EVEN_POWER);

    // Set the public key from the codomain curve
    copy_curve(&pk->curve, &sk->curve);
    pk->curve.is_A24_computed_and_normalized = false; // We don't send any precomputation

    assert(fp2_is_one(&pk->curve.C) == 0xFFFFFFFF);

    return found;
}
