/* SPDX-License-Identifier: Apache-2.0 */
#include "encode_signature.c"
#include "internal.h"
#include <mp.h>
#include "theta_schedule_de2.h"

#include <stdio.h>
#include <stdlib.h>

struct field_vector {
    unsigned char x[FP2_ENCODED_BYTES], y[FP2_ENCODED_BYTES];
    unsigned char sum[FP_ENCODED_BYTES], difference[FP_ENCODED_BYTES], product[FP_ENCODED_BYTES];
    unsigned char inverse[FP_ENCODED_BYTES], square[FP_ENCODED_BYTES];
    unsigned char product2[FP2_ENCODED_BYTES], inverse2[FP2_ENCODED_BYTES], square2[FP2_ENCODED_BYTES];
};

#include "field_vectors.h"

int
PQCLEAN_randombytes(uint8_t *output, size_t length)
{
    (void)output;
    (void)length;
    abort();
}

static void
require(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "encoding regression: %s\n", message);
        exit(EXIT_FAILURE);
    }
}

static void
check_secret_integer(void)
{
    unsigned char encoded[SECRET_IDEAL_BYTES];
    ibz_t input, decoded;
    ibz_mul_2exp(&input, &ibz_const_one, 8 * SECRET_IDEAL_BYTES - 1);
    require(ibz_to_bytes(encoded, &input, sizeof(encoded), true) == NULL,
            "positive signed overflow is rejected");
    ibz_neg(&input, &input);
    require(ibz_to_bytes(encoded, &input, sizeof(encoded), true) == encoded + sizeof(encoded),
            "most negative secret integer is accepted");
    require(encoded[0] == 0 && encoded[sizeof(encoded) - 1] == 0x80, "signed integer byte order");
    require(ibz_from_bytes(&decoded, encoded, sizeof(encoded), true) == encoded + sizeof(encoded) &&
            ibz_cmp(&input, &decoded) == 0, "signed secret integer sign extension");
    ibz_sub(&input, &input, &ibz_const_one);
    require(ibz_to_bytes(encoded, &input, sizeof(encoded), true) == NULL,
            "negative signed overflow is rejected");
}

static void
check_response_encoding(void)
{
    signature_t sig = {0}, decoded;
    unsigned char bytes[SIGNATURE_BYTES];
    static const uint16_t test_values[] = {
        0,
#if RESPONSE_EXPONENT_BYTES == 2
        255, 256,
#endif
        SQIsign_response_length
    };
    for (size_t i = 0; i < sizeof(test_values) / sizeof(*test_values); ++i) {
        sig.backtracking = test_values[i];
        sig.two_resp_length = SQIsign_response_length - sig.backtracking;
        signature_to_bytes(bytes, &sig);
        require(signature_from_bytes(&decoded, bytes) && decoded.backtracking == sig.backtracking &&
                decoded.two_resp_length == sig.two_resp_length, "response exponent round trip");
        require(bytes[FP2_ENCODED_BYTES] == (unsigned char)sig.backtracking, "low exponent byte");
#if RESPONSE_EXPONENT_BYTES == 2
        require(bytes[FP2_ENCODED_BYTES + 1] == (sig.backtracking >> 8), "high exponent byte");
#endif
    }
    bytes[FP2_ENCODED_BYTES + RESPONSE_EXPONENT_BYTES] = 1;
    require(!signature_from_bytes(&decoded, bytes), "out-of-range exponent sum is rejected");
    memset(bytes, 0xff, FP2_ENCODED_BYTES);
    require(!signature_from_bytes(&decoded, bytes), "noncanonical field encoding is rejected");
}

static void
check_field_boundaries(void)
{
    unsigned char encoded[FP_ENCODED_BYTES], round_trip[FP_ENCODED_BYTES];
    fp_t value, one;
    ibz_t below_prime;
    fp_set_one(&one);
    require(ibz_to_bytes(encoded, &QUATALG_PINFTY.p, sizeof(encoded), false) != NULL &&
            !fp_decode(&value, encoded), "the modulus is not canonical");
    ibz_sub(&below_prime, &QUATALG_PINFTY.p, &ibz_const_one);
    require(ibz_to_bytes(encoded, &below_prime, sizeof(encoded), false) != NULL &&
            fp_decode(&value, encoded), "p minus one decodes");
    fp_encode(round_trip, &value);
    require(memcmp(encoded, round_trip, sizeof(encoded)) == 0, "p minus one round trip");
    fp_mul(&value, &value, &value);
    require(fp_is_equal(&value, &one), "in-place multiplication of p minus one");
    memset(encoded, 0, sizeof(encoded));
    require(fp_decode(&value, encoded) && fp_is_zero(&value), "zero decodes");
    encoded[0] = 1;
    require(fp_decode(&value, encoded) && fp_is_equal(&value, &one), "one decodes");
}

static void
check_scalar_inverse(void)
{
    scalar_t a = {1}, b = {2}, c = {3}, d = {5}, expected = {0}, zero = {0};
    const scalar_t one = {1}, five = {5};
    mp_invert_matrix(a, b, c, d, TORSION_EVEN_POWER, NWORDS_ORDER);
    mp_mod_2exp(a, TORSION_EVEN_POWER, NWORDS_ORDER);
    mp_mod_2exp(b, TORSION_EVEN_POWER, NWORDS_ORDER);
    mp_mod_2exp(c, TORSION_EVEN_POWER, NWORDS_ORDER);
    mp_mod_2exp(d, TORSION_EVEN_POWER, NWORDS_ORDER);
    mp_sub(expected, zero, five, NWORDS_ORDER);
    mp_mod_2exp(expected, TORSION_EVEN_POWER, NWORDS_ORDER);
    require(mp_compare(a, expected, NWORDS_ORDER) == 0 && b[0] == 2 && c[0] == 3,
            "full-width scalar inverse first entries");
    mp_sub(expected, zero, one, NWORDS_ORDER);
    mp_mod_2exp(expected, TORSION_EVEN_POWER, NWORDS_ORDER);
    require(mp_compare(d, expected, NWORDS_ORDER) == 0, "full-width scalar inverse last entry");
}

static void
check_field_oracle(void)
{
    unsigned char encoded[FP2_ENCODED_BYTES];
    for (size_t i = 0; i < sizeof(field_vectors) / sizeof(field_vectors[0]); ++i) {
        const struct field_vector *test = &field_vectors[i];
        fp_t x, y, result;
        require(fp_decode(&x, test->x) && fp_decode(&y, test->y), "oracle inputs in Fp");
        fp_add(&result, &x, &y);
        fp_encode(encoded, &result);
        require(memcmp(encoded, test->sum, FP_ENCODED_BYTES) == 0, "Fp sum against oracle");
        fp_sub(&result, &x, &y);
        fp_encode(encoded, &result);
        require(memcmp(encoded, test->difference, FP_ENCODED_BYTES) == 0, "Fp difference against oracle");
        fp_mul(&y, &x, &y);
        fp_encode(encoded, &y);
        require(memcmp(encoded, test->product, FP_ENCODED_BYTES) == 0, "aliased Fp product against oracle");
        fp_copy(&result, &x);
        fp_inv(&result);
        fp_encode(encoded, &result);
        require(memcmp(encoded, test->inverse, FP_ENCODED_BYTES) == 0, "Fp inverse against oracle");
        fp_sqr(&result, &x);
        fp_encode(encoded, &result);
        require(memcmp(encoded, test->square, FP_ENCODED_BYTES) == 0, "Fp square against oracle");
        fp_sqrt(&result);
        fp_sqr(&result, &result);
        fp_encode(encoded, &result);
        require(memcmp(encoded, test->square, FP_ENCODED_BYTES) == 0, "Fp square root against oracle");

        fp2_t x2, y2, result2;
        require(fp2_decode(&x2, test->x) && fp2_decode(&y2, test->y), "oracle inputs in Fp2");
        fp2_mul(&y2, &x2, &y2);
        fp2_encode(encoded, &y2);
        require(memcmp(encoded, test->product2, sizeof(encoded)) == 0, "aliased Fp2 product against oracle");
        fp2_copy(&result2, &x2);
        fp2_inv(&result2);
        fp2_encode(encoded, &result2);
        require(memcmp(encoded, test->inverse2, sizeof(encoded)) == 0, "Fp2 inverse against oracle");
        fp2_sqr(&result2, &x2);
        fp2_encode(encoded, &result2);
        require(memcmp(encoded, test->square2, sizeof(encoded)) == 0, "Fp2 square against oracle");
        require(fp2_sqrt_verify(&result2), "Fp2 oracle square has a square root");
        fp2_sqr(&result2, &result2);
        fp2_encode(encoded, &result2);
        require(memcmp(encoded, test->square2, sizeof(encoded)) == 0, "Fp2 square root against oracle");
    }
}

static void
create_scalar_ideal(quat_left_ideal_t *ideal, int scalar)
{
    quat_alg_elem_t generator;
    quat_alg_elem_init(&generator);
    quat_alg_elem_set(&generator, 1, scalar, 0, 0, 0);
    require(quat_lideal_create_principal(ideal, &generator, &MAXORD_O0, &QUATALG_PINFTY),
            "scalar principal ideal");
    quat_alg_elem_finalize(&generator);
}

static void
check_compact_ideals(void)
{
    quat_left_ideal_t left, right, expected;
    quat_left_ideal_init(&left);
    quat_left_ideal_init(&right);
    quat_left_ideal_init(&expected);

    create_scalar_ideal(&left, 6);
    create_scalar_ideal(&right, 10);
    create_scalar_ideal(&expected, 30);
    require(quat_lideal_inter(&left, &left, &right, &QUATALG_PINFTY), "aliased compact intersection");
    require(ibz_cmp_int32(&left.norm, 900) == 0 && quat_lattice_equal(&left.lattice, &expected.lattice),
            "trace-gcd intersection is lcm(6,10) times the order");
    ibz_t index;
    require(quat_lattice_index(&index, &left.lattice, &MAXORD_O0) && ibz_cmp_int32(&index, 810000) == 0,
            "intersection norm and actual index agree");

    create_scalar_ideal(&left, 3);
    create_scalar_ideal(&right, 5);
    create_scalar_ideal(&expected, 15);
    require(quat_lattice_mul_mlll(&right.lattice, &left.lattice, &right.lattice, &QUATALG_PINFTY) &&
            quat_lattice_equal(&right.lattice, &expected.lattice), "aliased coprime ordered product");

    quat_left_ideal_finalize(&left);
    quat_left_ideal_finalize(&right);
    quat_left_ideal_finalize(&expected);
}

static void
check_theta_schedule(void)
{
    for (unsigned extra = 0; extra <= 2; extra += 2) {
        for (unsigned length = extra ? 1 : 4; length + extra <= TORSION_EVEN_POWER; ++length) {
            uint16_t todo[THETA_SCHEDULE_DE2_MAX_SLOTS] = {length - 2 + extra};
            int current = 0;
            while (todo[current] != 1) {
                const unsigned remaining = todo[current];
                require(current + 1 < 8, "gluing schedule exceeds the Jacobian slots");
                const unsigned split = remaining >= 16 ? remaining / 2 : remaining - 1;
                todo[++current] = remaining - split;
            }
            unsigned leaves = 0;
            while (current >= 0) {
                while (todo[current] != 1) {
                    uint16_t split;
                    require(current + 1 < THETA_SCHEDULE_DE2_MAX_SLOTS &&
                            theta_schedule_de2_lookup((unsigned)current, todo[current], &split),
                            "ordinary schedule does not fit the allocated slots");
                    todo[current + 1] = todo[current] - split;
                    ++current;
                }
                ++leaves;
                for (int i = 0; i < current; ++i) {
                    require(todo[i] > 1, "schedule discarded a live ancestor");
                    --todo[i];
                }
                --current;
            }
            require(leaves == length - 2 + extra, "schedule omitted an isogeny step");
        }
    }
}

int
main(void)
{
    check_secret_integer();
    check_response_encoding();
    check_field_boundaries();
    check_field_oracle();
    check_compact_ideals();
    check_scalar_inverse();
    check_theta_schedule();
    puts("encoding, field and ideal oracles, scalar inverse and theta schedule: PASS");
    return EXIT_SUCCESS;
}
