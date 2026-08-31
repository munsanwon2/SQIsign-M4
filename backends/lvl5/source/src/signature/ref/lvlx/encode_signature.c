#include <signature.h>
#include <string.h>
#include <tutil.h>
#include <fp2.h>
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <encoded_sizes.h>
#include <assert.h>

typedef unsigned char byte_t;

// digits

static void
encode_digits(byte_t *enc, const ibz_digit_t *x, size_t nbytes)
{
#ifdef TARGET_BIG_ENDIAN
    const size_t ndigits = nbytes / sizeof(ibz_digit_t);
    const size_t rem = nbytes % sizeof(ibz_digit_t);

    for (size_t i = 0; i < ndigits; i++) {
        ibz_digit_t word = BSWAP64(x[i]);
        memcpy(enc + i * sizeof(ibz_digit_t), &word, sizeof(word));
    }
    if (rem) {
        ibz_digit_t ld = BSWAP64(x[ndigits]);
        memcpy(enc + ndigits * sizeof(ibz_digit_t), (byte_t *)&ld, rem);
    }
#else
    memcpy(enc, (const byte_t *)x, nbytes);
#endif
}

static void
decode_digits(ibz_digit_t *x, const byte_t *enc, size_t nbytes, size_t ndigits)
{
    assert(nbytes <= ndigits * sizeof(ibz_digit_t));
    memcpy((byte_t *)x, enc, nbytes);
    memset((byte_t *)x + nbytes, 0, ndigits * sizeof(ibz_digit_t) - nbytes);

#ifdef TARGET_BIG_ENDIAN
    for (size_t i = 0; i < ndigits; i++)
        x[i] = BSWAP64(x[i]);
#endif
}

// ibz_t

static byte_t *
ibz_to_bytes(byte_t *enc, const ibz_t *x, size_t nbytes, bool sgn)
{
    if (nbytes == 0 || nbytes > IBZ_LIMBS * sizeof(ibz_digit_t) || nbytes > SIZE_MAX / 8) {
        return NULL;
    }

    const size_t encoded_bits = 8 * nbytes;
    const int negative = ibz_is_negative(x);
    const int magnitude_bits = ibz_bitsize(x);

    if ((!sgn && negative) ||
        (!negative && (size_t)magnitude_bits > encoded_bits - (sgn ? 1 : 0))) {
        return NULL;
    }

    if (negative && encoded_bits < IBZ_BITS && (size_t)magnitude_bits >= encoded_bits) {
        // A signed n-byte integer may have magnitude 2^(8n-1) only for its
        // most-negative value. No larger magnitude is representable.
        ibz_t magnitude, limit;
        ibz_init(&magnitude);
        ibz_init(&limit);
        ibz_abs(&magnitude, x);
        ibz_pow(&limit, &ibz_const_two, (uint32_t)(encoded_bits - 1));
        int fits = ibz_cmp(&magnitude, &limit) == 0;
        ibz_finalize(&limit);
        ibz_finalize(&magnitude);
        if (!fits) {
            return NULL;
        }
    }

    // Use full intbig capacity so conversion cannot overrun a field-sized VLA
    // before the encoded-size check has a chance to reject the value.
    ibz_digit_t d[IBZ_LIMBS];
    memset(d, 0, sizeof(d));
    if (ibz_cmp(x, &ibz_const_zero) >= 0) {
        // non-negative, straightforward.
        if (!ibz_to_u64_digits_checked(d, IBZ_LIMBS, x)) {
            return NULL;
        }
    } else {
        // negative; use two's complement.
        ibz_t tmp;
        ibz_init(&tmp);
        ibz_neg(&tmp, x);
        ibz_sub(&tmp, &tmp, &ibz_const_one);
        if (!ibz_to_u64_digits_checked(d, IBZ_LIMBS, &tmp)) {
            ibz_finalize(&tmp);
            return NULL;
        }
        for (size_t i = 0; i < IBZ_LIMBS; ++i)
            d[i] = ~d[i];
        ibz_finalize(&tmp);
    }
    encode_digits(enc, d, nbytes);
    return enc + nbytes;
}

static const byte_t *
ibz_from_bytes(ibz_t *x, const byte_t *enc, size_t nbytes, bool sgn)
{
    assert(nbytes > 0);
    const size_t ndigits = (nbytes + sizeof(ibz_digit_t) - 1) / sizeof(ibz_digit_t);
    assert(ndigits > 0);
    if (ndigits > IBZ_LIMBS)
        return NULL;
    ibz_digit_t d[IBZ_LIMBS];
    memset(d, 0, sizeof(d));
    decode_digits(d, enc, nbytes, ndigits);
    if (sgn && enc[nbytes - 1] >> 7) {
        // negative, decode two's complement
        const size_t rem = nbytes % sizeof(ibz_digit_t);
        if (rem != 0) {
            d[ndigits - 1] |= (ibz_digit_t)-1 << (8 * rem);
        }
        for (size_t i = 0; i < ndigits; ++i)
            d[i] = ~d[i];
        ibz_copy_u64_digits(x, d, ndigits);
        ibz_add(x, x, &ibz_const_one);
        ibz_neg(x, x);
    } else {
        // non-negative
        ibz_copy_u64_digits(x, d, ndigits);
    }
    return enc + nbytes;
}

// public API

int
secret_key_to_bytes(byte_t *enc, const secret_key_t *sk, const public_key_t *pk)
{
    byte_t *const start = enc;

    enc = public_key_to_bytes(enc, pk);

    {
        fp2_t lhs, rhs;
        fp2_mul(&lhs, &sk->curve.A, &pk->curve.C);
        fp2_mul(&rhs, &sk->curve.C, &pk->curve.A);
        if (!fp2_is_equal(&lhs, &rhs)) {
            goto fail;
        }
    }

    enc = ibz_to_bytes(enc, &sk->secret_ideal.norm, FP_ENCODED_BYTES, false);
    if (enc == NULL) {
        goto fail;
    }
    {
        quat_alg_elem_t gen;
        quat_alg_elem_init(&gen);
        if (!quat_lideal_generator(&gen, &sk->secret_ideal, &QUATALG_PINFTY)) {
            quat_alg_elem_finalize(&gen);
            goto fail;
        }
        // we skip encoding the denominator since it won't change the generated ideal
        {
            // let's make sure that the denominator is indeed coprime to the norm of the ideal
            ibz_t gcd;
            ibz_init(&gcd);
            ibz_gcd(&gcd, &gen.denom, &sk->secret_ideal.norm);
            int coprime = ibz_is_one(&gcd);
            ibz_finalize(&gcd);
            if (!coprime) {
                quat_alg_elem_finalize(&gen);
                goto fail;
            }
        }
        enc = ibz_to_bytes(enc, &gen.coord[0], FP_ENCODED_BYTES, true);
        if (enc == NULL) {
            quat_alg_elem_finalize(&gen);
            goto fail;
        }
        enc = ibz_to_bytes(enc, &gen.coord[1], FP_ENCODED_BYTES, true);
        if (enc == NULL) {
            quat_alg_elem_finalize(&gen);
            goto fail;
        }
        enc = ibz_to_bytes(enc, &gen.coord[2], FP_ENCODED_BYTES, true);
        if (enc == NULL) {
            quat_alg_elem_finalize(&gen);
            goto fail;
        }
        enc = ibz_to_bytes(enc, &gen.coord[3], FP_ENCODED_BYTES, true);
        if (enc == NULL) {
            quat_alg_elem_finalize(&gen);
            goto fail;
        }
        quat_alg_elem_finalize(&gen);
    }

    for (int row = 0; row < 2; row++) {
        for (int col = 0; col < 2; col++) {
            if (ibz_is_negative(&sk->mat_BAcan_to_BA0_two[row][col]) ||
                ibz_bitsize(&sk->mat_BAcan_to_BA0_two[row][col]) > TORSION_EVEN_POWER)
                goto fail;
        }
    }

    enc = ibz_to_bytes(enc, &sk->mat_BAcan_to_BA0_two[0][0], TORSION_2POWER_BYTES, false);
    if (enc == NULL) {
        goto fail;
    }
    enc = ibz_to_bytes(enc, &sk->mat_BAcan_to_BA0_two[0][1], TORSION_2POWER_BYTES, false);
    if (enc == NULL) {
        goto fail;
    }
    enc = ibz_to_bytes(enc, &sk->mat_BAcan_to_BA0_two[1][0], TORSION_2POWER_BYTES, false);
    if (enc == NULL) {
        goto fail;
    }
    enc = ibz_to_bytes(enc, &sk->mat_BAcan_to_BA0_two[1][1], TORSION_2POWER_BYTES, false);
    if (enc == NULL) {
        goto fail;
    }

    if ((size_t)(enc - start) != SECRETKEY_BYTES) {
        goto fail;
    }
    return 1;

fail:
    memset(start, 0, SECRETKEY_BYTES);
    return 0;
}

int
secret_key_from_bytes(secret_key_t *sk, public_key_t *pk, const byte_t *enc)
{
#ifndef NDEBUG
    const byte_t *const start = enc;
#endif

    enc = public_key_from_bytes(pk, enc);
    if (enc == NULL) {
        return 0;
    }

    {
        ibz_t norm;
        ibz_init(&norm);
        quat_alg_elem_t gen;
        quat_alg_elem_init(&gen);
        enc = ibz_from_bytes(&norm, enc, FP_ENCODED_BYTES, false);
        enc = ibz_from_bytes(&gen.coord[0], enc, FP_ENCODED_BYTES, true);
        enc = ibz_from_bytes(&gen.coord[1], enc, FP_ENCODED_BYTES, true);
        enc = ibz_from_bytes(&gen.coord[2], enc, FP_ENCODED_BYTES, true);
        enc = ibz_from_bytes(&gen.coord[3], enc, FP_ENCODED_BYTES, true);
        /* paper Issue 9: norm is read from the encoded sk where it was stored
         * as the actual nrd of secret_ideal during keygen serialization. */
        int ok = quat_lideal_create_with_norm(
            &sk->secret_ideal, &gen, &norm, &MAXORD_O0, &QUATALG_PINFTY);
        ibz_finalize(&norm);
        quat_alg_elem_finalize(&gen);
        if (!ok)
            return 0;
    }

    enc = ibz_from_bytes(&sk->mat_BAcan_to_BA0_two[0][0], enc, TORSION_2POWER_BYTES, false);
    enc = ibz_from_bytes(&sk->mat_BAcan_to_BA0_two[0][1], enc, TORSION_2POWER_BYTES, false);
    enc = ibz_from_bytes(&sk->mat_BAcan_to_BA0_two[1][0], enc, TORSION_2POWER_BYTES, false);
    enc = ibz_from_bytes(&sk->mat_BAcan_to_BA0_two[1][1], enc, TORSION_2POWER_BYTES, false);

    for (int row = 0; row < 2; row++) {
        for (int col = 0; col < 2; col++) {
            if (ibz_bitsize(&sk->mat_BAcan_to_BA0_two[row][col]) > TORSION_EVEN_POWER)
                return 0;
        }
    }

    assert(enc - start == SECRETKEY_BYTES);

    sk->curve = pk->curve;
    if (!ec_curve_to_basis_2f_from_hint(
            &sk->canonical_basis, &sk->curve, TORSION_EVEN_POWER, pk->hint_pk)) {
        return 0;
    }
    return 1;
}
