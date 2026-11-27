#include "dpe.h"

int
dpe_get_z(ibz_t *x, const dpe_t y)
{
    dpe_t integer;
    dpe_round(integer, y);
    const int exponent = DPE_EXP(integer);
    const double mantissa = DPE_MANT(integer);

    if (mantissa == 0.0) {
        ibz_set(x, 0);
        return 1;
    }
    if (exponent >= IBZ_BITS) {
        if (exponent != IBZ_BITS || mantissa != -0.5)
            return 0;
        ibz_set(x, 0);
        (*x)[IBZ_LIMBS - 1] = UINT64_C(1) << 63;
        return 1;
    }

    const int kept = exponent < DPE_BITSIZE ? exponent : DPE_BITSIZE;
    const uint64_t significand = (uint64_t)ldexp(fabs(mantissa), kept);
    ibz_set_u64(x, significand);
    ibz_mul_2exp(x, x, exponent - kept);
    if (mantissa < 0.0)
        ibz_neg(x, x);
    return 1;
}

void
dpe_set_z(dpe_t x, const ibz_t *y)
{
    int bits = ibz_bitsize(y);
    if (bits == 0) {
        dpe_set_d(x, 0.0);
        return;
    }

    ibz_t magnitude;
    ibz_abs(&magnitude, y);
    const int kept = bits < DPE_BITSIZE ? bits : DPE_BITSIZE;
    const int discarded = bits - kept;
    const int word = discarded / 64;
    const unsigned offset = discarded % 64;
    uint64_t significand = magnitude[word] >> offset;
    if (offset != 0 && word + 1 < IBZ_LIMBS)
        significand |= magnitude[word + 1] << (64 - offset);

    if (discarded != 0) {
        const int guard_word = (discarded - 1) / 64;
        const uint64_t guard = UINT64_C(1) << ((discarded - 1) % 64);
        int sticky = (magnitude[guard_word] & (guard - 1)) != 0;
        for (int i = 0; i < guard_word; ++i)
            sticky |= magnitude[i] != 0;
        if ((magnitude[guard_word] & guard) && (sticky || (significand & 1)))
            significand++;
    }

    double mantissa = ldexp((double)significand, -kept);
    if (mantissa == 1.0) {
        mantissa = 0.5;
        bits++;
    }
    DPE_MANT(x) = ibz_is_negative(y) ? -mantissa : mantissa;
    DPE_EXP(x) = bits;
}
