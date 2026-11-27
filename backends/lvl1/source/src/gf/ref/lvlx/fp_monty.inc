// SPDX-License-Identifier: Apache-2.0
// The selected Monty kernel defines Radix, Nlimbs, Nbits and Nbytes.
#include <fp.h>
#include <encoded_sizes.h>

_Static_assert(Nlimbs == NWORDS_FIELD && Wordlength == RADIX, "field configuration mismatch");
_Static_assert(Nbytes <= FP_ENCODED_BYTES, "field encoding is too short");

const digit_t ZERO[NWORDS_FIELD] = { 0 };
const digit_t ONE[NWORDS_FIELD] = MONT_ONE;
static const fp_t TWO_INV = MONT_TWO_INV;
static const fp_t THREE_INV = MONT_THREE_INV;

void
fp_set_small(fp_t *x, digit_t value)
{
    fp_t raw = { 0 };
    for (unsigned offset = 0; offset < Wordlength; offset += Radix)
        raw[offset / Radix] = (value >> offset) & (((digit_t)1 << Radix) - 1);
    nres(raw, *x);
}

void
fp_mul_small(fp_t *x, const fp_t *a, uint32_t value)
{
    fp_t multiplier;
    fp_set_small(&multiplier, value);
    fp_mul(x, a, &multiplier);
}

void fp_set_zero(fp_t *x) { modzer(*x); }
void fp_set_one(fp_t *x) { modcpy(ONE, *x); }
uint32_t fp_is_equal(const fp_t *a, const fp_t *b) { return -(uint32_t)modcmp(*a, *b); }
uint32_t fp_is_zero(const fp_t *a) { return -(uint32_t)modis0(*a); }
void fp_copy(fp_t *x, const fp_t *a) { modcpy(*a, *x); }
void fp_cswap(fp_t *a, fp_t *b, uint32_t ctl) { modcsw(ctl & 1, *a, *b); }
void fp_add(fp_t *x, const fp_t *a, const fp_t *b) { modadd(*a, *b, *x); }
void fp_sub(fp_t *x, const fp_t *a, const fp_t *b) { modsub(*a, *b, *x); }
void fp_neg(fp_t *x, const fp_t *a) { modneg(*a, *x); }
void fp_sqr(fp_t *x, const fp_t *a) { modsqr(*a, *x); }
void fp_mul(fp_t *x, const fp_t *a, const fp_t *b) { modmul(*a, *b, *x); }
void fp_inv(fp_t *x) { modinv(*x, NULL, *x); }
uint32_t fp_is_square(const fp_t *a) { return -(uint32_t)modqr(NULL, *a); }
void fp_sqrt(fp_t *x) { modsqrt(*x, NULL, *x); }
void fp_half(fp_t *x, const fp_t *a) { fp_mul(x, a, &TWO_INV); }
void fp_div3(fp_t *x, const fp_t *a) { fp_mul(x, a, &THREE_INV); }
void fp_exp3div4(fp_t *x, const fp_t *a) { modpro(*a, *x); }

void
fp_encode(void *dst, const fp_t *a)
{
    uint8_t *bytes = dst;
    fp_t canonical;
    redc(*a, canonical);
    for (unsigned i = 0; i < Nbytes; ++i) {
        bytes[i] = (uint8_t)canonical[0];
        modshr(8, canonical);
    }
    memset(bytes + Nbytes, 0, FP_ENCODED_BYTES - Nbytes);
}

uint32_t
fp_decode(fp_t *d, const void *src)
{
    const uint8_t *bytes = src;
    const uint8_t top_mask = UINT8_MAX >> (8 * Nbytes - Nbits);
    uint32_t excess = bytes[Nbytes - 1] & (uint8_t)~top_mask;
    for (unsigned i = Nbytes; i < FP_ENCODED_BYTES; ++i)
        excess |= bytes[i];

    fp_t raw = { 0 };
    raw[0] = bytes[Nbytes - 1] & top_mask;
    for (int i = Nbytes - 2; i >= 0; --i) {
        modshl(8, raw);
        raw[0] |= bytes[i];
    }
    const uint32_t canonical = modfsb(raw) & (excess == 0);
    nres(raw, *d);
    for (unsigned i = 0; i < Nlimbs; ++i)
        (*d)[i] &= (digit_t)0 - canonical;
    return (uint32_t)0 - canonical;
}

void
fp_decode_reduce(fp_t *d, const void *src, size_t len)
{
    const uint8_t *bytes = src;
    fp_t acc = { 0 }, digit, radix;
    fp_set_small(&radix, 256);
    for (size_t i = len; i-- > 0;) {
        fp_set_small(&digit, bytes[i]);
        fp_mul(&acc, &acc, &radix);
        fp_add(&acc, &acc, &digit);
    }
    fp_copy(d, &acc);
}
