#include "sqisign_fp505_19x27.h"

#define NLIMBS SQISIGN_FP505_LIMBS
#define RBITS SQISIGN_FP505_RADIX_BITS
#define RMASK SQISIGN_FP505_MASK

/* p = 27*2^500 - 1 in radix 2^27. */
const fp505_t fp505_modulus = {{
    UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff),
    UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff),
    UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff),
    UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff),
    UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x0006bfff)
}};

const fp505_t fp505_zero = {{
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000),
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000),
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000),
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000),
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000)
}};

/* R mod p, where R=2^513. */
const fp505_t fp505_one = {{
    UINT32_C(0x0000012f), UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000),
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000),
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000),
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000),
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x0002c000)
}};

/* R^2 mod p, represented as an ordinary radix-2^27 integer. */
const fp505_t fp505_r2 = {{
    UINT32_C(0x0343c668), UINT32_C(0x04bda12f), UINT32_C(0x03425ed0), UINT32_C(0x04bda12f),
    UINT32_C(0x03425ed0), UINT32_C(0x04bda12f), UINT32_C(0x03425ed0), UINT32_C(0x04bda12f),
    UINT32_C(0x03425ed0), UINT32_C(0x04bda12f), UINT32_C(0x03425ed0), UINT32_C(0x04bda12f),
    UINT32_C(0x03425ed0), UINT32_C(0x04bda12f), UINT32_C(0x03425ed0), UINT32_C(0x04bda12f),
    UINT32_C(0x03425ed0), UINT32_C(0x04bda12f), UINT32_C(0x00045ed0)
}};

/* Standard-domain integer one. */
static const fp505_t fp505_std_one = {{
    UINT32_C(0x00000001), UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000),
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000),
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000),
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000),
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000)
}};

/* Exponents in little-endian radix-2^27, standard integer form. */
static const uint32_t exponent_p_minus_2[NLIMBS] = {
    UINT32_C(0x07fffffd), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff),
    UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff),
    UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff),
    UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff),
    UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x0006bfff)
};

static const uint32_t exponent_sqrt[NLIMBS] = {
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000),
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000),
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000),
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x00000000),
    UINT32_C(0x00000000), UINT32_C(0x00000000), UINT32_C(0x0001b000)
};

static const uint32_t exponent_legendre[NLIMBS] = {
    UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff),
    UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff),
    UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff),
    UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x07ffffff),
    UINT32_C(0x07ffffff), UINT32_C(0x07ffffff), UINT32_C(0x00035fff)
};

static uint32_t ct_mask_u32(uint32_t bit)
{
    return UINT32_C(0) - (bit & UINT32_C(1));
}

static uint32_t ct_is_nonzero_u32(uint32_t x)
{
    return (x | (UINT32_C(0) - x)) >> 31;
}

static void select_limbs(uint32_t out[NLIMBS],
                         const uint32_t x[NLIMBS],
                         const uint32_t y[NLIMBS],
                         uint32_t choose_x)
{
    const uint32_t m = ct_mask_u32(choose_x);
    size_t i;
    for (i = 0; i < NLIMBS; i++) {
        out[i] = (x[i] & m) | (y[i] & ~m);
    }
}

/* Compare canonical radix-2^27 values. Returns 1 iff a >= b. */
static uint32_t limbs_ge(const uint32_t a[NLIMBS],
                         const uint32_t b[NLIMBS])
{
    uint64_t borrow = 0;
    size_t i;
    for (i = 0; i < NLIMBS; i++) {
        const uint64_t bi = (uint64_t)b[i] + borrow;
        borrow = ((uint64_t)a[i] < bi);
    }
    return (uint32_t)(borrow ^ UINT64_C(1));
}

/* r = a-p if a>=p, otherwise a. Input must fit in 19 radix limbs. */
static void conditional_sub_p(uint32_t r[NLIMBS],
                              const uint32_t a[NLIMBS])
{
    uint32_t d[NLIMBS];
    uint64_t borrow = 0;
    size_t i;
    for (i = 0; i < NLIMBS; i++) {
        const uint64_t pi = (uint64_t)fp505_modulus.limb[i] + borrow;
        const uint64_t ai = (uint64_t)a[i];
        d[i] = (uint32_t)(ai - pi) & RMASK;
        borrow = (ai < pi);
    }
    select_limbs(r, d, a, (uint32_t)(borrow ^ UINT64_C(1)));
}

/*
 * CIOS Montgomery multiplication in radix 2^27.
 *
 * Because p[0] = 2^27-1 == -1 mod 2^27, the Montgomery constant
 * -p^{-1} mod 2^27 is exactly 1.  Thus m=t[0] at each iteration.
 * Inputs and output are canonical 19-limb residues.
 */
static void montgomery_mul(uint32_t out[NLIMBS],
                           const uint32_t a[NLIMBS],
                           const uint32_t b[NLIMBS])
{
    uint32_t t[NLIMBS + 1u] = {0};
    size_t i, j;

    for (i = 0; i < NLIMBS; i++) {
        uint64_t z;
        uint64_t carry;
        uint32_t m;

        z = (uint64_t)t[0] + (uint64_t)a[0] * (uint64_t)b[i];
        m = (uint32_t)z & RMASK;

        /* The low limb cancels modulo 2^27. */
        z += (uint64_t)m * (uint64_t)fp505_modulus.limb[0];
        carry = z >> RBITS;

        for (j = 1; j < NLIMBS; j++) {
            z = (uint64_t)t[j]
              + (uint64_t)a[j] * (uint64_t)b[i]
              + (uint64_t)m * (uint64_t)fp505_modulus.limb[j]
              + carry;
            t[j - 1u] = (uint32_t)z & RMASK;
            carry = z >> RBITS;
        }

        z = (uint64_t)t[NLIMBS] + carry;
        t[NLIMBS - 1u] = (uint32_t)z & RMASK;
        t[NLIMBS] = (uint32_t)(z >> RBITS);
    }

    /*
     * With canonical inputs the CIOS result is below 2p and t[NLIMBS] is
     * zero because 2p << R.  Retain a defensive fold for malformed internal
     * values: one high radix limb equals R and R mod p is fp505_one.
     */
    {
        /*
         * Fold the defensive high radix limb without a data-dependent branch.
         * It is zero for canonical public-API inputs, but the fixed schedule
         * also keeps internal fault behaviour deterministic.
         */
        uint64_t carry = 0;
        const uint32_t high = t[NLIMBS];
        for (i = 0; i < NLIMBS; i++) {
            const uint64_t z = (uint64_t)t[i]
                             + (uint64_t)high * fp505_one.limb[i]
                             + carry;
            t[i] = (uint32_t)z & RMASK;
            carry = z >> RBITS;
        }
        (void)carry;
    }

    conditional_sub_p(out, t);
}

static void pow_fixed(fp505_t *r, const fp505_t *a,
                      const uint32_t exponent[NLIMBS], unsigned top_bit)
{
    fp505_t acc = fp505_one;
    fp505_t base = *a;
    int limb = (int)(top_bit / RBITS);
    int off = (int)(top_bit - (unsigned)limb * RBITS);

    /*
     * Walk radix limbs directly.  This avoids a variable integer division by
     * 29 in every exponent bit; loop bounds depend only on fixed exponents.
     */
    for (; limb >= 0; limb--) {
        for (; off >= 0; off--) {
            fp505_t sq;
            fp505_t mul;
            const uint32_t ebit =
                (exponent[(unsigned)limb] >> (unsigned)off) & 1u;

            montgomery_mul(sq.limb, acc.limb, acc.limb);
            montgomery_mul(mul.limb, sq.limb, base.limb);
            select_limbs(acc.limb, mul.limb, sq.limb, ebit);
        }
        off = (int)RBITS - 1;
    }

    *r = acc;
    fp505_secure_clear(&base);
    fp505_secure_clear(&acc);
}

/* Scalar incomplete reduction for the 18x27 + 19 layout. */
static void incomplete_reduce_limbs(uint32_t x[NLIMBS])
{
    uint32_t carry = 0;
    size_t i;
    for (i = 0; i < NLIMBS - 1u; i++) {
        const uint64_t z = (uint64_t)x[i] + carry;
        x[i] = (uint32_t)z & RMASK;
        carry = (uint32_t)(z >> RBITS);
    }
    x[NLIMBS - 1u] += carry;
}

/* Exact floor(x/27) without a hardware division instruction. */
static uint32_t div_coeff_u32(uint32_t x)
{
    const uint32_t reciprocal = UINT32_C(0x097b425e);
    uint32_t q = (uint32_t)(((uint64_t)x * reciprocal) >> 32);
    const uint32_t rem = x - q * UINT32_C(27);
    q += (uint32_t)(rem >= UINT32_C(27));
    return q;
}

void fp505_loose_from_canonical(fp505_loose_t *r, const fp505_t *a)
{
    size_t i;
    for (i = 0; i < NLIMBS; i++) r->limb[i] = a->limb[i];
}

void fp505_loose_add_bound4(fp505_loose_t *r, const fp505_loose_t *a,
                            const fp505_loose_t *b)
{
    size_t i;
    for (i = 0; i < NLIMBS; i++) r->limb[i] = a->limb[i] + b->limb[i];
}

void fp505_loose_sub_bound4(fp505_loose_t *r, const fp505_loose_t *a,
                            const fp505_loose_t *b)
{
    size_t i;
    /* Component-wise 4p prevents underflow for inputs below bound 4. */
    for (i = 0; i < NLIMBS; i++) {
        r->limb[i] = a->limb[i]
                   + UINT32_C(4) * fp505_modulus.limb[i]
                   - b->limb[i];
    }
}

void fp505_loose_incomplete_reduce(fp505_loose_t *r,
                                   const fp505_loose_t *a)
{
    size_t i;
    for (i = 0; i < NLIMBS; i++) r->limb[i] = a->limb[i];
    incomplete_reduce_limbs(r->limb);
}

void fp505_loose_complete_reduce(fp505_t *r, const fp505_loose_t *a)
{
    uint32_t x[NLIMBS];
    uint32_t q;
    size_t i;

    for (i = 0; i < NLIMBS; i++) x[i] = a->limb[i];
    incomplete_reduce_limbs(x);

    /*
     * p+1 = 27*2^14*(2^27)^18.  Split the top limb into
     * x[top] = q*(27*2^14) + rem and fold q into limb 0.
     * This is the scalar counterpart of the paper's fixed Barrett fold.
     */
    q = div_coeff_u32(x[NLIMBS - 1u] >> 14);
    x[NLIMBS - 1u] -= q * UINT32_C(0x0006c000);
    x[0] += q;
    incomplete_reduce_limbs(x);

    conditional_sub_p(r->limb, x);
    for (i = 0; i < NLIMBS; i++) x[i] = 0;
}

void fp505_set_zero(fp505_t *r) { *r = fp505_zero; }
void fp505_set_one(fp505_t *r) { *r = fp505_one; }
void fp505_copy(fp505_t *r, const fp505_t *a) { *r = *a; }

void fp505_from_u32(fp505_t *r, uint32_t x)
{
    fp505_t s = fp505_zero;
    s.limb[0] = x & RMASK;
    s.limb[1] = x >> RBITS;
    montgomery_mul(r->limb, s.limb, fp505_r2.limb);
}

int fp505_decode(fp505_t *r, const uint8_t in[SQISIGN_FP505_ENCODED_BYTES])
{
    fp505_t s = fp505_zero;
    uint64_t acc = 0;
    unsigned acc_bits = 0;
    size_t limb_index = 0;
    size_t i;

    for (i = 0; i < SQISIGN_FP505_ENCODED_BYTES; i++) {
        acc |= (uint64_t)in[i] << acc_bits;
        acc_bits += 8u;
        if (acc_bits >= RBITS && limb_index < NLIMBS) {
            s.limb[limb_index++] = (uint32_t)acc & RMASK;
            acc >>= RBITS;
            acc_bits -= RBITS;
        }
    }
    if (limb_index < NLIMBS) {
        s.limb[limb_index++] = (uint32_t)acc & RMASK;
        acc >>= RBITS;
        acc_bits = (acc_bits > RBITS) ? (acc_bits - RBITS) : 0u;
    }

    /* Any residual high bits, or s>=p, are non-canonical. */
    if (acc != 0u || limb_index != NLIMBS ||
        limbs_ge(s.limb, fp505_modulus.limb)) {
        fp505_set_zero(r);
        fp505_secure_clear(&s);
        return 0;
    }

    montgomery_mul(r->limb, s.limb, fp505_r2.limb);
    fp505_secure_clear(&s);
    return 1;
}

void fp505_encode(uint8_t out[SQISIGN_FP505_ENCODED_BYTES], const fp505_t *a)
{
    fp505_t s;
    uint64_t acc = 0;
    unsigned acc_bits = 0;
    size_t limb_index = 0;
    size_t i;

    montgomery_mul(s.limb, a->limb, fp505_std_one.limb);

    for (i = 0; i < SQISIGN_FP505_ENCODED_BYTES; i++) {
        while (acc_bits < 8u && limb_index < NLIMBS) {
            acc |= (uint64_t)s.limb[limb_index++] << acc_bits;
            acc_bits += RBITS;
        }
        out[i] = (uint8_t)acc;
        acc >>= 8;
        acc_bits -= 8u;
    }
    fp505_secure_clear(&s);
}

void fp505_add(fp505_t *r, const fp505_t *a, const fp505_t *b)
{
    uint32_t sum[NLIMBS];
    uint64_t carry = 0;
    size_t i;
    for (i = 0; i < NLIMBS; i++) {
        const uint64_t z = (uint64_t)a->limb[i] + b->limb[i] + carry;
        sum[i] = (uint32_t)z & RMASK;
        carry = z >> RBITS;
    }

    /* Since a+b<2p<R, carry is zero. */
    conditional_sub_p(r->limb, sum);
}

void fp505_sub(fp505_t *r, const fp505_t *a, const fp505_t *b)
{
    uint32_t d[NLIMBS];
    uint32_t corrected[NLIMBS];
    uint64_t borrow = 0;
    uint64_t carry = 0;
    size_t i;

    for (i = 0; i < NLIMBS; i++) {
        const uint64_t bi = (uint64_t)b->limb[i] + borrow;
        const uint64_t ai = a->limb[i];
        d[i] = (uint32_t)(ai - bi) & RMASK;
        borrow = (ai < bi);
    }

    /* Add p iff the radix subtraction borrowed. */
    {
        const uint32_t m = ct_mask_u32((uint32_t)borrow);
        for (i = 0; i < NLIMBS; i++) {
            const uint64_t z = (uint64_t)d[i]
                             + (fp505_modulus.limb[i] & m)
                             + carry;
            corrected[i] = (uint32_t)z & RMASK;
            carry = z >> RBITS;
        }
    }
    for (i = 0; i < NLIMBS; i++) r->limb[i] = corrected[i];
}

void fp505_neg(fp505_t *r, const fp505_t *a)
{
    fp505_sub(r, &fp505_zero, a);
}

void fp505_double(fp505_t *r, const fp505_t *a)
{
    fp505_add(r, a, a);
}

void fp505_mul(fp505_t *r, const fp505_t *a, const fp505_t *b)
{
    uint32_t out[NLIMBS];
    montgomery_mul(out, a->limb, b->limb);
    {
        size_t i;
        for (i = 0; i < NLIMBS; i++) r->limb[i] = out[i];
    }
}

void fp505_sqr(fp505_t *r, const fp505_t *a)
{
    fp505_mul(r, a, a);
}

int fp505_inv(fp505_t *r, const fp505_t *a)
{
    fp505_t candidate;
    const uint32_t ok = (uint32_t)(fp505_iszero(a) ^ 1);
    pow_fixed(&candidate, a, exponent_p_minus_2, 504u);
    select_limbs(r->limb, candidate.limb, fp505_zero.limb, ok);
    fp505_secure_clear(&candidate);
    return (int)ok;
}

int fp505_sqrt(fp505_t *r, const fp505_t *a)
{
    fp505_t candidate;
    fp505_t check;
    uint32_t ok;
    pow_fixed(&candidate, a, exponent_sqrt, 502u);
    fp505_sqr(&check, &candidate);
    ok = (uint32_t)fp505_equal(&check, a);
    select_limbs(r->limb, candidate.limb, fp505_zero.limb, ok);
    fp505_secure_clear(&candidate);
    fp505_secure_clear(&check);
    return (int)ok;
}

int fp505_legendre(const fp505_t *a)
{
    fp505_t t;
    const uint32_t nz = (uint32_t)(fp505_iszero(a) ^ 1);
    uint32_t is_square;
    int result;
    pow_fixed(&t, a, exponent_legendre, 503u);
    is_square = (uint32_t)fp505_equal(&t, &fp505_one);
    /* 0 for zero, +1 for square, -1 for non-square. */
    result = (int)nz * ((int)(is_square << 1) - 1);
    fp505_secure_clear(&t);
    return result;
}

int fp505_iszero(const fp505_t *a)
{
    uint32_t x = 0;
    size_t i;
    for (i = 0; i < NLIMBS; i++) x |= a->limb[i];
    return (int)(ct_is_nonzero_u32(x) ^ 1u);
}

int fp505_equal(const fp505_t *a, const fp505_t *b)
{
    uint32_t x = 0;
    size_t i;
    for (i = 0; i < NLIMBS; i++) x |= a->limb[i] ^ b->limb[i];
    return (int)(ct_is_nonzero_u32(x) ^ 1u);
}

void fp505_secure_clear(fp505_t *a)
{
    volatile uint32_t *p = (volatile uint32_t *)a->limb;
    size_t i;
    for (i = 0; i < NLIMBS; i++) p[i] = 0;
}

void fp2_505_set_zero(fp2_505_t *r)
{
    fp505_set_zero(&r->re);
    fp505_set_zero(&r->im);
}

void fp2_505_set_one(fp2_505_t *r)
{
    fp505_set_one(&r->re);
    fp505_set_zero(&r->im);
}

void fp2_505_copy(fp2_505_t *r, const fp2_505_t *a) { *r = *a; }

void fp2_505_add(fp2_505_t *r, const fp2_505_t *a, const fp2_505_t *b)
{
    fp505_t re, im;
    fp505_add(&re, &a->re, &b->re);
    fp505_add(&im, &a->im, &b->im);
    r->re = re;
    r->im = im;
}

void fp2_505_sub(fp2_505_t *r, const fp2_505_t *a, const fp2_505_t *b)
{
    fp505_t re, im;
    fp505_sub(&re, &a->re, &b->re);
    fp505_sub(&im, &a->im, &b->im);
    r->re = re;
    r->im = im;
}

void fp2_505_neg(fp2_505_t *r, const fp2_505_t *a)
{
    fp505_t re, im;
    fp505_neg(&re, &a->re);
    fp505_neg(&im, &a->im);
    r->re = re;
    r->im = im;
}

void fp2_505_mul(fp2_505_t *r, const fp2_505_t *a, const fp2_505_t *b)
{
    fp505_t ac, bd, apb, cpd, cross, re, im;
    fp505_mul(&ac, &a->re, &b->re);
    fp505_mul(&bd, &a->im, &b->im);
    fp505_add(&apb, &a->re, &a->im);
    fp505_add(&cpd, &b->re, &b->im);
    fp505_mul(&cross, &apb, &cpd);
    fp505_sub(&re, &ac, &bd);
    fp505_sub(&im, &cross, &ac);
    fp505_sub(&im, &im, &bd);
    r->re = re;
    r->im = im;
    fp505_secure_clear(&ac);
    fp505_secure_clear(&bd);
    fp505_secure_clear(&apb);
    fp505_secure_clear(&cpd);
    fp505_secure_clear(&cross);
}

void fp2_505_sqr(fp2_505_t *r, const fp2_505_t *a)
{
    fp505_t sum, diff, re, im;
    fp505_add(&sum, &a->re, &a->im);
    fp505_sub(&diff, &a->re, &a->im);
    fp505_mul(&re, &sum, &diff);
    fp505_mul(&im, &a->re, &a->im);
    fp505_double(&im, &im);
    r->re = re;
    r->im = im;
    fp505_secure_clear(&sum);
    fp505_secure_clear(&diff);
}

int fp2_505_inv(fp2_505_t *r, const fp2_505_t *a)
{
    fp505_t rr, ii, den, den_inv, re, im;
    int ok;

    fp505_sqr(&rr, &a->re);
    fp505_sqr(&ii, &a->im);
    fp505_add(&den, &rr, &ii);

    /*
     * fp505_inv() writes zero for a zero denominator.  Continue through the
     * same multiplication schedule in both cases, so inversion failure does
     * not introduce a secret-dependent control-flow edge here.
     */
    ok = fp505_inv(&den_inv, &den);
    fp505_mul(&re, &a->re, &den_inv);
    fp505_mul(&im, &a->im, &den_inv);
    fp505_neg(&im, &im);
    r->re = re;
    r->im = im;

    fp505_secure_clear(&rr);
    fp505_secure_clear(&ii);
    fp505_secure_clear(&den);
    fp505_secure_clear(&den_inv);
    fp505_secure_clear(&re);
    fp505_secure_clear(&im);
    return ok;
}

int fp2_505_iszero(const fp2_505_t *a)
{
    return fp505_iszero(&a->re) & fp505_iszero(&a->im);
}

int fp2_505_equal(const fp2_505_t *a, const fp2_505_t *b)
{
    return fp505_equal(&a->re, &b->re) & fp505_equal(&a->im, &b->im);
}

void fp2_505_secure_clear(fp2_505_t *a)
{
    fp505_secure_clear(&a->re);
    fp505_secure_clear(&a->im);
}

void fp505_mul_2way(fp505_t out[2], const fp505_t a[2], const fp505_t b[2])
{
    fp505_t t0, t1;
    fp505_mul(&t0, &a[0], &b[0]);
    fp505_mul(&t1, &a[1], &b[1]);
    out[0] = t0;
    out[1] = t1;
}

void fp2_505_mul_2way(fp2_505_t out[2], const fp2_505_t a[2],
                      const fp2_505_t b[2])
{
    fp2_505_t t0, t1;
    fp2_505_mul(&t0, &a[0], &b[0]);
    fp2_505_mul(&t1, &a[1], &b[1]);
    out[0] = t0;
    out[1] = t1;
}
