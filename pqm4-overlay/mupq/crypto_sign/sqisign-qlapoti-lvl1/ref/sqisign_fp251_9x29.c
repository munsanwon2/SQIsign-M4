#include "sqisign_fp251_9x29.h"

#define NLIMBS SQISIGN_FP251_LIMBS
#define RBITS SQISIGN_FP251_RADIX_BITS
#define RMASK SQISIGN_FP251_MASK

/* p = 5*2^248 - 1 in radix 2^29. */
const fp251_t fp251_modulus = {{
    UINT32_C(0x1fffffff), UINT32_C(0x1fffffff),
    UINT32_C(0x1fffffff), UINT32_C(0x1fffffff),
    UINT32_C(0x1fffffff), UINT32_C(0x1fffffff),
    UINT32_C(0x1fffffff), UINT32_C(0x1fffffff),
    UINT32_C(0x0004ffff)
}};

const fp251_t fp251_zero = {{0,0,0,0,0,0,0,0,0}};

/* R mod p, where R=2^261. */
const fp251_t fp251_one = {{
    UINT32_C(0x00000666), 0, 0, 0, 0, 0, 0, 0,
    UINT32_C(0x00020000)
}};

/* R^2 mod p, represented as an ordinary radix-2^29 integer. */
const fp251_t fp251_r2 = {{
    UINT32_C(0x0cf5c28f), UINT32_C(0x06666666),
    UINT32_C(0x13333333), UINT32_C(0x19999999),
    UINT32_C(0x0ccccccc), UINT32_C(0x06666666),
    UINT32_C(0x13333333), UINT32_C(0x19999999),
    UINT32_C(0x0001cccc)
}};

/* Standard-domain integer one. */
static const fp251_t fp251_std_one = {{1,0,0,0,0,0,0,0,0}};

/* Exponents in little-endian radix-2^29, standard integer form. */
static const uint32_t exponent_p_minus_2[NLIMBS] = {
    UINT32_C(0x1ffffffd), UINT32_C(0x1fffffff),
    UINT32_C(0x1fffffff), UINT32_C(0x1fffffff),
    UINT32_C(0x1fffffff), UINT32_C(0x1fffffff),
    UINT32_C(0x1fffffff), UINT32_C(0x1fffffff),
    UINT32_C(0x0004ffff)
};

static const uint32_t exponent_sqrt[NLIMBS] = {
    0,0,0,0,0,0,0,0,UINT32_C(0x00014000)
};

static const uint32_t exponent_legendre[NLIMBS] = {
    UINT32_C(0x1fffffff), UINT32_C(0x1fffffff),
    UINT32_C(0x1fffffff), UINT32_C(0x1fffffff),
    UINT32_C(0x1fffffff), UINT32_C(0x1fffffff),
    UINT32_C(0x1fffffff), UINT32_C(0x1fffffff),
    UINT32_C(0x00027fff)
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

/* Compare canonical radix-2^29 values. Returns 1 iff a >= b. */
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

/* r = a-p if a>=p, otherwise a. Input must fit in nine radix limbs. */
static void conditional_sub_p(uint32_t r[NLIMBS],
                              const uint32_t a[NLIMBS])
{
    uint32_t d[NLIMBS];
    uint64_t borrow = 0;
    size_t i;
    for (i = 0; i < NLIMBS; i++) {
        const uint64_t pi = (uint64_t)fp251_modulus.limb[i] + borrow;
        const uint64_t ai = (uint64_t)a[i];
        d[i] = (uint32_t)(ai - pi) & RMASK;
        borrow = (ai < pi);
    }
    select_limbs(r, d, a, (uint32_t)(borrow ^ UINT64_C(1)));
}

/*
 * CIOS Montgomery multiplication in radix 2^29.
 *
 * Because p[0] = 2^29-1 == -1 mod 2^29, the Montgomery constant
 * -p^{-1} mod 2^29 is exactly 1.  Thus m=t[0] at each iteration.
 * Inputs and output are canonical nine-limb residues.
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

        /* The low limb cancels modulo 2^29. */
        z += (uint64_t)m * (uint64_t)fp251_modulus.limb[0];
        carry = z >> RBITS;

        for (j = 1; j < NLIMBS; j++) {
            z = (uint64_t)t[j]
              + (uint64_t)a[j] * (uint64_t)b[i]
              + (uint64_t)m * (uint64_t)fp251_modulus.limb[j]
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
     * values: one high radix limb equals R and R mod p is fp251_one.
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
                             + (uint64_t)high * fp251_one.limb[i]
                             + carry;
            t[i] = (uint32_t)z & RMASK;
            carry = z >> RBITS;
        }
        (void)carry;
    }

    conditional_sub_p(out, t);
}

static void pow_fixed(fp251_t *r, const fp251_t *a,
                      const uint32_t exponent[NLIMBS], unsigned top_bit)
{
    fp251_t acc = fp251_one;
    fp251_t base = *a;
    int limb = (int)(top_bit / RBITS);
    int off = (int)(top_bit - (unsigned)limb * RBITS);

    /*
     * Walk radix limbs directly.  This avoids a variable integer division by
     * 29 in every exponent bit; loop bounds depend only on fixed exponents.
     */
    for (; limb >= 0; limb--) {
        for (; off >= 0; off--) {
            fp251_t sq;
            fp251_t mul;
            const uint32_t ebit =
                (exponent[(unsigned)limb] >> (unsigned)off) & 1u;

            montgomery_mul(sq.limb, acc.limb, acc.limb);
            montgomery_mul(mul.limb, sq.limb, base.limb);
            select_limbs(acc.limb, mul.limb, sq.limb, ebit);
        }
        off = (int)RBITS - 1;
    }

    *r = acc;
    fp251_secure_clear(&base);
    fp251_secure_clear(&acc);
}

/* Scalar incomplete reduction for the 8x29 + 19 layout. */
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

/* Exact floor(x/5) for a 32-bit x, using a fixed reciprocal multiply. */
static uint32_t div5_u32(uint32_t x)
{
    return (uint32_t)(((uint64_t)x * UINT64_C(0xcccccccd)) >> 34);
}

void fp251_loose_from_canonical(fp251_loose_t *r, const fp251_t *a)
{
    size_t i;
    for (i = 0; i < NLIMBS; i++) r->limb[i] = a->limb[i];
}

void fp251_loose_add_bound4(fp251_loose_t *r, const fp251_loose_t *a,
                            const fp251_loose_t *b)
{
    size_t i;
    for (i = 0; i < NLIMBS; i++) r->limb[i] = a->limb[i] + b->limb[i];
}

void fp251_loose_sub_bound4(fp251_loose_t *r, const fp251_loose_t *a,
                            const fp251_loose_t *b)
{
    size_t i;
    /* Component-wise 4p prevents underflow for inputs below bound 4. */
    for (i = 0; i < NLIMBS; i++) {
        r->limb[i] = a->limb[i]
                   + UINT32_C(4) * fp251_modulus.limb[i]
                   - b->limb[i];
    }
}

void fp251_loose_incomplete_reduce(fp251_loose_t *r,
                                   const fp251_loose_t *a)
{
    size_t i;
    for (i = 0; i < NLIMBS; i++) r->limb[i] = a->limb[i];
    incomplete_reduce_limbs(r->limb);
}

void fp251_loose_complete_reduce(fp251_t *r, const fp251_loose_t *a)
{
    uint32_t x[NLIMBS];
    uint32_t q;
    size_t i;

    for (i = 0; i < NLIMBS; i++) x[i] = a->limb[i];
    incomplete_reduce_limbs(x);

    /*
     * p+1 = 5*2^16*(2^29)^8.  Split x[8] into
     * x[8] = q*(5*2^16) + rem, and fold q back into limb 0.
     * The paper computes q by Barrett reduction after shifting by 16.
     */
    q = div5_u32(x[NLIMBS - 1u] >> 16);
    x[NLIMBS - 1u] -= q * UINT32_C(0x00050000);
    x[0] += q;
    incomplete_reduce_limbs(x);

    conditional_sub_p(r->limb, x);
    for (i = 0; i < NLIMBS; i++) x[i] = 0;
}

void fp251_set_zero(fp251_t *r) { *r = fp251_zero; }
void fp251_set_one(fp251_t *r) { *r = fp251_one; }
void fp251_copy(fp251_t *r, const fp251_t *a) { *r = *a; }

void fp251_from_u32(fp251_t *r, uint32_t x)
{
    fp251_t s = fp251_zero;
    s.limb[0] = x & RMASK;
    s.limb[1] = x >> RBITS;
    montgomery_mul(r->limb, s.limb, fp251_r2.limb);
}

int fp251_decode(fp251_t *r, const uint8_t in[SQISIGN_FP251_ENCODED_BYTES])
{
    fp251_t s = fp251_zero;
    uint64_t acc = 0;
    unsigned acc_bits = 0;
    size_t limb_index = 0;
    size_t i;

    for (i = 0; i < SQISIGN_FP251_ENCODED_BYTES; i++) {
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
        limbs_ge(s.limb, fp251_modulus.limb)) {
        fp251_set_zero(r);
        fp251_secure_clear(&s);
        return 0;
    }

    montgomery_mul(r->limb, s.limb, fp251_r2.limb);
    fp251_secure_clear(&s);
    return 1;
}

void fp251_encode(uint8_t out[SQISIGN_FP251_ENCODED_BYTES], const fp251_t *a)
{
    fp251_t s;
    uint64_t acc = 0;
    unsigned acc_bits = 0;
    size_t limb_index = 0;
    size_t i;

    montgomery_mul(s.limb, a->limb, fp251_std_one.limb);

    for (i = 0; i < SQISIGN_FP251_ENCODED_BYTES; i++) {
        while (acc_bits < 8u && limb_index < NLIMBS) {
            acc |= (uint64_t)s.limb[limb_index++] << acc_bits;
            acc_bits += RBITS;
        }
        out[i] = (uint8_t)acc;
        acc >>= 8;
        acc_bits -= 8u;
    }
    fp251_secure_clear(&s);
}

void fp251_add(fp251_t *r, const fp251_t *a, const fp251_t *b)
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

void fp251_sub(fp251_t *r, const fp251_t *a, const fp251_t *b)
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
                             + (fp251_modulus.limb[i] & m)
                             + carry;
            corrected[i] = (uint32_t)z & RMASK;
            carry = z >> RBITS;
        }
    }
    for (i = 0; i < NLIMBS; i++) r->limb[i] = corrected[i];
}

void fp251_neg(fp251_t *r, const fp251_t *a)
{
    fp251_sub(r, &fp251_zero, a);
}

void fp251_double(fp251_t *r, const fp251_t *a)
{
    fp251_add(r, a, a);
}

void fp251_mul(fp251_t *r, const fp251_t *a, const fp251_t *b)
{
    uint32_t out[NLIMBS];
    montgomery_mul(out, a->limb, b->limb);
    {
        size_t i;
        for (i = 0; i < NLIMBS; i++) r->limb[i] = out[i];
    }
}

void fp251_sqr(fp251_t *r, const fp251_t *a)
{
    fp251_mul(r, a, a);
}

int fp251_inv(fp251_t *r, const fp251_t *a)
{
    fp251_t candidate;
    const uint32_t ok = (uint32_t)(fp251_iszero(a) ^ 1);
    pow_fixed(&candidate, a, exponent_p_minus_2, 250u);
    select_limbs(r->limb, candidate.limb, fp251_zero.limb, ok);
    fp251_secure_clear(&candidate);
    return (int)ok;
}

int fp251_sqrt(fp251_t *r, const fp251_t *a)
{
    fp251_t candidate;
    fp251_t check;
    uint32_t ok;
    pow_fixed(&candidate, a, exponent_sqrt, 248u);
    fp251_sqr(&check, &candidate);
    ok = (uint32_t)fp251_equal(&check, a);
    select_limbs(r->limb, candidate.limb, fp251_zero.limb, ok);
    fp251_secure_clear(&candidate);
    fp251_secure_clear(&check);
    return (int)ok;
}

int fp251_legendre(const fp251_t *a)
{
    fp251_t t;
    const uint32_t nz = (uint32_t)(fp251_iszero(a) ^ 1);
    uint32_t is_square;
    int result;
    pow_fixed(&t, a, exponent_legendre, 249u);
    is_square = (uint32_t)fp251_equal(&t, &fp251_one);
    /* 0 for zero, +1 for square, -1 for non-square. */
    result = (int)nz * ((int)(is_square << 1) - 1);
    fp251_secure_clear(&t);
    return result;
}

int fp251_iszero(const fp251_t *a)
{
    uint32_t x = 0;
    size_t i;
    for (i = 0; i < NLIMBS; i++) x |= a->limb[i];
    return (int)(ct_is_nonzero_u32(x) ^ 1u);
}

int fp251_equal(const fp251_t *a, const fp251_t *b)
{
    uint32_t x = 0;
    size_t i;
    for (i = 0; i < NLIMBS; i++) x |= a->limb[i] ^ b->limb[i];
    return (int)(ct_is_nonzero_u32(x) ^ 1u);
}

void fp251_secure_clear(fp251_t *a)
{
    volatile uint32_t *p = (volatile uint32_t *)a->limb;
    size_t i;
    for (i = 0; i < NLIMBS; i++) p[i] = 0;
}

void fp2_251_set_zero(fp2_251_t *r)
{
    fp251_set_zero(&r->re);
    fp251_set_zero(&r->im);
}

void fp2_251_set_one(fp2_251_t *r)
{
    fp251_set_one(&r->re);
    fp251_set_zero(&r->im);
}

void fp2_251_copy(fp2_251_t *r, const fp2_251_t *a) { *r = *a; }

void fp2_251_add(fp2_251_t *r, const fp2_251_t *a, const fp2_251_t *b)
{
    fp251_t re, im;
    fp251_add(&re, &a->re, &b->re);
    fp251_add(&im, &a->im, &b->im);
    r->re = re;
    r->im = im;
}

void fp2_251_sub(fp2_251_t *r, const fp2_251_t *a, const fp2_251_t *b)
{
    fp251_t re, im;
    fp251_sub(&re, &a->re, &b->re);
    fp251_sub(&im, &a->im, &b->im);
    r->re = re;
    r->im = im;
}

void fp2_251_neg(fp2_251_t *r, const fp2_251_t *a)
{
    fp251_t re, im;
    fp251_neg(&re, &a->re);
    fp251_neg(&im, &a->im);
    r->re = re;
    r->im = im;
}

void fp2_251_mul(fp2_251_t *r, const fp2_251_t *a, const fp2_251_t *b)
{
    fp251_t ac, bd, apb, cpd, cross, re, im;
    fp251_mul(&ac, &a->re, &b->re);
    fp251_mul(&bd, &a->im, &b->im);
    fp251_add(&apb, &a->re, &a->im);
    fp251_add(&cpd, &b->re, &b->im);
    fp251_mul(&cross, &apb, &cpd);
    fp251_sub(&re, &ac, &bd);
    fp251_sub(&im, &cross, &ac);
    fp251_sub(&im, &im, &bd);
    r->re = re;
    r->im = im;
    fp251_secure_clear(&ac);
    fp251_secure_clear(&bd);
    fp251_secure_clear(&apb);
    fp251_secure_clear(&cpd);
    fp251_secure_clear(&cross);
}

void fp2_251_sqr(fp2_251_t *r, const fp2_251_t *a)
{
    fp251_t sum, diff, re, im;
    fp251_add(&sum, &a->re, &a->im);
    fp251_sub(&diff, &a->re, &a->im);
    fp251_mul(&re, &sum, &diff);
    fp251_mul(&im, &a->re, &a->im);
    fp251_double(&im, &im);
    r->re = re;
    r->im = im;
    fp251_secure_clear(&sum);
    fp251_secure_clear(&diff);
}

int fp2_251_inv(fp2_251_t *r, const fp2_251_t *a)
{
    fp251_t rr, ii, den, den_inv, re, im;
    int ok;

    fp251_sqr(&rr, &a->re);
    fp251_sqr(&ii, &a->im);
    fp251_add(&den, &rr, &ii);

    /*
     * fp251_inv() writes zero for a zero denominator.  Continue through the
     * same multiplication schedule in both cases, so inversion failure does
     * not introduce a secret-dependent control-flow edge here.
     */
    ok = fp251_inv(&den_inv, &den);
    fp251_mul(&re, &a->re, &den_inv);
    fp251_mul(&im, &a->im, &den_inv);
    fp251_neg(&im, &im);
    r->re = re;
    r->im = im;

    fp251_secure_clear(&rr);
    fp251_secure_clear(&ii);
    fp251_secure_clear(&den);
    fp251_secure_clear(&den_inv);
    fp251_secure_clear(&re);
    fp251_secure_clear(&im);
    return ok;
}

int fp2_251_iszero(const fp2_251_t *a)
{
    return fp251_iszero(&a->re) & fp251_iszero(&a->im);
}

int fp2_251_equal(const fp2_251_t *a, const fp2_251_t *b)
{
    return fp251_equal(&a->re, &b->re) & fp251_equal(&a->im, &b->im);
}

void fp2_251_secure_clear(fp2_251_t *a)
{
    fp251_secure_clear(&a->re);
    fp251_secure_clear(&a->im);
}

void fp251_mul_2way(fp251_t out[2], const fp251_t a[2], const fp251_t b[2])
{
    fp251_t t0, t1;
    fp251_mul(&t0, &a[0], &b[0]);
    fp251_mul(&t1, &a[1], &b[1]);
    out[0] = t0;
    out[1] = t1;
}

void fp2_251_mul_2way(fp2_251_t out[2], const fp2_251_t a[2],
                      const fp2_251_t b[2])
{
    fp2_251_t t0, t1;
    fp2_251_mul(&t0, &a[0], &b[0]);
    fp2_251_mul(&t1, &a[1], &b[1]);
    out[0] = t0;
    out[1] = t1;
}
