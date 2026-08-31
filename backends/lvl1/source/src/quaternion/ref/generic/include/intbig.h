#ifndef INTBIG_H
#define INTBIG_H
#pragma once
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <rng.h>

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef IBZ_LIMBS
/* SQISIGN_VARIANT is defined by the build system as token `lvl1`/`lvl3`/`lvl5`.
 *
 * Note: a previous version compared tokens via `#if SQISIGN_VARIANT == SQISIGN_LVL1`,
 * but the C preprocessor evaluates undefined identifier tokens as 0 in `#if`
 * expressions, so `lvl1 == lvl1` becomes `0 == 0` and ALL THREE BRANCHES
 * matched — meaning every level got IBZ_LIMBS=110. This caused the L3 precomp
 * data file (encoded with 168-limb arrays) to overflow the 110-limb buffer
 * with 'excess elements in array initializer'.
 *
 * Fix: token concatenation. CONCAT(IBZ_LIMBS_, SQISIGN_VARIANT) expands to
 * IBZ_LIMBS_lvl1 / IBZ_LIMBS_lvl3 / IBZ_LIMBS_lvl5, each of which IS defined
 * to the correct number. This selects per-level width via real macro
 * substitution, not the broken token-equality #if. */
#define _IBZ_CONCAT(a, b) a##b
#define _IBZ_CONCAT2(a, b) _IBZ_CONCAT(a, b)

/* Revised worst-case magnitude bounds 1665 / 2521 / 3319, with one sign bit
 * and rounded up to complete 64-bit limbs. */
#define IBZ_LIMBS_lvl1 27  /* NIST-I:   1728 total / 1727 magnitude bits */
#define IBZ_LIMBS_lvl3 40  /* NIST-III: 2560 total / 2559 magnitude bits */
#define IBZ_LIMBS_lvl5 52  /* NIST-V:   3328 total / 3327 magnitude bits */

#ifdef SQISIGN_VARIANT
#define IBZ_LIMBS _IBZ_CONCAT2(IBZ_LIMBS_, SQISIGN_VARIANT)
#else
#define IBZ_LIMBS IBZ_LIMBS_lvl5 /* Conservative standalone default: Level V */
#endif
#endif

/* Keep every route cap within the configured signed storage width. */
#define IBZ_HNF_ROUTE_BITS_lvl1 (27 * 64)
#define IBZ_HNF_ROUTE_BITS_lvl3 (40 * 64)
#define IBZ_HNF_ROUTE_BITS_lvl5 (52 * 64)
#ifdef SQISIGN_VARIANT
#define IBZ_HNF_ROUTE_BITS _IBZ_CONCAT2(IBZ_HNF_ROUTE_BITS_, SQISIGN_VARIANT)
#else
#define IBZ_HNF_ROUTE_BITS IBZ_HNF_ROUTE_BITS_lvl5
#endif

#ifndef IBZ_BITS
#define IBZ_BITS (IBZ_LIMBS * 64)
#endif
    typedef uint64_t ibz_t[IBZ_LIMBS];

    typedef uint64_t limb_t;
    typedef uint32_t half_limb_t;
    /* Fixed-precision integers always use 64-bit storage limbs.  Keep this
     * name distinct from tutil.h's target-dependent `digit_t` macro: the
     * latter is 32 bit on Cortex-M4 and is used by the field/MP backend. */
    typedef uint64_t ibz_digit_t;

    void ibz_init(ibz_t *x);

    void ibz_finalize(ibz_t *x);


    void ibz_copy(ibz_t *target, const ibz_t *value);
    void ibz_swap(ibz_t *a, ibz_t *b);


    void ibz_neg(ibz_t *neg, const ibz_t *a);
    /** Absolute value.  For the sole unrepresentable value -2^(IBZ_BITS-1),
     * its unsigned magnitude bit pattern (the same ibz_t value) is retained. */
    void ibz_abs(ibz_t *abs, const ibz_t *a);

    void ibz_add(ibz_t *sum, const ibz_t *a, const ibz_t *b);
    void ibz_sub(ibz_t *diff, const ibz_t *a, const ibz_t *b);
    void ibz_mul(ibz_t *prod, const ibz_t *a, const ibz_t *b);

    /** Divide by 2^exp, truncating toward zero. */
    void ibz_div_2exp(ibz_t *quotient, const ibz_t *a, uint32_t exp);
    void ibz_mul_2exp(ibz_t *result, const ibz_t *a, size_t shift);


    int ibz_cmp(const ibz_t *a, const ibz_t *b);
    int ibz_is_zero(const ibz_t *x);
    int ibz_is_one(const ibz_t *x);
    int ibz_is_even(const ibz_t *x);
    int ibz_is_odd(const ibz_t *x);
    int ibz_is_negative(const ibz_t *x);

    void ibz_set(ibz_t *i, int32_t x);
    void ibz_set_u64(ibz_t *i, uint64_t x);
    int32_t ibz_get(const ibz_t *i);

    int ibz_cmp_int32(const ibz_t *x, int32_t y);
    int ibz_bitsize(const ibz_t *a);
    /** Number of base-2^64 digits needed for |a| (one for zero). */
    size_t ibz_digits_required(const ibz_t *a);
    /** Export |a| to exactly target_len available digits, clearing unused
     * digits.  Returns zero without a partial value when the buffer is too
     * short, and one on success. */
    int ibz_to_u64_digits_checked(ibz_digit_t *digits, size_t target_len, const ibz_t *a);
    /** Legacy magnitude export.  The caller must provide at least
     * ibz_digits_required(a) digits. */
    void ibz_to_u64_digits(ibz_digit_t *digits, const ibz_t *a);

    /** Export |a| as little-endian 32-bit digits.  This is the conversion
     * boundary used by the Cortex-M4 field/MP backend. */
    int ibz_to_u32_digits_checked(uint32_t *digits, size_t target_len, const ibz_t *a);

    /** @brief generate a random value in the inclusive interval [a, b]
     *  The interval must be ordered and its width must fit the positive
     *  signed ibz_t range.
     * @returns 1 on success, 0 on failure
     */
    int ibz_rand_interval(ibz_t *rand, const ibz_t *a, const ibz_t *b);

    /** @brief generate random value in [-m, m]
     *  assumed that m > 0 and bitlength of m < 32 bit
     * @returns 1 on success, 0 on failiure
     */
    int ibz_rand_interval_minm_m(ibz_t *rand, int32_t m);

#if defined(RADIX_32)
#define ibz_to_mp_digits_checked(T, N, I) ibz_to_u32_digits_checked((uint32_t *)(T), (N), (I))
#define ibz_copy_mp_digits(I, T, N) ibz_copy_u32_digits((I), (const uint32_t *)(T), (N))
#else
#define ibz_to_mp_digits_checked(T, N, I) ibz_to_u64_digits_checked((ibz_digit_t *)(T), (N), (I))
#define ibz_copy_mp_digits(I, T, N) ibz_copy_u64_digits((I), (const ibz_digit_t *)(T), (N))
#endif

#define ibz_to_digit_array(T, I)                                                                                       \
    do {                                                                                                               \
        memset((T), 0, sizeof(T));                                                                                     \
        (void)ibz_to_mp_digits_checked((T), sizeof(T) / sizeof(*(T)), (I));                                             \
    } while (0)
    void ibz_copy_u64_digits(ibz_t *a, const ibz_digit_t *digits, size_t len);
    void ibz_copy_u32_digits(ibz_t *a, const uint32_t *digits, size_t len);
#define ibz_copy_digit_array(I, T)                                                                                     \
    do {                                                                                                               \
        ibz_copy_mp_digits((I), (T), sizeof(T) / sizeof(*(T)));                                                        \
    } while (0)
    size_t ibz_size_in_base(const ibz_t *a, int base);

    /** Position of the first set bit; returns 0 for zero. */
    int ibz_two_adic(const ibz_t *pow);

    /** Truncating division.  Division by zero terminates the process. */
    void ibz_div(ibz_t *quotient, ibz_t *remainder, const ibz_t *a, const ibz_t *b);

    /** Floor division; the nonzero remainder has the divisor's sign.
     * Division by zero terminates the process. */
    void ibz_div_floor(ibz_t *q, ibz_t *r, const ibz_t *n, const ibz_t *d);

    /** Floor remainder.  A zero modulus terminates the process. */
    void ibz_mod(ibz_t *r, const ibz_t *a, const ibz_t *b);
    unsigned long ibz_mod_ui(const ibz_t *n, unsigned long d);

    void ibz_pow(ibz_t *pow, const ibz_t *x, uint32_t e);
    /** Modular exponentiation.  The modulus must be nonzero and the exponent
     * nonnegative; violation of either precondition terminates the process. */
    void ibz_pow_mod(ibz_t *pow, const ibz_t *x, const ibz_t *e, const ibz_t *m);


    /** Greatest common divisor.  If the mathematical result is
     * 2^(IBZ_BITS-1), its unsigned magnitude bit pattern is returned. */
    void ibz_gcd(ibz_t *gcd, const ibz_t *a, const ibz_t *b);
    void ibz_gcdext(ibz_t *gcd, ibz_t *x, ibz_t *y, const ibz_t *a, const ibz_t *b);

    int ibz_invmod(ibz_t *inv, const ibz_t *a, const ibz_t *mod);

    int ibz_divides(const ibz_t *a, const ibz_t *b);

    int ibz_sqrt(ibz_t *sqrt, const ibz_t *a);

    void ibz_sqrt_floor(ibz_t *sqrt, const ibz_t *a);

    int ibz_legendre(const ibz_t *a, const ibz_t *p);

    int ibz_sqrt_mod_p(ibz_t *sqrt, const ibz_t *a, const ibz_t *p);

    extern const uint64_t ibz_const_zero[IBZ_LIMBS];
    extern const uint64_t ibz_const_one[IBZ_LIMBS];
    extern const uint64_t ibz_const_two[IBZ_LIMBS];
    extern const uint64_t ibz_const_three[IBZ_LIMBS];

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* IBZ_H */
