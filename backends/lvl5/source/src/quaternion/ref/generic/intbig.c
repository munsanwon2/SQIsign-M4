/*
 * Fixed-precision integer arithmetic implementation
 * Using 2's complement representation for signed integers
 * Uses the per-variant limb budgets selected in intbig.h.
 */

#include "intbig.h"
#if !defined(TARGET_ARM)
#include <stdlib.h>
#endif

static _Noreturn void
ibz_fatal_error(const char *category, const char *operation)
{
#if defined(TARGET_ARM)
    /* Legacy intbig arithmetic has void APIs for invariant violations.  Keep
     * the embedded fatal path freestanding; the HardFault handler may record
     * and reset the device. */
    (void)category;
    (void)operation;
    __builtin_trap();
    for (;;) {
    }
#else
    (void)category;
    (void)operation;
    abort();
#endif
}

#if defined(SQISIGN_INTBIG_OVERFLOW_CHECK)
static void
ibz_overflow_abort(const char *operation)
{
    ibz_fatal_error("OVERFLOW", operation);
}
#endif

static void
ibz_division_by_zero_abort(const char *operation)
{
    ibz_fatal_error("DIVZERO", operation);
}

static void
ibz_invalid_argument_abort(const char *operation)
{
    ibz_fatal_error("INVALID", operation);
}

const uint64_t ibz_const_zero[IBZ_LIMBS] = { 0 };
const uint64_t ibz_const_one[IBZ_LIMBS] = { 1 };
const uint64_t ibz_const_two[IBZ_LIMBS] = { 2 };
const uint64_t ibz_const_three[IBZ_LIMBS] = { 3 };

#define KARATSUBA_THRESHOLD 32

// Helper macros
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#if defined(__GNUC__) || defined(__clang__)
#define INTBIG_NOINLINE __attribute__((noinline))
#define INTBIG_ALWAYS_INLINE __attribute__((always_inline))
#else
#define INTBIG_NOINLINE
#define INTBIG_ALWAYS_INLINE
#endif

// Bit manipulation helpers
static inline int
clz64(uint64_t x)
{
    if (x == 0)
        return 64;
#if defined(TARGET_ARM) && (defined(__GNUC__) || defined(__clang__))
    const uint32_t high = (uint32_t)(x >> 32);
    return high != 0 ? __builtin_clz(high)
                     : 32 + __builtin_clz((uint32_t)x);
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_clzll(x);
#else
    int n = 0;
    if (x <= 0x00000000FFFFFFFFULL) {
        n += 32;
        x <<= 32;
    }
    if (x <= 0x0000FFFFFFFFFFFFULL) {
        n += 16;
        x <<= 16;
    }
    if (x <= 0x00FFFFFFFFFFFFFFULL) {
        n += 8;
        x <<= 8;
    }
    if (x <= 0x0FFFFFFFFFFFFFFFULL) {
        n += 4;
        x <<= 4;
    }
    if (x <= 0x3FFFFFFFFFFFFFFFULL) {
        n += 2;
        x <<= 2;
    }
    if (x <= 0x7FFFFFFFFFFFFFFFULL) {
        n += 1;
    }
    return n;
#endif
}

static inline int
ctz64(uint64_t x)
{
    if (x == 0)
        return 64;
#if defined(TARGET_ARM) && (defined(__GNUC__) || defined(__clang__))
    const uint32_t low = (uint32_t)x;
    return low != 0 ? __builtin_ctz(low)
                    : 32 + __builtin_ctz((uint32_t)(x >> 32));
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(x);
#else
    int n = 0;
    if ((x & 0x00000000FFFFFFFFULL) == 0) {
        n += 32;
        x >>= 32;
    }
    if ((x & 0x000000000000FFFFULL) == 0) {
        n += 16;
        x >>= 16;
    }
    if ((x & 0x00000000000000FFULL) == 0) {
        n += 8;
        x >>= 8;
    }
    if ((x & 0x000000000000000FULL) == 0) {
        n += 4;
        x >>= 4;
    }
    if ((x & 0x0000000000000003ULL) == 0) {
        n += 2;
        x >>= 2;
    }
    if ((x & 0x0000000000000001ULL) == 0) {
        n += 1;
    }
    return n;
#endif
}

// 64x64 -> 128 bit multiplication, with a portable 32-bit-half fallback
static void
mul64_128(uint64_t a, uint64_t b, uint64_t *hi, uint64_t *lo)
{
#if defined(HAVE_UINT128)
    __uint128_t product = (__uint128_t)a * b;
    *lo = (uint64_t)product;
    *hi = (uint64_t)(product >> 64);
#else
    uint32_t a_lo = (uint32_t)a;
    uint32_t a_hi = (uint32_t)(a >> 32);
    uint32_t b_lo = (uint32_t)b;
    uint32_t b_hi = (uint32_t)(b >> 32);

    uint64_t p00 = (uint64_t)a_lo * b_lo;
    uint64_t p01 = (uint64_t)a_lo * b_hi;
    uint64_t p10 = (uint64_t)a_hi * b_lo;
    uint64_t p11 = (uint64_t)a_hi * b_hi;

    uint64_t middle = p01 + p10;
    uint64_t carry = (middle < p01) ? (1ULL << 32) : 0;

    *lo = p00 + (middle << 32);
    carry += (*lo < p00) ? 1 : 0;
    *hi = p11 + (middle >> 32) + carry;
#endif
}

static size_t
ibz_used_limbs(const uint64_t *value, size_t capacity)
{
    while (capacity > 0 && value[capacity - 1] == 0)
        --capacity;
    return capacity;
}

/* Full unsigned product.  The output capacity is exactly twice the fixed
 * precision, which is also enough for two maximum-width magnitudes. */
static void
ibz_mul_unsigned_wide(uint64_t product[2 * IBZ_LIMBS],
                      const uint64_t *a,
                      size_t a_size,
                      const uint64_t *b,
                      size_t b_size)
{
    memset(product, 0, 2 * IBZ_LIMBS * sizeof(*product));

    for (size_t i = 0; i < a_size; ++i) {
        if (a[i] == 0)
            continue;

#if defined(HAVE_UINT128)
        uint64_t carry = 0;
        for (size_t j = 0; j < b_size; ++j) {
            const size_t index = i + j;
            __uint128_t accumulator = (__uint128_t)a[i] * b[j] +
                                      product[index] + carry;
            product[index] = (uint64_t)accumulator;
            carry = (uint64_t)(accumulator >> 64);
        }

        size_t index = i + b_size;
        while (carry != 0 && index < 2 * IBZ_LIMBS) {
            __uint128_t accumulator = (__uint128_t)product[index] + carry;
            product[index] = (uint64_t)accumulator;
            carry = (uint64_t)(accumulator >> 64);
            ++index;
        }
#else
        for (size_t j = 0; j < b_size; ++j) {
            if (b[j] == 0)
                continue;

            const size_t index = i + j;
            uint64_t hi, lo;
            mul64_128(a[i], b[j], &hi, &lo);

            uint64_t old = product[index];
            product[index] = old + lo;
            uint64_t carry = product[index] < old;

            size_t k = index + 1;
            if (k < 2 * IBZ_LIMBS) {
                old = product[k];
                product[k] = old + hi;
                uint64_t next_carry = product[k] < old;
                old = product[k];
                product[k] = old + carry;
                next_carry |= product[k] < old;
                carry = next_carry;
                ++k;

                while (carry != 0 && k < 2 * IBZ_LIMBS) {
                    old = product[k];
                    product[k] = old + 1;
                    carry = product[k] == 0;
                    ++k;
                }
            }
        }
#endif
    }
}

#if !defined(SQISIGN_INTBIG_OVERFLOW_CHECK)
/* Low half of an unsigned product.  Production arithmetic is defined at the
 * fixed precision and consumes only these limbs; computing and storing the
 * discarded upper triangle doubled both work and stack on Cortex-M4. */
static void
ibz_mul_unsigned_low(ibz_t *product,
                     const uint64_t *a,
                     size_t a_size,
                     const uint64_t *b,
                     size_t b_size)
{
    ibz_init(product);
    for (size_t i = 0; i < a_size; ++i) {
        if (a[i] == 0)
            continue;
        const size_t row_size = MIN(b_size, (size_t)IBZ_LIMBS - i);
        for (size_t j = 0; j < row_size; ++j) {
            if (b[j] == 0)
                continue;

            const size_t index = i + j;
            uint64_t hi, lo;
            mul64_128(a[i], b[j], &hi, &lo);

            uint64_t old = (*product)[index];
            (*product)[index] = old + lo;
            uint64_t carry = (*product)[index] < old;

            size_t k = index + 1;
            if (k < IBZ_LIMBS) {
                old = (*product)[k];
                (*product)[k] = old + hi;
                uint64_t next_carry = (*product)[k] < old;
                old = (*product)[k];
                (*product)[k] = old + carry;
                next_carry |= (*product)[k] < old;
                carry = next_carry;
                ++k;

                while (carry != 0 && k < IBZ_LIMBS) {
                    old = (*product)[k];
                    (*product)[k] = old + 1;
                    carry = (*product)[k] == 0;
                    ++k;
                }
            }
        }
    }
}
#endif

#if defined(HAVE_UINT128)
/* Knuth-style normalized long division in base 2^64.  The dividend may be a
 * regular fixed-width value or a full 2N-limb multiplication result.  Only
 * the remainder is required by wide modular multiplication, so quotient may
 * be NULL. */
static void
ibz_divrem_unsigned_wide(uint64_t quotient[2 * IBZ_LIMBS],
                         uint64_t remainder[IBZ_LIMBS],
                         const uint64_t *dividend,
                         size_t dividend_capacity,
                         const uint64_t divisor[IBZ_LIMBS])
{
    uint64_t u[2 * IBZ_LIMBS + 1];
    uint64_t v[IBZ_LIMBS];
    const size_t dividend_size = ibz_used_limbs(dividend, dividend_capacity);
    const size_t divisor_size = ibz_used_limbs(divisor, IBZ_LIMBS);

    if (divisor_size == 0)
        ibz_division_by_zero_abort("division");

    if (quotient != NULL)
        memset(quotient, 0, 2 * IBZ_LIMBS * sizeof(*quotient));
    memset(remainder, 0, IBZ_LIMBS * sizeof(*remainder));

    if (dividend_size == 0)
        return;
    if (dividend_size < divisor_size) {
        memcpy(remainder, dividend, dividend_size * sizeof(*remainder));
        return;
    }

    if (divisor_size == 1) {
        uint64_t carry = 0;
        for (size_t i = dividend_size; i-- > 0;) {
            __uint128_t numerator = ((__uint128_t)carry << 64) | dividend[i];
            if (quotient != NULL)
                quotient[i] = (uint64_t)(numerator / divisor[0]);
            carry = (uint64_t)(numerator % divisor[0]);
        }
        remainder[0] = carry;
        return;
    }

    memset(u, 0, sizeof(u));
    memset(v, 0, sizeof(v));

    const unsigned int normalization = (unsigned int)clz64(divisor[divisor_size - 1]);
    if (normalization == 0) {
        memcpy(v, divisor, divisor_size * sizeof(*v));
        memcpy(u, dividend, dividend_size * sizeof(*u));
    } else {
        uint64_t carry = 0;
        for (size_t i = 0; i < divisor_size; ++i) {
            const uint64_t limb = divisor[i];
            v[i] = (limb << normalization) | carry;
            carry = limb >> (64 - normalization);
        }

        carry = 0;
        for (size_t i = 0; i < dividend_size; ++i) {
            const uint64_t limb = dividend[i];
            u[i] = (limb << normalization) | carry;
            carry = limb >> (64 - normalization);
        }
        u[dividend_size] = carry;
    }

    const size_t quotient_size = dividend_size - divisor_size + 1;
    const uint64_t divisor_top = v[divisor_size - 1];
    for (size_t position = quotient_size; position-- > 0;) {
        const size_t top = position + divisor_size;
        uint64_t qhat;
        uint64_t rhat;
        int rhat_overflow = 0;

        if (u[top] >= divisor_top) {
            qhat = UINT64_MAX;
            __uint128_t sum = (__uint128_t)u[top - 1] + divisor_top;
            rhat = (uint64_t)sum;
            rhat_overflow = (int)(sum >> 64);
        } else {
            __uint128_t numerator = ((__uint128_t)u[top] << 64) | u[top - 1];
            qhat = (uint64_t)(numerator / divisor_top);
            rhat = (uint64_t)(numerator % divisor_top);
        }

        while (!rhat_overflow &&
               (__uint128_t)qhat * v[divisor_size - 2] >
                   ((__uint128_t)rhat << 64) + u[top - 2]) {
            --qhat;
            __uint128_t sum = (__uint128_t)rhat + divisor_top;
            rhat = (uint64_t)sum;
            rhat_overflow = (int)(sum >> 64);
        }

        uint64_t borrow = 0;
        for (size_t i = 0; i < divisor_size; ++i) {
            __uint128_t subtrahend = (__uint128_t)qhat * v[i] + borrow;
            const uint64_t low = (uint64_t)subtrahend;
            const uint64_t high = (uint64_t)(subtrahend >> 64);
            const uint64_t old = u[position + i];
            u[position + i] = old - low;
            borrow = high + (old < low);
        }

        const uint64_t old_top = u[top];
        u[top] = old_top - borrow;
        if (old_top < borrow) {
            --qhat;
            uint64_t carry = 0;
            for (size_t i = 0; i < divisor_size; ++i) {
                __uint128_t sum = (__uint128_t)u[position + i] + v[i] + carry;
                u[position + i] = (uint64_t)sum;
                carry = (uint64_t)(sum >> 64);
            }
            u[top] += carry;
        }

        if (quotient != NULL)
            quotient[position] = qhat;
    }

    if (normalization == 0) {
        memcpy(remainder, u, divisor_size * sizeof(*remainder));
    } else {
        for (size_t i = 0; i < divisor_size; ++i)
            remainder[i] = (u[i] >> normalization) |
                           (u[i + 1] << (64 - normalization));
    }
}
#endif

// Initialize/finalize
void
ibz_init(ibz_t *x)
{
    memset(*x, 0, sizeof(ibz_t));
}

void
ibz_finalize(ibz_t *x)
{
    /* Fixed-precision temporaries frequently contain secret ideal and
     * sampling state.  Volatile stores make finalization an actual wipe even
     * when the object becomes dead immediately afterwards. */
    volatile uint64_t *limbs = (volatile uint64_t *)(void *)*x;
    for (size_t i = 0; i < IBZ_LIMBS; ++i)
        limbs[i] = 0;
}

// Copy and swap
void
ibz_copy(ibz_t *target, const ibz_t *value)
{
    memcpy(*target, *value, sizeof(ibz_t));
}

void
ibz_swap(ibz_t *a, ibz_t *b)
{
    for (int i = 0; i < IBZ_LIMBS; ++i) {
        const uint64_t tmp = (*a)[i];
        (*a)[i] = (*b)[i];
        (*b)[i] = tmp;
    }
}

// Check if negative (2's complement)
int
ibz_is_negative(const ibz_t *x)
{
    return ((*x)[IBZ_LIMBS - 1] >> 63) & 1;
}

// Raw two's-complement negation. Multiplication also uses this on an unsigned
// magnitude before assigning its sign, so overflow policy belongs in the
// public signed wrapper below.
static void
ibz_neg_raw(ibz_t *neg, const ibz_t *a)
{
    uint64_t carry = 1;
    for (int i = 0; i < IBZ_LIMBS; i++) {
        uint64_t tmp = ~(*a)[i];
        (*neg)[i] = tmp + carry;
        carry = ((*neg)[i] < tmp) ? 1 : 0;
    }
}

/*
 * Several internal algorithms need an unsigned magnitude.  In particular,
 * the magnitude of the most-negative signed value is 2^(IBZ_BITS-1), which
 * has no positive ibz_t representation.  Keep that value as an unsigned limb
 * array instead of feeding it to the signed comparison/arithmetic helpers.
 */
static void
ibz_abs_unsigned(ibz_t *magnitude, const ibz_t *a)
{
    if (ibz_is_negative(a))
        ibz_neg_raw(magnitude, a);
    else
        ibz_copy(magnitude, a);
}

static int
ibz_cmp_unsigned(const ibz_t *a, const ibz_t *b)
{
    for (int i = IBZ_LIMBS - 1; i >= 0; --i) {
        if ((*a)[i] > (*b)[i])
            return 1;
        if ((*a)[i] < (*b)[i])
            return -1;
    }
    return 0;
}

static int
ibz_is_min_value(const ibz_t *a)
{
    if ((*a)[IBZ_LIMBS - 1] != (UINT64_C(1) << 63))
        return 0;
    for (int i = 0; i < IBZ_LIMBS - 1; ++i)
        if ((*a)[i] != 0)
            return 0;
    return 1;
}

static int
ibz_bitsize_unsigned(const ibz_t *a)
{
    for (int i = IBZ_LIMBS - 1; i >= 0; --i) {
        if ((*a)[i] != 0)
            return i * 64 + (64 - clz64((*a)[i]));
    }
    return 0;
}

static void
ibz_sub_unsigned(ibz_t *difference, const ibz_t *a, const ibz_t *b)
{
    uint64_t borrow = 0;
    for (int i = 0; i < IBZ_LIMBS; ++i) {
        uint64_t ai = (*a)[i];
        uint64_t first = ai - borrow;
        uint64_t first_borrow = first > ai;
        uint64_t second = first - (*b)[i];
        uint64_t second_borrow = second > first;
        (*difference)[i] = second;
        borrow = first_borrow | second_borrow;
    }
}

static void
ibz_shift_right_unsigned(ibz_t *result, const ibz_t *a, uint32_t shift)
{
    if (shift >= IBZ_BITS) {
        ibz_init(result);
        return;
    }

    const uint32_t limb_shift = shift / 64;
    const uint32_t bit_shift = shift % 64;
    const uint32_t output_limbs = IBZ_LIMBS - limb_shift;

    /* Low-to-high publication leaves every not-yet-read source limb intact,
     * including when result == a.  Avoiding the former full-width input copy
     * matters especially inside portable division, which shifts once per
     * quotient bit. */
    for (uint32_t i = 0; i < output_limbs; ++i) {
        uint64_t value = (*a)[i + limb_shift] >> bit_shift;
        if (bit_shift != 0 && i + limb_shift + 1 < IBZ_LIMBS)
            value |= (*a)[i + limb_shift + 1] << (64 - bit_shift);
        (*result)[i] = value;
    }
    for (uint32_t i = output_limbs; i < IBZ_LIMBS; ++i)
        (*result)[i] = 0;
}

static void
ibz_shift_left_unsigned(ibz_t *result, const ibz_t *a, uint32_t shift)
{
    if (shift >= IBZ_BITS) {
        ibz_init(result);
        return;
    }

    const uint32_t limb_shift = shift / 64;
    const uint32_t bit_shift = shift % 64;

    /* High-to-low is the alias-safe direction for a left shift. */
    for (int i = IBZ_LIMBS - 1; i >= (int)limb_shift; --i) {
        uint64_t value = (*a)[i - limb_shift] << bit_shift;
        if (bit_shift != 0 && i > (int)limb_shift)
            value |= (*a)[i - limb_shift - 1] >> (64 - bit_shift);
        (*result)[i] = value;
    }
    for (uint32_t i = 0; i < limb_shift; ++i)
        (*result)[i] = 0;
}

#if !defined(HAVE_UINT128)
/* Return a base-2^32 digit without aliasing ibz_t through a uint32_t pointer.
 * Splitting the numeric uint64_t limb also makes this independent of host
 * byte order. */
static inline INTBIG_ALWAYS_INLINE uint32_t
ibz_word32(const ibz_t *value, size_t word)
{
    const uint64_t limb = (*value)[word >> 1];
    return (word & 1u) ? (uint32_t)(limb >> 32) : (uint32_t)limb;
}

/* GCC's may_alias type makes this a strict-aliasing-safe aligned 32-bit store
 * on little-endian Cortex-M4.  The numeric fallback is byte-order independent
 * and avoids relying on a freestanding compiler to inline four-byte memcpy. */
#if !defined(TARGET_BIG_ENDIAN) && \
    (defined(__GNUC__) || defined(__clang__))
typedef uint32_t ibz_alias_u32_t __attribute__((__may_alias__));
#endif

static inline INTBIG_ALWAYS_INLINE void
ibz_set_word32(ibz_t *value, size_t word, uint32_t digit)
{
#if !defined(TARGET_BIG_ENDIAN) && \
    (defined(__GNUC__) || defined(__clang__))
    ((ibz_alias_u32_t *)(void *)*value)[word] = digit;
#else
    const size_t limb = word >> 1;
    if ((word & 1u) != 0)
        (*value)[limb] = ((*value)[limb] & UINT32_MAX) |
                         ((uint64_t)digit << 32);
    else
        (*value)[limb] = ((*value)[limb] & (UINT64_MAX << 32)) | digit;
#endif
}

static size_t
ibz_used_words32(const ibz_t *value)
{
    size_t words = 2u * (size_t)IBZ_LIMBS;
    while (words != 0 && ibz_word32(value, words - 1) == 0)
        --words;
    return words;
}
#endif

static uint32_t ibz_div_unsigned_small(ibz_t *quotient,
                                       const ibz_t *dividend,
                                       uint32_t divisor);
static uint32_t ibz_mod_unsigned_small(const ibz_t *dividend,
                                       uint32_t divisor);

/* Unsigned fixed-width division.  Targets with 128-bit intermediates use
 * normalized base-2^64 long division.  The portable path uses the same
 * algorithm in base 2^32, where every digit product fits uint64_t.  Both
 * paths publish outputs only after consuming inputs, so every output/input
 * aliasing pattern remains valid. */
static void
ibz_div_unsigned(ibz_t *quotient,
                 ibz_t *remainder,
                 const ibz_t *dividend,
                 const ibz_t *divisor)
{
    const size_t divisor_used = ibz_used_limbs(*divisor, IBZ_LIMBS);
    if (divisor_used == 0)
        ibz_division_by_zero_abort("division");

    /* Exact powers of two occur throughout the Qlapoti/id2iso path.  A
     * logical shift and low-bit mask avoid both the base-2^16 UDIV loop for
     * tiny powers and the bitwise long divider for wide powers.  Preserve a
     * requested remainder before publishing an aliasing quotient. */
    const size_t power_limb = divisor_used - 1;
    const uint64_t power_word = (*divisor)[power_limb];
    int divisor_is_power_of_two =
        (power_word & (power_word - UINT64_C(1))) == 0;
    for (size_t i = 0; i < power_limb && divisor_is_power_of_two; ++i)
        divisor_is_power_of_two = (*divisor)[i] == 0;
    if (divisor_is_power_of_two) {
        const uint32_t shift =
            (uint32_t)(power_limb * 64 + (63 - clz64(power_word)));
        ibz_t power_remainder;
        if (remainder != NULL) {
            ibz_copy(&power_remainder, dividend);
            const uint32_t whole_limbs = shift / 64;
            const uint32_t remaining_bits = shift % 64;
            if (remaining_bits == 0) {
                for (uint32_t i = whole_limbs; i < IBZ_LIMBS; ++i)
                    power_remainder[i] = 0;
            } else {
                power_remainder[whole_limbs] &=
                    (UINT64_C(1) << remaining_bits) - UINT64_C(1);
                for (uint32_t i = whole_limbs + 1; i < IBZ_LIMBS; ++i)
                    power_remainder[i] = 0;
            }
        }
        if (quotient != NULL)
            ibz_shift_right_unsigned(quotient, dividend, shift);
        if (remainder != NULL)
            ibz_copy(remainder, &power_remainder);
        return;
    }

    /* Qlapoti's Cornacchia path divides by its 101-entry small-prime table.
     * Use two 32-bit chunks per 64-bit limb instead of running one full-width
     * compare/subtract iteration per quotient bit. */
    if (divisor_used == 1 && (*divisor)[0] <= UINT32_MAX) {
        const uint32_t small_divisor = (uint32_t)(*divisor)[0];
        const uint32_t small_remainder =
            quotient != NULL
                ? ibz_div_unsigned_small(quotient, dividend, small_divisor)
                : ibz_mod_unsigned_small(dividend, small_divisor);
        if (remainder != NULL) {
            ibz_init(remainder);
            (*remainder)[0] = small_remainder;
        }
        return;
    }

#if defined(HAVE_UINT128)
    uint64_t q_wide[2 * IBZ_LIMBS];
    uint64_t r_words[IBZ_LIMBS];
    ibz_divrem_unsigned_wide(quotient != NULL ? q_wide : NULL,
                             r_words,
                             *dividend,
                             IBZ_LIMBS,
                             *divisor);
    if (quotient != NULL)
        memcpy(*quotient, q_wide, sizeof(ibz_t));
    if (remainder != NULL)
        memcpy(*remainder, r_words, sizeof(ibz_t));
#else
    /* Knuth division in B=2^32.  The normalized divisor occupies one
     * full-width local object.  The caller-provided remainder is also the
     * normalized dividend workspace; all current callers already own it, so
     * this removes one full-width temporary relative to the old fallback.
     * One possible extra dividend digit is kept in a scalar. */
    uint32_t v[2 * IBZ_LIMBS];
    const size_t dividend_size = ibz_used_words32(dividend);
    const size_t divisor_size = ibz_used_words32(divisor);
    ibz_t *quotient_output = quotient == remainder ? NULL : quotient;

    if (remainder == NULL)
        ibz_invalid_argument_abort("unsigned division remainder workspace");

    if (dividend_size < divisor_size) {
        /* Copy before a quotient output can overwrite either input. */
        if (remainder != dividend)
            ibz_copy(remainder, dividend);
        if (quotient_output != NULL)
            ibz_init(quotient_output);
        return;
    }

    const uint32_t divisor_top = ibz_word32(divisor, divisor_size - 1);
    const unsigned int normalization =
        (unsigned int)(clz64((uint64_t)divisor_top) - 32);
    uint32_t carry = 0;
    for (size_t i = 0; i < divisor_size; ++i) {
        const uint64_t value =
            ((uint64_t)ibz_word32(divisor, i) << normalization) | carry;
        v[i] = (uint32_t)value;
        carry = (uint32_t)(value >> 32);
    }

    /* v now preserves a divisor that aliases the remainder. */
    if (remainder != dividend)
        ibz_copy(remainder, dividend);
    carry = 0;
    for (size_t i = 0; i < dividend_size; ++i) {
        const uint64_t value =
            ((uint64_t)ibz_word32(remainder, i) << normalization) | carry;
        ibz_set_word32(remainder, i, (uint32_t)value);
        carry = (uint32_t)(value >> 32);
    }
    uint32_t u_top = carry;

    /* Inputs are fully consumed now, so a quotient may safely alias either. */
    if (quotient_output != NULL)
        ibz_init(quotient_output);

    const size_t quotient_size = dividend_size - divisor_size + 1;
    const uint32_t normalized_divisor_top = v[divisor_size - 1];
    for (size_t position = quotient_size; position-- > 0;) {
        const size_t top = position + divisor_size;
        uint32_t top_word =
            top == dividend_size ? u_top : ibz_word32(remainder, top);
        uint32_t qhat;
        uint64_t rhat;

        if (top_word >= normalized_divisor_top) {
            qhat = UINT32_MAX;
            rhat = (uint64_t)ibz_word32(remainder, top - 1) +
                   normalized_divisor_top;
        } else {
            const uint64_t numerator =
                ((uint64_t)top_word << 32) |
                ibz_word32(remainder, top - 1);
            qhat = (uint32_t)(numerator / normalized_divisor_top);
            rhat = numerator % normalized_divisor_top;
        }

        while (rhat <= UINT32_MAX &&
               (uint64_t)qhat * v[divisor_size - 2] >
                   (rhat << 32) + ibz_word32(remainder, top - 2)) {
            --qhat;
            rhat += normalized_divisor_top;
        }

        uint64_t borrow = 0;
        for (size_t i = 0; i < divisor_size; ++i) {
            const uint64_t product = (uint64_t)qhat * v[i] + borrow;
            const uint32_t product_low = (uint32_t)product;
            const uint32_t old = ibz_word32(remainder, position + i);
            ibz_set_word32(remainder, position + i, old - product_low);
            borrow = (product >> 32) + (old < product_low);
        }

        const int underflow = (uint64_t)top_word < borrow;
        top_word -= (uint32_t)borrow;
        if (underflow) {
            --qhat;
            uint64_t add_carry = 0;
            for (size_t i = 0; i < divisor_size; ++i) {
                const uint64_t sum = (uint64_t)ibz_word32(
                                         remainder, position + i) +
                                     v[i] + add_carry;
                ibz_set_word32(remainder, position + i, (uint32_t)sum);
                add_carry = sum >> 32;
            }
            top_word += (uint32_t)add_carry;
        }

        if (top == dividend_size)
            u_top = top_word;
        else
            ibz_set_word32(remainder, top, top_word);
        if (quotient_output != NULL)
            (*quotient_output)[position >> 1] |=
                (uint64_t)qhat << ((position & 1u) * 32u);
    }

    if (normalization != 0) {
        for (size_t i = 0; i < divisor_size; ++i) {
            const uint32_t word = ibz_word32(remainder, i);
            const uint32_t next = i + 1 == dividend_size
                                      ? u_top
                                      : ibz_word32(remainder, i + 1);
            ibz_set_word32(remainder,
                           i,
                           (word >> normalization) |
                               (next << (32 - normalization)));
        }
    }
    for (size_t i = divisor_size; i < 2u * (size_t)IBZ_LIMBS; ++i)
        ibz_set_word32(remainder, i, 0);
#endif
}

static uint32_t
ibz_div_unsigned_small(ibz_t *quotient, const ibz_t *dividend, uint32_t divisor)
{
    if (divisor <= UINT16_MAX) {
        uint32_t remainder = 0;
        for (int i = IBZ_LIMBS - 1; i >= 0; --i) {
            const uint64_t word = (*dividend)[i];
            uint64_t quotient_word = 0;
            for (int chunk = 3; chunk >= 0; --chunk) {
                const uint32_t numerator =
                    (remainder << 16) |
                    (uint32_t)((word >> (16 * chunk)) & UINT64_C(0xffff));
                const uint32_t quotient_chunk = numerator / divisor;
                remainder = numerator % divisor;
                quotient_word |= (uint64_t)quotient_chunk << (16 * chunk);
            }
            (*quotient)[i] = quotient_word;
        }
        return remainder;
    }

    uint64_t remainder = 0;

    for (int i = IBZ_LIMBS - 1; i >= 0; --i) {
        uint64_t hi, lo;
        /* remainder < divisor <= UINT32_MAX, so each shifted half fits u64. */
        hi = (remainder << 32) | ((*dividend)[i] >> 32);
        const uint64_t high_quotient = hi / divisor;
        remainder = hi % divisor;
        lo = (remainder << 32) | ((*dividend)[i] & UINT64_C(0xffffffff));
        (*quotient)[i] = (high_quotient << 32) | (lo / divisor);
        remainder = lo % divisor;
    }

    return (uint32_t)remainder;
}

static uint32_t
ibz_mod_unsigned_small(const ibz_t *dividend, uint32_t divisor)
{
    if (divisor <= UINT16_MAX) {
        uint32_t remainder = 0;
        size_t size = ibz_used_limbs(*dividend, IBZ_LIMBS);
        while (size-- > 0) {
            const uint64_t word = (*dividend)[size];
            for (int chunk = 3; chunk >= 0; --chunk) {
                const uint32_t numerator =
                    (remainder << 16) |
                    (uint32_t)((word >> (16 * chunk)) & UINT64_C(0xffff));
                remainder = numerator % divisor;
            }
        }
        return remainder;
    }

    uint64_t remainder = 0;
    size_t size = ibz_used_limbs(*dividend, IBZ_LIMBS);

    while (size-- > 0) {
        uint64_t half = (remainder << 32) | ((*dividend)[size] >> 32);
        remainder = half % divisor;
        half = (remainder << 32) |
               ((*dividend)[size] & UINT64_C(0xffffffff));
        remainder = half % divisor;
    }
    return (uint32_t)remainder;
}

// Negation (2's complement)
void
ibz_neg(ibz_t *neg, const ibz_t *a)
{
#if defined(SQISIGN_INTBIG_OVERFLOW_CHECK)
    int input_negative = ibz_is_negative(a);
#endif
    ibz_neg_raw(neg, a);
#if defined(SQISIGN_INTBIG_OVERFLOW_CHECK)
    /* The only negative input whose two's-complement negation is still
     * negative is -2^(IBZ_BITS-1), which has no positive counterpart. */
    if (input_negative && ibz_is_negative(neg))
        ibz_overflow_abort("negation");
#endif
}

// Absolute value
void
ibz_abs(ibz_t *abs, const ibz_t *a)
{
    /* The magnitude of INTBIG_MIN is not representable as a positive ibz_t.
     * Preserve its magnitude bit pattern (and therefore INTBIG_MIN itself),
     * just as fixed-width two's-complement absolute-value instructions do. */
    ibz_abs_unsigned(abs, a);
}

// Addition
void
ibz_add(ibz_t *sum, const ibz_t *a, const ibz_t *b)
{
#if defined(SQISIGN_INTBIG_OVERFLOW_CHECK)
    int a_negative = ibz_is_negative(a);
    int b_negative = ibz_is_negative(b);
#endif
    uint64_t carry = 0;
    for (int i = 0; i < IBZ_LIMBS; i++) {
        uint64_t tmp = (*a)[i] + carry;
        carry = (tmp < (*a)[i]) ? 1 : 0;
        (*sum)[i] = tmp + (*b)[i];
        carry += ((*sum)[i] < tmp) ? 1 : 0;
    }
#if defined(SQISIGN_INTBIG_OVERFLOW_CHECK)
    if (a_negative == b_negative && ibz_is_negative(sum) != a_negative)
        ibz_overflow_abort("addition");
#endif
}

// Subtraction
void
ibz_sub(ibz_t *diff, const ibz_t *a, const ibz_t *b)
{
#if defined(SQISIGN_INTBIG_OVERFLOW_CHECK)
    int a_negative = ibz_is_negative(a);
    int b_negative = ibz_is_negative(b);
#endif
    uint64_t borrow = 0;
    for (int i = 0; i < IBZ_LIMBS; i++) {
        uint64_t tmp = (*a)[i] - borrow;
        borrow = (tmp > (*a)[i]) ? 1 : 0;
        uint64_t tmp2 = tmp - (*b)[i];
        borrow += (tmp2 > tmp) ? 1 : 0;
        (*diff)[i] = tmp2;
    }
#if defined(SQISIGN_INTBIG_OVERFLOW_CHECK)
    if (a_negative != b_negative && ibz_is_negative(diff) != a_negative)
        ibz_overflow_abort("subtraction");
#endif
}

// static void
// ibz_mul_karatsuba_internal(uint64_t *result,
//                            const uint64_t *a,
//                            int a_size,
//                            const uint64_t *b,
//                            int b_size,
//                            uint64_t *workspace)
// {
//     // Base case: use schoolbook for small sizes
//     if (a_size < KARATSUBA_THRESHOLD || b_size < KARATSUBA_THRESHOLD) {
//         // Schoolbook multiplication
//         memset(result, 0, (a_size + b_size) * sizeof(uint64_t));

//         for (int i = 0; i < a_size; i++) {
//             if (a[i] == 0)
//                 continue;

//             uint64_t carry = 0;
//             for (int j = 0; j < b_size; j++) {
//                 uint64_t hi, lo;
//                 mul64_128(a[i], b[j], &hi, &lo);

//                 // Add lo to result[i+j]
//                 result[i + j] += lo;
//                 uint64_t c1 = (result[i + j] < lo) ? 1 : 0;

//                 // Add previous carry
//                 result[i + j] += carry;
//                 c1 += (result[i + j] < carry) ? 1 : 0;

//                 // New carry = hi + c1
//                 carry = hi + c1;
//             }
//             if (carry && (i + b_size < a_size + b_size)) {
//                 result[i + b_size] = carry;
//             }
//         }
//         return;
//     }

//     // Karatsuba: split numbers in half
//     int m = MAX(a_size, b_size) / 2;

//     // a = a1*B^m + a0
//     // b = b1*B^m + b0
//     // result = a1*b1*B^(2m) + ((a1+a0)*(b1+b0) - a1*b1 - a0*b0)*B^m + a0*b0

//     int a0_size = MIN(m, a_size);
//     int a1_size = a_size > m ? a_size - m : 0;
//     int b0_size = MIN(m, b_size);
//     int b1_size = b_size > m ? b_size - m : 0;

//     // Recursive calls for a0*b0 and a1*b1
//     ibz_mul_karatsuba_internal(workspace, a, a0_size, b, b0_size, workspace + 2 * (a_size + b_size));

//     if (a1_size > 0 && b1_size > 0) {
//         ibz_mul_karatsuba_internal(
//             workspace + 2 * m, a + m, a1_size, b + m, b1_size, workspace + 2 * (a_size + b_size));
//     }

//     // Copy results to final location
//     memcpy(result, workspace, (a_size + b_size) * sizeof(uint64_t));
// }

// // OPTIMIZED: Multiplication with Karatsuba for large numbers
// void
// ibz_mul(ibz_t *prod, const ibz_t *a, const ibz_t *b)
// {
//     // Find actual sizes
//     int a_size = IBZ_LIMBS;
//     int b_size = IBZ_LIMBS;

//     while (a_size > 1 && (*a)[a_size - 1] == 0)
//         a_size--;
//     while (b_size > 1 && (*b)[b_size - 1] == 0)
//         b_size--;

//     // // OPTIMIZATION: Use Karatsuba for large multiplications
//     // if (a_size > KARATSUBA_THRESHOLD && b_size > KARATSUBA_THRESHOLD) {
//     //     static uint64_t workspace[6 * IBZ_LIMBS]; // Working space for Karatsuba
//     //     static uint64_t temp_result[2 * IBZ_LIMBS];

//     //     memset(temp_result, 0, sizeof(temp_result));
//     //     ibz_mul_karatsuba_internal(temp_result, *a, a_size, *b, b_size, workspace);

//     //     // Copy result
//     //     for (int i = 0; i < IBZ_LIMBS; i++) {
//     //         (*prod)[i] = temp_result[i];
//     //     }
//     //     return;
//     // }

//     // Fallback to schoolbook for smaller numbers
//     static uint64_t temp_result[2 * IBZ_LIMBS];
//     memset(temp_result, 0, sizeof(temp_result));

//     for (int i = 0; i < a_size; i++) {
//         uint64_t a_limb = (*a)[i];
//         if (a_limb == 0)
//             continue;

//         for (int j = 0; j < b_size; j++) {
//             uint64_t b_limb = (*b)[j];
//             if (b_limb == 0)
//                 continue;

//             if (i + j >= 2 * IBZ_LIMBS)
//                 continue;

//             uint64_t hi, lo;
//             mul64_128(a_limb, b_limb, &hi, &lo);

//             int k = i + j;
//             temp_result[k] += lo;
//             uint64_t carry = (temp_result[k] < lo) ? 1 : 0;

//             k++;
//             if (k < 2 * IBZ_LIMBS) {
//                 temp_result[k] += hi;
//                 uint64_t c2 = (temp_result[k] < hi) ? 1 : 0;
//                 temp_result[k] += carry;
//                 c2 += (temp_result[k] < carry) ? 1 : 0;
//                 carry = c2;

//                 k++;
//                 while (carry && k < 2 * IBZ_LIMBS) {
//                     temp_result[k] += carry;
//                     carry = (temp_result[k] == 0) ? 1 : 0;
//                     k++;
//                 }
//             }
//         }
//     }

//     for (int i = 0; i < IBZ_LIMBS; i++) {
//         (*prod)[i] = temp_result[i];
//     }
// }

void
ibz_mul(ibz_t *prod, const ibz_t *a, const ibz_t *b)
{
    /* aa and bb consume both inputs before prod is published, so they are
     * also the alias copies when prod == a or prod == b. */
    int a_neg = ibz_is_negative(a);
    int b_neg = ibz_is_negative(b);
    int neg = (a_neg != b_neg);

    ibz_t aa, bb;
    ibz_abs(&aa, a);
    ibz_abs(&bb, b);

    // Find actual sizes on |a|, |b|
    size_t a_size = IBZ_LIMBS;
    size_t b_size = IBZ_LIMBS;

    while (a_size > 1 && aa[a_size - 1] == 0) a_size--;
    while (b_size > 1 && bb[b_size - 1] == 0) b_size--;

#if defined(SQISIGN_INTBIG_OVERFLOW_CHECK)
    /* Checked builds need the high half to decide whether the mathematical
     * signed result fits. */
    uint64_t temp_result[2 * IBZ_LIMBS];
    ibz_mul_unsigned_wide(temp_result, aa, a_size, bb, b_size);

    {
        int overflow = 0;
        for (int i = IBZ_LIMBS; i < 2 * IBZ_LIMBS; i++)
            overflow |= temp_result[i] != 0;

        if (!neg) {
            overflow |= (temp_result[IBZ_LIMBS - 1] >> 63) != 0;
        } else if (temp_result[IBZ_LIMBS - 1] > (UINT64_C(1) << 63)) {
            overflow = 1;
        } else if (temp_result[IBZ_LIMBS - 1] == (UINT64_C(1) << 63)) {
            for (int i = 0; i < IBZ_LIMBS - 1; i++)
                overflow |= temp_result[i] != 0;
        }

        if (overflow)
            ibz_overflow_abort("multiplication");
    }

    // write low limbs
    for (int i = 0; i < IBZ_LIMBS; i++) {
        (*prod)[i] = temp_result[i];
    }
#else
    ibz_mul_unsigned_low(prod, aa, a_size, bb, b_size);
#endif

    // apply sign
    if (neg && !ibz_is_zero(prod)) {
        ibz_neg_raw(prod, prod);
    }
}

// void
// ibz_mul_2exp(ibz_t* result, const ibz_t* a, size_t shift)
// {
//     if (shift == 0) {
//         ibz_copy(result, a);
//         return;
//     }

//     size_t limb_shift = shift / 64;
//     size_t bit_shift = shift % 64;

//     ibz_init(result);

//     if (bit_shift == 0) {
//         // limb ���� ����Ʈ��
//         for (size_t i = 0; i < IBZ_LIMBS - limb_shift; i++) {
//             (*result)[i + limb_shift] = (*a)[i];
//         }
//     } else {
//         // limb + bit ����Ʈ
//         for (size_t i = 0; i < IBZ_LIMBS - limb_shift; i++) {
//             (*result)[i + limb_shift] |= (*a)[i] << bit_shift;
//             if (i + limb_shift + 1 < IBZ_LIMBS) {
//                 (*result)[i + limb_shift + 1] |= (*a)[i] >> (64 - bit_shift);
//             }
//         }
//     }
// }
void
ibz_mul_2exp(ibz_t *result, const ibz_t *a, size_t shift)
{
    /* The private left-shift helper publishes high-to-low, hence it is
     * already safe when result == a.  Avoid reserving a full-width alias
     * copy in every production call. */
#if defined(SQISIGN_INTBIG_OVERFLOW_CHECK)
    ibz_t a_copy;
    const ibz_t *src = a;
    if (result == a) {
        ibz_copy(&a_copy, a);
        src = &a_copy;
    }
#else
    const ibz_t *src = a;
#endif

    if (shift >= IBZ_BITS) {
#if defined(SQISIGN_INTBIG_OVERFLOW_CHECK)
        if (!ibz_is_zero(src))
            ibz_overflow_abort("left shift");
#endif
        ibz_init(result);
        return;
    }

    ibz_shift_left_unsigned(result, src, (uint32_t)shift);

#if defined(SQISIGN_INTBIG_OVERFLOW_CHECK)
    /* A signed left shift is exact iff shifting the result arithmetically
     * back by the same amount recovers the input. */
    {
        ibz_t roundtrip;
        ibz_init(&roundtrip);
        ibz_div_2exp(&roundtrip, result, (uint32_t)shift);
        if (ibz_cmp(&roundtrip, src) != 0)
            ibz_overflow_abort("left shift");
    }
#endif
}

// Division by power of 2
void
ibz_div_2exp(ibz_t *quotient, const ibz_t *a, uint32_t exp)
{
    const int negative = ibz_is_negative(a);

    /* Match truncating integer division, not an arithmetic right shift:
     * -3 / 2 is -1 and every finite value divided by 2^IBZ_BITS is zero. */
    if (!negative) {
        ibz_shift_right_unsigned(quotient, a, exp);
        return;
    }

    if (exp >= IBZ_BITS) {
        ibz_init(quotient);
        return;
    }

    /* Arithmetic right shift rounds a negative value down.  Add one exactly
     * when discarded low bits are nonzero to obtain truncation toward zero.
     * Determine that condition before an in-place publication. */
    int discarded = 0;
    const uint32_t discarded_limbs = exp / 64;
    const uint32_t discarded_bits = exp % 64;
    for (uint32_t limb = 0; limb < discarded_limbs; ++limb)
        discarded |= (*a)[limb] != 0;
    if (discarded_bits != 0) {
        const uint64_t mask =
            (UINT64_C(1) << discarded_bits) - UINT64_C(1);
        discarded |= ((*a)[discarded_limbs] & mask) != 0;
    }

    ibz_shift_right_unsigned(quotient, a, exp);
    const uint32_t sign_limbs = exp / 64;
    const uint32_t sign_bits = exp % 64;
    for (uint32_t limb = IBZ_LIMBS - sign_limbs;
         limb < IBZ_LIMBS;
         ++limb)
        (*quotient)[limb] = UINT64_MAX;
    if (sign_bits != 0) {
        const uint32_t limb = IBZ_LIMBS - sign_limbs - 1;
        (*quotient)[limb] |= UINT64_MAX << (64 - sign_bits);
    }
    if (discarded)
        ibz_add(quotient, quotient, (const ibz_t *)&ibz_const_one);
}

// Comparison
int
ibz_cmp(const ibz_t *a, const ibz_t *b)
{
    int a_neg = ibz_is_negative(a);
    int b_neg = ibz_is_negative(b);

    if (a_neg != b_neg) {
        return a_neg ? -1 : 1;
    }

    for (int i = IBZ_LIMBS - 1; i >= 0; i--) {
        if ((*a)[i] > (*b)[i])
            return 1;
        if ((*a)[i] < (*b)[i])
            return -1;
    }
    return 0;
}

// Basic predicates
int
ibz_is_zero(const ibz_t *x)
{
    for (int i = 0; i < IBZ_LIMBS; i++) {
        if ((*x)[i] != 0)
            return 0;
    }
    return 1;
}

int
ibz_is_one(const ibz_t *x)
{
    if ((*x)[0] != 1)
        return 0;
    for (int i = 1; i < IBZ_LIMBS; i++) {
        if ((*x)[i] != 0)
            return 0;
    }
    return 1;
}

int
ibz_is_even(const ibz_t *x)
{
    return ((*x)[0] & 1) == 0;
}

int
ibz_is_odd(const ibz_t *x)
{
    return ((*x)[0] & 1) == 1;
}

// Set from int32
void
ibz_set(ibz_t *i, int32_t x)
{
    if (x >= 0) {
        (*i)[0] = (uint64_t)x;
        memset(&(*i)[1], 0, (IBZ_LIMBS - 1) * sizeof(uint64_t));
    } else {
        (*i)[0] = (uint64_t)x;
        memset(&(*i)[1], 0xFF, (IBZ_LIMBS - 1) * sizeof(uint64_t));
    }
}

void
ibz_set_u64(ibz_t *i, uint64_t x)
{
    (*i)[0] = x;
    memset(&(*i)[1], 0, (IBZ_LIMBS - 1) * sizeof(uint64_t));
}

// int
// ibz_convert_to_str(const ibz_t *i, char *str, int base)
// {
//     if (!str || (base != 10 && base != 16))
//         return 0;

//     ibz_t abs_i, base_ibz, q, r;
//     ibz_abs(&abs_i, i);
//     ibz_set(&base_ibz, base);

//     char temp[4096];
//     int pos = 0;

//     if (ibz_is_zero(i)) {
//         str[0] = '0';
//         str[1] = '\0';
//         return 1;
//     }

//     ibz_copy(&q, &abs_i);
//     while (!ibz_is_zero(&q)) {
//         ibz_div(&q, &r, &q, &base_ibz);
//         int digit = (int)r[0];
//         temp[pos++] = (digit < 10) ? ('0' + digit) : ('a' + digit - 10);
//     }

//     int offset = 0;
//     if (ibz_is_negative(i)) {
//         str[offset++] = '-';
//     }

//     for (int ind = 0; ind < pos; ind++) {
//         str[offset + ind] = temp[pos - 1 - ind];
//     }
//     str[offset + pos] = '\0';

//     return 1;
// }
int
ibz_convert_to_str(const ibz_t *i, char *str, int base)
{
#if defined(TARGET_ARM)
    /* Integer formatting is a host diagnostic.  Avoid its IBZ_BITS-byte
     * scratch frame and all stdio dependencies in the Cortex-M4 build. */
    (void)i;
    (void)str;
    (void)base;
    return 0;
#else
    if (!str || (base != 10 && base != 16))
        return 0;

    ibz_t magnitude, q;
    ibz_abs_unsigned(&magnitude, i);

    /* Base 10 is the longest supported representation and always uses fewer
     * than IBZ_BITS digits.  Unlike the former fixed 4096-byte buffer, this
     * remains correct for every configured IBZ_LIMBS value. */
    char temp[IBZ_BITS + 1];
    size_t pos = 0;

    if (ibz_is_zero(i)) {
        str[0] = '0';
        str[1] = '\0';
        return 1;
    }

    ibz_copy(&q, &magnitude);
    while (!ibz_is_zero(&q)) {
        ibz_t q_next;
        int digit = (int)ibz_div_unsigned_small(&q_next, &q, (uint32_t)base);

        temp[pos++] = (digit < 10)
                          ? ('0' + digit)
                          : ('a' + digit - 10);

        ibz_copy(&q, &q_next);
    }

    size_t offset = 0;
    if (ibz_is_negative(i)) {
        str[offset++] = '-';
    }

    for (size_t ind = 0; ind < pos; ind++) {
        str[offset + ind] = temp[pos - 1 - ind];
    }
    str[offset + pos] = '\0';

    return 1;
#endif
}


int
ibz_set_from_str(ibz_t *i, const char *str, int base)
{
    if (!str || (base != 10 && base != 16))
        return 0;

    ibz_init(i);

    int is_negative = 0;
    int pos = 0;

    if (str[0] == '-') {
        is_negative = 1;
        pos = 1;
    } else if (str[0] == '+') {
        pos = 1;
    }

    if (str[pos] == '\0')
        return 0;

    ibz_t magnitude, limit;
    ibz_init(&magnitude);
    if (is_negative) {
        ibz_init(&limit);
        limit[IBZ_LIMBS - 1] = UINT64_C(1) << 63;
    } else {
        for (int limb = 0; limb < IBZ_LIMBS; ++limb)
            limit[limb] = UINT64_MAX;
        limit[IBZ_LIMBS - 1] >>= 1;
    }

    while (str[pos] != '\0') {
        char c = str[pos];
        int digit;

        if (c >= '0' && c <= '9') {
            digit = c - '0';
        } else if (base == 16 && c >= 'a' && c <= 'f') {
            digit = c - 'a' + 10;
        } else if (base == 16 && c >= 'A' && c <= 'F') {
            digit = c - 'A' + 10;
        } else {
            return 0; // Invalid character
        }

        if (digit >= base) {
            return 0; // Invalid digit for base
        }

        /* magnitude = magnitude * base + digit, with an explicit unsigned
         * overflow check so parsing never wraps into a different value. */
        uint64_t carry = (uint64_t)digit;
        for (int limb = 0; limb < IBZ_LIMBS; ++limb) {
            uint64_t hi, lo;
            mul64_128(magnitude[limb], (uint64_t)base, &hi, &lo);
            uint64_t sum = lo + carry;
            uint64_t add_carry = sum < lo;
            magnitude[limb] = sum;
            if (hi > UINT64_MAX - add_carry) {
                ibz_init(i);
                return 0;
            }
            carry = hi + add_carry;
        }
        if (carry != 0 || ibz_cmp_unsigned(&magnitude, &limit) > 0) {
            ibz_init(i);
            return 0;
        }

        pos++;
    }

    if (is_negative && !ibz_is_zero(&magnitude))
        ibz_neg_raw(i, &magnitude);
    else
        ibz_copy(i, &magnitude);

    return 1;
}

// Get as int32
int32_t
ibz_get(const ibz_t *i)
{
    return (int32_t)((*i)[0] & 0xFFFFFFFF);
}

int
ibz_rand_interval(ibz_t *rand, const ibz_t *a, const ibz_t *b)
{
    const int order = ibz_cmp(a, b);
    if (order > 0) {
        ibz_init(rand);
        return 0;
    }
    if (order == 0) {
        ibz_copy(rand, a);
        return 1;
    }

    ibz_t range;
    if (ibz_is_negative(a) && !ibz_is_negative(b)) {
        /* b-a equals b+|a| and cannot wrap at this signed input range.  Its
         * sign bit therefore detects exactly the intervals whose inclusive
         * width is not representable, without three full-width temporaries. */
        ibz_sub_unsigned(&range, b, a);
        if (ibz_is_negative(&range)) {
            ibz_init(rand);
            return 0;
        }
    } else {
        ibz_sub(&range, b, a);
    }

    int len_bits = ibz_bitsize(&range);
    int len_bytes = (len_bits + 7) / 8;
    int len_limbs = (len_bytes + 7) / 8;

    uint64_t mask = (len_bits % 64 == 0) ? UINT64_MAX : ((1ULL << (len_bits % 64)) - 1);

    // Rejection sampling
    for (int tries = 0; tries < 1000; tries++) {
#if defined(TARGET_ARM) && !defined(TARGET_BIG_ENDIAN)
        /* STM32F407 is little-endian.  Decode directly into the output so
         * Cortex-M4 does not reserve a second IBZ-sized entropy buffer. */
        unsigned char *bytes = (unsigned char *)(void *)*rand;
        ibz_init(rand);
#else
        unsigned char bytes[IBZ_LIMBS * sizeof(uint64_t)] = { 0 };
        ibz_init(rand);
#endif

        if (sqisign_randombytes(bytes, len_bytes) != 0) {
#if !defined(TARGET_ARM) || defined(TARGET_BIG_ENDIAN)
            sqisign_randombytes_wipe(bytes, sizeof(bytes));
#else
            ibz_init(rand);
#endif
            return 0;
        }

#if !defined(TARGET_ARM) || defined(TARGET_BIG_ENDIAN)
        /* Decode the entropy explicitly as a little-endian unsigned integer.
         * Casting rand to bytes made the partial most-significant limb and its
         * mask host-endian-dependent. */
        for (int byte = 0; byte < len_bytes; ++byte)
            (*rand)[byte / 8] |= (uint64_t)bytes[byte] << (8 * (byte % 8));
        sqisign_randombytes_wipe(bytes, sizeof(bytes));
#endif

        if (len_limbs > 0 && len_limbs <= IBZ_LIMBS) {
            (*rand)[len_limbs - 1] &= mask;
        }

        /* range = b - a, so accepting range itself is required for the
         * documented inclusive interval [a, b]. */
        if (ibz_cmp(rand, &range) <= 0) {
            ibz_add(rand, rand, a);
            return 1;
        }
    }

    return 0;
}

int
ibz_rand_interval_i(ibz_t *rand, int32_t a, int32_t b)
{
    ibz_t a_ibz, b_ibz;
    ibz_set(&a_ibz, a);
    ibz_set(&b_ibz, b);
    return ibz_rand_interval(rand, &a_ibz, &b_ibz);
}

int
ibz_rand_interval_minm_m(ibz_t *rand, int32_t m)
{
    if (m < 0) {
        ibz_init(rand);
        return 0;
    }

    const uint32_t range = 2U * (uint32_t)m;
    if (range == 0) {
        ibz_init(rand);
        return 1;
    }
    const int len_bits = 64 - clz64(range);
    const int len_bytes = (len_bits + 7) / 8;
    const uint32_t mask = len_bits == 32
                              ? UINT32_MAX
                              : (UINT32_C(1) << len_bits) - UINT32_C(1);

    for (int tries = 0; tries < 1000; ++tries) {
        unsigned char bytes[sizeof(uint32_t)] = { 0 };
        if (sqisign_randombytes(bytes, (size_t)len_bytes) != 0) {
            sqisign_randombytes_wipe(bytes, sizeof(bytes));
            ibz_init(rand);
            return 0;
        }
        uint32_t sample = 0;
        for (int byte = 0; byte < len_bytes; ++byte)
            sample |= (uint32_t)bytes[byte] << (8 * byte);
        sqisign_randombytes_wipe(bytes, sizeof(bytes));
        sample &= mask;
        if (sample <= range) {
            const int32_t signed_sample =
                (int32_t)((int64_t)sample - (int64_t)m);
            ibz_set(rand, signed_sample);
            return 1;
        }
    }

    ibz_init(rand);
    return 0;
}

int
ibz_rand_interval_bits(ibz_t *rand, uint32_t m)
{
    if (m >= IBZ_BITS - 1) {
        ibz_init(rand);
        return 0;
    }

    /* Sampling [0,2^(m+1)] and subtracting 2^m is exactly the generic
     * inclusive interval algorithm, but needs no min/max/range objects. */
    const uint32_t len_bits = m + 2;
    const uint32_t len_bytes = (len_bits + 7) / 8;
    const uint32_t len_limbs = (len_bytes + 7) / 8;
    const uint32_t top_bits = len_bits % 64;
    const uint64_t top_mask = top_bits == 0
                                  ? UINT64_MAX
                                  : (UINT64_C(1) << top_bits) - UINT64_C(1);
    const uint32_t boundary_bit = m + 1;
    const uint32_t boundary_limb = boundary_bit / 64;
    const uint64_t boundary_mask =
        UINT64_C(1) << (boundary_bit % 64);

    for (int tries = 0; tries < 1000; ++tries) {
#if defined(TARGET_ARM) && !defined(TARGET_BIG_ENDIAN)
        unsigned char *bytes = (unsigned char *)(void *)*rand;
        ibz_init(rand);
#else
        unsigned char bytes[IBZ_LIMBS * sizeof(uint64_t)] = { 0 };
        ibz_init(rand);
#endif
        if (sqisign_randombytes(bytes, len_bytes) != 0) {
#if !defined(TARGET_ARM) || defined(TARGET_BIG_ENDIAN)
            sqisign_randombytes_wipe(bytes, sizeof(bytes));
#else
            ibz_init(rand);
#endif
            return 0;
        }
#if !defined(TARGET_ARM) || defined(TARGET_BIG_ENDIAN)
        for (uint32_t byte = 0; byte < len_bytes; ++byte)
            (*rand)[byte / 8] |=
                (uint64_t)bytes[byte] << (8 * (byte % 8));
        sqisign_randombytes_wipe(bytes, sizeof(bytes));
#endif
        (*rand)[len_limbs - 1] &= top_mask;

        int accept = ((*rand)[boundary_limb] & boundary_mask) == 0;
        if (!accept && (*rand)[boundary_limb] == boundary_mask) {
            accept = 1;
            for (uint32_t limb = 0; limb < boundary_limb; ++limb)
                accept &= (*rand)[limb] == 0;
        }
        if (!accept)
            continue;

        const uint32_t subtract_limb = m / 64;
        const uint64_t subtract_mask = UINT64_C(1) << (m % 64);
        const uint64_t original = (*rand)[subtract_limb];
        (*rand)[subtract_limb] = original - subtract_mask;
        uint64_t borrow = original < subtract_mask;
        for (uint32_t limb = subtract_limb + 1;
             limb < IBZ_LIMBS;
             ++limb) {
            const uint64_t word = (*rand)[limb];
            (*rand)[limb] = word - borrow;
            borrow = word < borrow;
        }
        return 1;
    }

    ibz_init(rand);
    return 0;
}

// Compare with int32
int
ibz_cmp_int32(const ibz_t *x, int32_t y)
{
    const uint64_t extension = y < 0 ? UINT64_MAX : 0;
    const uint64_t low = (uint64_t)(int64_t)y;
    const int x_negative = ibz_is_negative(x);
    const int y_negative = y < 0;

    if (x_negative != y_negative)
        return x_negative ? -1 : 1;

    for (int i = IBZ_LIMBS - 1; i >= 0; --i) {
        const uint64_t yi = i == 0 ? low : extension;
        if ((*x)[i] > yi)
            return 1;
        if ((*x)[i] < yi)
            return -1;
    }
    return 0;
}

// Get bit size (absolute value)
int
ibz_bitsize(const ibz_t *a)
{
    if (!ibz_is_negative(a))
        return ibz_bitsize_unsigned(a);

    /* Stream the two's-complement negation and retain only its highest
     * nonzero limb.  Materializing |a| served no other purpose here. */
    uint64_t carry = 1;
    uint64_t highest = 0;
    int highest_index = -1;
    for (int i = 0; i < IBZ_LIMBS; ++i) {
        const uint64_t inverted = ~(*a)[i];
        const uint64_t limb = inverted + carry;
        carry = limb < inverted;
        if (limb != 0) {
            highest = limb;
            highest_index = i;
        }
    }
    return highest_index < 0
               ? 0
               : highest_index * 64 + (64 - clz64(highest));
}

// Get size in base
size_t
ibz_size_in_base(const ibz_t *a, int base)
{
    if (base != 10 && base != 16)
        return 0;
    if (ibz_is_zero(a)) {
        return 1;
    }

    ibz_t temp;
    ibz_abs_unsigned(&temp, a);

    size_t count = 0;

    while (!ibz_is_zero(&temp)) {
        ibz_t q;
        (void)ibz_div_unsigned_small(&q, &temp, (uint32_t)base);
        ibz_copy(&temp, &q);
        count++;
    }

    return count > 0 ? count : 1;
}

size_t
ibz_digits_required(const ibz_t *a)
{
    const int bits = ibz_bitsize(a);
    return bits == 0 ? 1 : ((size_t)bits + 63) / 64;
}

int
ibz_to_u64_digits_checked(ibz_digit_t *target, size_t target_len, const ibz_t *a)
{
    if (target == NULL || target_len == 0 ||
        target_len > SIZE_MAX / sizeof(*target))
        return 0;

    ibz_t magnitude;
    ibz_abs_unsigned(&magnitude, a);
    const size_t required = ibz_digits_required(a);
    memset(target, 0, target_len * sizeof(*target));
    if (required > target_len)
        return 0;

    memcpy(target, magnitude, required * sizeof(*target));
    return 1;
}

// Convert the unsigned magnitude to a digit array.  The legacy entry point
// requires at least ibz_digits_required(a) output elements; new code should
// use ibz_to_u64_digits_checked when the destination is externally sized.
void
ibz_to_u64_digits(ibz_digit_t *target, const ibz_t *a)
{
    const size_t required = ibz_digits_required(a);
    (void)ibz_to_u64_digits_checked(target, required, a);
}

int
ibz_to_u32_digits_checked(uint32_t *target, size_t target_len, const ibz_t *a)
{
    if (target == NULL || target_len == 0 ||
        target_len > SIZE_MAX / sizeof(*target))
        return 0;

    ibz_t magnitude;
    ibz_abs_unsigned(&magnitude, a);

    size_t required = 1;
    for (int i = IBZ_LIMBS - 1; i >= 0; --i) {
        if (magnitude[i] != 0) {
            required = 2u * (size_t)i + ((magnitude[i] >> 32) != 0 ? 2u : 1u);
            break;
        }
    }

    memset(target, 0, target_len * sizeof(*target));
    if (required > target_len)
        return 0;

    for (size_t i = 0; i < required; ++i) {
        const uint64_t limb = magnitude[i / 2u];
        target[i] = (uint32_t)(limb >> (32u * (i & 1u)));
    }
    return 1;
}

// Copy from a little-endian 64-bit digit array.
void
ibz_copy_u64_digits(ibz_t *a, const ibz_digit_t *digits, size_t len)
{
    ibz_init(a);
    for (size_t i = 0; i < len && i < IBZ_LIMBS; i++) {
        (*a)[i] = digits[i];
    }
}

// Copy from a little-endian 32-bit field/MP digit array.
void
ibz_copy_u32_digits(ibz_t *a, const uint32_t *digits, size_t len)
{
    ibz_init(a);
    const size_t max_digits = 2u * (size_t)IBZ_LIMBS;
    if (len > max_digits)
        len = max_digits;
    for (size_t i = 0; i < len; ++i)
        (*a)[i / 2u] |= (uint64_t)digits[i] << (32u * (i & 1u));
}

// Get 2-adic valuation (trailing zeros)
int
ibz_two_adic(const ibz_t *pow)
{
    for (int i = 0; i < IBZ_LIMBS; i++) {
        if ((*pow)[i] != 0) {
            return i * 64 + ctz64((*pow)[i]);
        }
    }
    return IBZ_BITS;
}

// void
// ibz_div(ibz_t *quotient, ibz_t *remainder, const ibz_t *a, const ibz_t *b)
// {
//     if (ibz_is_zero(b)) {
//         ibz_init(quotient);
//         ibz_init(remainder);
//         return;
//     }

//     ibz_init(quotient);
//     ibz_init(remainder);

//     int a_neg = ibz_is_negative(a);
//     int b_neg = ibz_is_negative(b);
//     int quot_neg = (a_neg != b_neg);

//     ibz_t dividend, divisor;
//     ibz_abs(&dividend, a);
//     ibz_abs(&divisor, b);

//     if (ibz_cmp(&dividend, &divisor) < 0) {
//         ibz_copy(remainder, a);
//         return;
//     }

//     // OPTIMIZATION: Use word-by-word division for better performance
//     ibz_t q, r;
//     ibz_init(&q);
//     ibz_copy(&r, &dividend);

//     int divisor_bits = ibz_bitsize(&divisor);
//     int dividend_bits = ibz_bitsize(&dividend);

//     // Shift divisor to align with dividend
//     int shift = dividend_bits - divisor_bits;
//     ibz_t shifted_divisor;
//     ibz_mul_2exp(&shifted_divisor, &divisor, shift);

//     // Division loop
//     for (int i = shift; i >= 0; i--) {
//         if (ibz_cmp(&r, &shifted_divisor) >= 0) {
//             ibz_sub(&r, &r, &shifted_divisor);
//             q[i / 64] |= (1ULL << (i % 64));
//         }
//         if (i > 0) {
//             ibz_div_2exp(&shifted_divisor, &shifted_divisor, 1);
//         }
//     }

//     if (quot_neg && !ibz_is_zero(&q)) {
//         ibz_neg(&q, &q);
//     }

//     if (a_neg && !ibz_is_zero(&r)) {
//         ibz_neg(&r, &r);
//     }

//     ibz_copy(quotient, &q);
//     ibz_copy(remainder, &r);
// }

void
ibz_div(ibz_t *quotient, ibz_t *remainder, const ibz_t *a, const ibz_t *b)
{
    ibz_t q, r, dividend, divisor;

    /* q/r are local until the final publication, so inputs remain intact for
     * every supported output/input alias pattern without separate copies. */
    if (ibz_is_zero(b))
        ibz_division_by_zero_abort("division");

    int a_neg   = ibz_is_negative(a);
    int b_neg   = ibz_is_negative(b);
    int quot_neg = (a_neg != b_neg);

    ibz_abs_unsigned(&dividend, a);
    ibz_abs_unsigned(&divisor, b);
    ibz_div_unsigned(&q, &r, &dividend, &divisor);

    if (quot_neg && !ibz_is_zero(&q))
        ibz_neg_raw(&q, &q);
    if (a_neg && !ibz_is_zero(&r))
        ibz_neg_raw(&r, &r);

#if defined(SQISIGN_INTBIG_OVERFLOW_CHECK)
    /* INTBIG_MIN / -1 is the sole division result outside the signed range. */
    if (!quot_neg && ibz_is_negative(&q))
        ibz_overflow_abort("division");
#endif

    // emit the result only at the end (alias-safe)
    if (quotient)  ibz_copy(quotient,  &q);
    if (remainder) ibz_copy(remainder, &r);
}


// Floor division
void
ibz_div_floor(ibz_t *q, ibz_t *r, const ibz_t *n, const ibz_t *d)
{
    ibz_t q_tmp, r_tmp;
    const int opposite_signs = ibz_is_negative(n) != ibz_is_negative(d);
    ibz_div(&q_tmp, &r_tmp, n, d);

    /* Truncating division differs from floor exactly when n/d is negative
     * and non-integral.  The floor remainder then has the divisor's sign. */
    if (!ibz_is_zero(&r_tmp) && opposite_signs) {
        ibz_sub(&q_tmp, &q_tmp, (const ibz_t *)&ibz_const_one);
        ibz_add(&r_tmp, &r_tmp, d);
    }

    if (q)
        ibz_copy(q, &q_tmp);
    if (r)
        ibz_copy(r, &r_tmp);
}

// Modulo
void
ibz_mod(ibz_t *r, const ibz_t *a, const ibz_t *b)
{
    ibz_t modulus;
    const int a_negative = ibz_is_negative(a);
    const int b_negative = ibz_is_negative(b);

    /* Only the floor remainder is requested.  Reducing magnitudes directly
     * avoids constructing a quotient and the nested public-division frames.
     * Preserve b first, then use the public output as the remainder scratch;
     * this keeps every r/a/b alias pattern valid with only one local ibz_t. */
    ibz_abs_unsigned(&modulus, b);
    if (a_negative)
        ibz_neg_raw(r, a);
    else if (r != a)
        ibz_copy(r, a);
    ibz_div_unsigned(NULL, r, r, &modulus);

    if (a_negative && !ibz_is_zero(r))
        ibz_sub_unsigned(r, &modulus, r);
    if (b_negative && !ibz_is_zero(r))
        ibz_sub_unsigned(r, r, &modulus);
}

/* Reduce a signed value modulo a positive unsigned magnitude.  The modulus
 * may be 2^(IBZ_BITS-1), which is not a positive signed ibz_t. */
static void
ibz_mod_positive_magnitude(ibz_t *r, const ibz_t *a, const ibz_t *modulus)
{
    ibz_t reduced;
    ibz_abs_unsigned(&reduced, a);
    ibz_div_unsigned(NULL, &reduced, &reduced, modulus);
    if (ibz_is_negative(a) && !ibz_is_zero(&reduced))
        ibz_sub_unsigned(&reduced, modulus, &reduced);
    ibz_copy(r, &reduced);
}

#if !defined(HAVE_UINT128)
/* a,b are in [0, modulus); modulus is an unsigned positive magnitude. */
static void
ibz_add_mod_positive(ibz_t *sum,
                     const ibz_t *a,
                     const ibz_t *b,
                     const ibz_t *modulus)
{
    uint64_t carry = 0;
    for (int i = 0; i < IBZ_LIMBS; ++i) {
        const uint64_t ai = (*a)[i];
        const uint64_t bi = (*b)[i];
        const uint64_t partial = ai + carry;
        const uint64_t first_carry = partial < ai;
        const uint64_t value = partial + bi;
        carry = first_carry | (value < partial);
        (*sum)[i] = value;
    }

    /* a,b < modulus <= 2^(IBZ_BITS-1), so the full unsigned sum cannot
     * carry out of the fixed representation and one subtraction suffices. */
    if (ibz_cmp_unsigned(sum, modulus) >= 0)
        ibz_sub_unsigned(sum, sum, modulus);
}

/* CIOS Montgomery multiplication with radix B=2^32 and R=B^words.
 *
 * Preconditions: modulus is odd and greater than one, and a,b are unsigned
 * representatives in [0,modulus).  The accumulator array is exactly one
 * ibz_t (W bytes), which is smaller than the portable modular-product frame.
 * Results are published only after the final modulus read, so out may alias
 * a, b, or modulus.
 */
static INTBIG_NOINLINE void
ibz_mont_mul_positive(ibz_t *out,
                      const ibz_t *a,
                      const ibz_t *b,
                      const ibz_t *modulus)
{
    uint32_t accumulator[2 * IBZ_LIMBS];
    const size_t words = ibz_used_words32(modulus);
    const uint32_t n0 = ibz_word32(modulus, 0);
    uint32_t inverse = 1;
    uint32_t top = 0;

    /* For odd n0, five Newton steps lift the inverse from 1 to 32 bits. */
    for (unsigned int step = 0; step < 5; ++step)
        inverse *= 2u - n0 * inverse;
    const uint32_t n0_prime = 0u - inverse;

    for (size_t j = 0; j < words; ++j)
        accumulator[j] = 0;

    for (size_t i = 0; i < words; ++i) {
        const uint32_t bi = ibz_word32(b, i);
        uint64_t carry = 0;

        /* Add a*b_i.  top+carry can equal B, so preserve all 33 bits. */
        for (size_t j = 0; j < words; ++j) {
            const uint64_t value =
                (uint64_t)ibz_word32(a, j) * bi + accumulator[j] + carry;
            accumulator[j] = (uint32_t)value;
            carry = value >> 32;
        }
        uint64_t high = (uint64_t)top + carry;

        /* q cancels the low word because n0*n0_prime == -1 modulo B. */
        const uint32_t q = accumulator[0] * n0_prime;
        carry = 0;
        for (size_t j = 0; j < words; ++j) {
            const uint64_t value =
                (uint64_t)q * ibz_word32(modulus, j) +
                accumulator[j] + carry;
#if defined(SQISIGN_INTBIG_OVERFLOW_CHECK)
            if (j == 0 && (uint32_t)value != 0)
                ibz_overflow_abort("Montgomery cancellation");
#endif
            if (j != 0)
                accumulator[j - 1] = (uint32_t)value;
            carry = value >> 32;
        }

        /* This is the integrated exact division by B.  If A<2m before this
         * iteration, then A+a*b_i+q*m<2Bm, hence the new top is at most one. */
        high += carry; /* high <= 2*B-1, so uint64_t is sufficient. */
        accumulator[words - 1] = (uint32_t)high;
        top = (uint32_t)(high >> 32);
#if defined(SQISIGN_INTBIG_OVERFLOW_CHECK)
        if (top > 1)
            ibz_overflow_abort("Montgomery top carry");
#endif
    }

    /* The candidate is below 2*modulus, so at most one subtraction is needed.
     * A nonzero top necessarily makes the candidate at least R>modulus. */
    int subtract = top != 0;
    if (!subtract) {
        /* Equality also requires subtraction; only a first differing word
         * below the modulus clears this default. */
        subtract = 1;
        for (size_t j = words; j-- != 0;) {
            const uint32_t tj = accumulator[j];
            const uint32_t nj = ibz_word32(modulus, j);
            if (tj != nj) {
                subtract = tj > nj;
                break;
            }
        }
    }

    if (subtract) {
        uint32_t borrow = 0;
        for (size_t j = 0; j < words; ++j) {
            const uint32_t tj = accumulator[j];
            const uint32_t nj = ibz_word32(modulus, j);
            const uint32_t difference = tj - nj;
            const uint32_t next_borrow =
                (tj < nj) | ((tj == nj) & borrow);
            accumulator[j] = difference - borrow;
            borrow = next_borrow;
        }
        top -= borrow;
    }

#if defined(SQISIGN_INTBIG_OVERFLOW_CHECK)
    if (top != 0)
        ibz_overflow_abort("Montgomery final subtraction");
#endif

    for (size_t limb = 0; limb < IBZ_LIMBS; ++limb) {
        const size_t low_word = 2 * limb;
        uint64_t value = 0;
        if (low_word < words)
            value = accumulator[low_word];
        if (low_word + 1 < words)
            value |= (uint64_t)accumulator[low_word + 1] << 32;
        (*out)[limb] = value;
    }
}
#endif

/* Overflow-free modular multiplication for fixed integers.  On 128-bit
 * capable targets the complete 2N-limb product is reduced directly, so no
 * high product limbs are lost at the fixed-width boundary.  The portable
 * fallback retains add-and-double for products that do not fit. */
static void
ibz_mul_mod_positive(ibz_t *product,
                     const ibz_t *a,
                     const ibz_t *b,
                     const ibz_t *modulus)
{
#if defined(HAVE_UINT128)
    const size_t a_size = ibz_used_limbs(*a, IBZ_LIMBS);
    const size_t b_size = ibz_used_limbs(*b, IBZ_LIMBS);
    if (a_size == 0 || b_size == 0) {
        ibz_init(product);
        return;
    }

    uint64_t wide_product[2 * IBZ_LIMBS];
    uint64_t reduced[IBZ_LIMBS];
    ibz_mul_unsigned_wide(wide_product, *a, a_size, *b, b_size);
    ibz_divrem_unsigned_wide(NULL,
                             reduced,
                             wide_product,
                             a_size + b_size,
                             *modulus);
    memcpy(*product, reduced, sizeof(ibz_t));
#else
    const int a_bits = ibz_bitsize_unsigned(a);
    const int b_bits = ibz_bitsize_unsigned(b);
    if (a_bits == 0 || b_bits == 0) {
        ibz_init(product);
        return;
    }

    if (a_bits <= (IBZ_BITS - 1) - b_bits) {
        ibz_t unreduced;
        ibz_mul(&unreduced, a, b);
        ibz_div_unsigned(NULL, product, &unreduced, modulus);
        return;
    }

    const ibz_t *multiplicand = a;
    const ibz_t *multiplier = b;
    int multiplier_bits = b_bits;
    if (a_bits < b_bits) {
        multiplicand = b;
        multiplier = a;
        multiplier_bits = a_bits;
    }

    ibz_t result, addend;
    ibz_init(&result);
    ibz_copy(&addend, multiplicand);

    for (int bit = 0; bit < multiplier_bits; ++bit) {
        if ((((*multiplier)[bit / 64] >> (bit % 64)) & UINT64_C(1)) != 0)
            ibz_add_mod_positive(&result, &result, &addend, modulus);
        if (bit + 1 < multiplier_bits)
            ibz_add_mod_positive(&addend, &addend, &addend, modulus);
    }
    ibz_copy(product, &result);
#endif
}

// Power
void
ibz_pow(ibz_t *pow, const ibz_t *x, uint32_t e)
{
    if (e == 0) {
        ibz_set(pow, 1);
        return;
    }

    ibz_t result;
    ibz_copy(&result, x);
    uint32_t bit = UINT32_C(1) << 31;
    while ((e & bit) == 0)
        bit >>= 1;
    for (bit >>= 1; bit != 0; bit >>= 1) {
        ibz_mul(&result, &result, &result);
        if ((e & bit) != 0)
            ibz_mul(&result, &result, x);
    }

    ibz_copy(pow, &result);
}

// Modulo for unsigned long
unsigned long
ibz_mod_ui(const ibz_t *n, unsigned long d)
{
    if (d == 0)
        ibz_division_by_zero_abort("modulo");

    /* This is the path used by primality trial division.  Reducing two
     * 32-bit halves per active limb avoids constructing an unused quotient. */
    if (d <= UINT32_MAX) {
        if (!ibz_is_negative(n))
            return (unsigned long)ibz_mod_unsigned_small(n, (uint32_t)d);

        ibz_t magnitude;
        ibz_abs_unsigned(&magnitude, n);
        const unsigned long remainder =
            (unsigned long)ibz_mod_unsigned_small(&magnitude, (uint32_t)d);
        return remainder == 0 ? 0 : d - remainder;
    }

    ibz_t magnitude;
    ibz_abs_unsigned(&magnitude, n);
    ibz_t divisor, remainder;
    ibz_set_u64(&divisor, d);
    ibz_div_unsigned(NULL, &remainder, &magnitude, &divisor);
    if (ibz_is_negative(n) && !ibz_is_zero(&remainder))
        ibz_sub_unsigned(&remainder, &divisor, &remainder);

    return (unsigned long)remainder[0];
}

// Probabilistic primality test (Miller-Rabin)
int
ibz_probab_prime(const ibz_t *n, int reps)
{
    if (reps <= 0)
        return 0;
    if (ibz_cmp_int32(n, 2) == 0)
        return 1;
    if (ibz_cmp_int32(n, 3) == 0)
        return 1;
    if (ibz_is_even(n))
        return 0;
    if (ibz_cmp_int32(n, 1) <= 0)
        return 0;

    /* Match the fixed-precision implementation's cheap composite filter.
     * Most candidates produced by quat_represent_integer have a small prime
     * factor; rejecting them here avoids entering a full Miller--Rabin
     * modular exponentiation. */
    static const unsigned long small_primes[50] = {
          3,   5,   7,  11,  13,  17,  19,  23,  29,  31,
         37,  41,  43,  47,  53,  59,  61,  67,  71,  73,
         79,  83,  89,  97, 101, 103, 107, 109, 113, 127,
        131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
        181, 191, 193, 197, 199, 211, 223, 227, 229, 233
    };
    for (size_t i = 0; i < sizeof(small_primes) / sizeof(small_primes[0]); ++i) {
        const unsigned long prime = small_primes[i];
        if (ibz_cmp_int32(n, (int32_t)prime) == 0)
            return 1;
        if (ibz_mod_ui(n, prime) == 0)
            return 0;
    }

    ibz_t n_minus_1, d, a, x;
    ibz_sub(&n_minus_1, n, (const ibz_t *)&ibz_const_one);
    ibz_copy(&d, &n_minus_1);

    const int s = ibz_two_adic(&d);
    ibz_div_2exp(&d, &d, (uint32_t)s);

    for (int i = 0; i < reps; i++) {
        /* x is dead until modular exponentiation, so use it for n-2 instead
         * of reserving two more full-width objects in this hot frame. */
        ibz_sub(&x, n, (const ibz_t *)&ibz_const_two);
        if (!ibz_rand_interval(&a, (const ibz_t *)&ibz_const_two, &x))
            return 0;

        // x = a^d mod n
        ibz_pow_mod(&x, &a, &d, n);

        if (ibz_is_one(&x) || ibz_cmp(&x, &n_minus_1) == 0) {
            continue;
        }

        int composite = 1;
        for (int round = 1; round < s; ++round) {
            ibz_mul_mod_positive(&x, &x, &x, n);

            if (ibz_is_one(&x)) {
                return 0; // Composite
            }
            if (ibz_cmp(&x, &n_minus_1) == 0) {
                composite = 0;
                break;
            }
        }

        if (composite) {
            return 0;
        }
    }

    return 1; // Probably prime
}

// Modular exponentiation
void
ibz_pow_mod(ibz_t *pow, const ibz_t *x, const ibz_t *e, const ibz_t *m)
{
    ibz_t result, base, modulus;
    ibz_abs_unsigned(&modulus, m);
    if (ibz_is_zero(&modulus))
        ibz_division_by_zero_abort("modular exponentiation");
    if (ibz_is_negative(e))
        ibz_invalid_argument_abort("negative modular exponent");

    const int exponent_bits = ibz_bitsize_unsigned(e);
    if (ibz_is_one(&modulus)) {
        /* This also defines 0^0 mod 1 consistently with the old reduction. */
        ibz_init(&result);
        ibz_copy(pow, &result);
        return;
    }

    if (!ibz_is_negative(x) && ibz_cmp_unsigned(x, &modulus) < 0)
        ibz_copy(&base, x);
    else
        ibz_mod_positive_magnitude(&base, x, &modulus);

    if (exponent_bits == 0) {
        ibz_set(&result, 1);
    }
#if !defined(HAVE_UINT128)
    else if (ibz_is_odd(&modulus)) {
        const size_t words = ibz_used_words32(&modulus);

        /* Reuse result to form R^2 mod modulus without another full-width
         * object: R=2^(32*words), so R^2 requires 64*words doublings. */
        ibz_set(&result, 1);
        for (size_t doubling = 0; doubling < 64u * words; ++doubling)
            ibz_add_mod_positive(&result, &result, &result, &modulus);

        /* Mont(base,R^2)=base*R.  Starting from the mandatory top exponent
         * bit avoids materializing Montgomery one in another ibz_t. */
        ibz_mont_mul_positive(&base, &base, &result, &modulus);
        ibz_copy(&result, &base);
        for (int bit = exponent_bits - 2; bit >= 0; --bit) {
            ibz_mont_mul_positive(&result, &result, &result, &modulus);
            if (((*e)[bit / 64] & (UINT64_C(1) << (bit % 64))) != 0)
                ibz_mont_mul_positive(&result, &result, &base, &modulus);
        }

        /* Mont(result,1) removes the remaining factor R. */
        ibz_mont_mul_positive(
            &result, &result, (const ibz_t *)&ibz_const_one, &modulus);
    }
#endif
    else {
        /* Native-uint128 targets and even moduli keep the established path. */
        /* Left-to-right exponentiation starts from the mandatory top bit.
         * It removes the exponent copy and one modular multiplication. */
        ibz_copy(&result, &base);
        for (int bit = exponent_bits - 2; bit >= 0; --bit) {
            ibz_mul_mod_positive(&result, &result, &result, &modulus);
            if (((*e)[bit / 64] & (UINT64_C(1) << (bit % 64))) != 0)
                ibz_mul_mod_positive(&result, &result, &base, &modulus);
        }
    }

    ibz_copy(pow, &result);
}

// GCD using binary algorithm
void
ibz_gcd(ibz_t *gcd, const ibz_t *a, const ibz_t *b)
{
    ibz_t u, v;
    ibz_abs_unsigned(&u, a);
    ibz_abs_unsigned(&v, b);

    if (ibz_is_zero(&u)) {
        ibz_copy(gcd, &v);
        return;
    }
    if (ibz_is_zero(&v)) {
        ibz_copy(gcd, &u);
        return;
    }
    if (ibz_is_one(&u) || ibz_is_one(&v)) {
        ibz_set(gcd, 1);
        return;
    }

    int shift = MIN(ibz_two_adic(&u), ibz_two_adic(&v));
    ibz_shift_right_unsigned(&u, &u, (uint32_t)shift);
    ibz_shift_right_unsigned(&v, &v, (uint32_t)shift);

    ibz_t *u_work = &u;
    ibz_t *v_work = &v;
    while (!ibz_is_zero(u_work)) {
        ibz_shift_right_unsigned(u_work,
                                 u_work,
                                 (uint32_t)ibz_two_adic(u_work));
        ibz_shift_right_unsigned(v_work,
                                 v_work,
                                 (uint32_t)ibz_two_adic(v_work));

        if (ibz_cmp_unsigned(u_work, v_work) > 0) {
            ibz_t *swap = u_work;
            u_work = v_work;
            v_work = swap;
        }

        ibz_sub_unsigned(v_work, v_work, u_work);
    }

    // Restore common factors of 2 in the unsigned magnitude domain.
    ibz_shift_left_unsigned(v_work, v_work, (uint32_t)shift);

    ibz_copy(gcd, v_work);
}

void
ibz_gcdext(ibz_t *gcd, ibz_t *x, ibz_t *y,
           const ibz_t *a, const ibz_t *b)
{
    // special case: a == 0
    if (ibz_is_zero(a)) {
        ibz_abs_unsigned(gcd, b);
        ibz_set(x, 0);
        if (ibz_is_negative(b) && !ibz_is_min_value(b)) {
            ibz_set(y, -1);
        } else {
            ibz_set(y, 1);
        }
        return;
    }

    // special case: b == 0
    if (ibz_is_zero(b)) {
        ibz_abs_unsigned(gcd, a);
        if (ibz_is_negative(a) && !ibz_is_min_value(a)) {
            ibz_set(x, -1);
        } else {
            ibz_set(x, 1);
        }
        ibz_set(y, 0);
        return;
    }

    ibz_t aa, bb;
    ibz_t x0, x1, y0, y1;
    ibz_t q, r;
    ibz_t tmp1, tmp2;

    // copy a, b locally since they may be reused later
    ibz_copy(&aa, a);
    ibz_copy(&bb, b);

    // (x0, y0) = (1, 0), (x1, y1) = (0, 1)
    ibz_set(&x0, 1);
    ibz_set(&y0, 0);
    ibz_set(&x1, 0);
    ibz_set(&y1, 1);

    // standard extended Euclidean algorithm
    while (!ibz_is_zero(&bb)) {
        // aa = q * bb + r  (ibz_div: q, r with signs matching aa, bb)
        ibz_div(&q, &r, &aa, &bb);

        // (aa, bb) <- (bb, r)
        ibz_copy(&aa, &bb);
        ibz_copy(&bb, &r);

        // (x0, x1) <- (x1, x0 - q*x1)
        ibz_mul(&tmp1, &q, &x1);    // tmp1 = q * x1
        ibz_sub(&tmp2, &x0, &tmp1); // tmp2 = x0 - q*x1
        ibz_copy(&x0, &x1);
        ibz_copy(&x1, &tmp2);

        // (y0, y1) <- (y1, y0 - q*y1)
        ibz_mul(&tmp1, &q, &y1);    // tmp1 = q * y1
        ibz_sub(&tmp2, &y0, &tmp1); // tmp2 = y0 - q*y1
        ibz_copy(&y0, &y1);
        ibz_copy(&y1, &tmp2);
    }

    // now aa == gcd(a, b) (signed)
    ibz_copy(gcd, &aa);
    ibz_copy(x, &x0);
    ibz_copy(y, &y0);

    // normalize gcd to positive: if gcd < 0 flip the sign of gcd, x, y
    if (ibz_is_negative(gcd) && !ibz_is_min_value(gcd)) {
        ibz_neg(gcd, gcd);
        ibz_neg(x, x);
        ibz_neg(y, y);
    }
}


// void
// ibz_gcdext(ibz_t *gcd, ibz_t *x, ibz_t *y, const ibz_t *a, const ibz_t *b)
// {
//     if (ibz_is_zero(a)) {
//         ibz_abs(gcd, b);
//         ibz_set(x, 0);
//         if (ibz_is_negative(b)) {
//             ibz_set(y, -1);
//         } else {
//             ibz_set(y, 1);
//         }
//         return;
//     }

//     if (ibz_is_zero(b)) {
//         ibz_abs(gcd, a);
//         if (ibz_is_negative(a)) {
//             ibz_set(x, -1);
//         } else {
//             ibz_set(x, 1);
//         }
//         ibz_set(y, 0);
//         return;
//     }

//     int a_was_neg = ibz_is_negative(a);
//     int b_was_neg = ibz_is_negative(b);

//     ibz_t u, v, A, B, C, D, temp;

//     ibz_abs(&u, a);
//     ibz_abs(&v, b);

//     ibz_set(&A, 1);
//     ibz_set(&B, 0);
//     ibz_set(&C, 0);
//     ibz_set(&D, 1);

//     int shift = 0;
//     while (ibz_is_even(&u) && ibz_is_even(&v)) {
//         ibz_div_2exp(&u, &u, 1);
//         ibz_div_2exp(&v, &v, 1);
//         shift++;
//     }

//     ibz_t orig_u, orig_v;
//     ibz_copy(&orig_u, &u);
//     ibz_copy(&orig_v, &v);

//     // Binary Extended GCD
//     while (!ibz_is_zero(&u)) {
//         while (ibz_is_even(&u)) {
//             ibz_div_2exp(&u, &u, 1);

//             if (ibz_is_even(&A) && ibz_is_even(&B)) {
//                 ibz_div_2exp(&A, &A, 1);
//                 ibz_div_2exp(&B, &B, 1);
//             } else {
//                 // A = (A + orig_v) / 2
//                 ibz_add(&temp, &A, &orig_v);
//                 ibz_div_2exp(&A, &temp, 1);
//                 // B = (B - orig_u) / 2
//                 ibz_sub(&temp, &B, &orig_u);
//                 ibz_div_2exp(&B, &temp, 1);
//             }
//         }

//         while (ibz_is_even(&v)) {
//             ibz_div_2exp(&v, &v, 1);

//             if (ibz_is_even(&C) && ibz_is_even(&D)) {
//                 ibz_div_2exp(&C, &C, 1);
//                 ibz_div_2exp(&D, &D, 1);
//             } else {
//                 // C = (C + orig_v) / 2
//                 ibz_add(&temp, &C, &orig_v);
//                 ibz_div_2exp(&C, &temp, 1);
//                 // D = (D - orig_u) / 2
//                 ibz_sub(&temp, &D, &orig_u);
//                 ibz_div_2exp(&D, &temp, 1);
//             }
//         }

//         if (ibz_cmp(&u, &v) >= 0) {
//             ibz_sub(&u, &u, &v);
//             ibz_sub(&A, &A, &C);
//             ibz_sub(&B, &B, &D);
//         } else {
//             // v = v - u, C = C - A, D = D - B
//             ibz_sub(&v, &v, &u);
//             ibz_sub(&C, &C, &A);
//             ibz_sub(&D, &D, &B);
//         }
//     }

//     // gcd = v * 2^shift
//     ibz_copy(gcd, &v);
//     ibz_mul_2exp(gcd, gcd, shift);

//     ibz_copy(x, &C);
//     ibz_copy(y, &D);

//     if (a_was_neg) {
//         ibz_neg(x, x);
//     }
//     if (b_was_neg) {
//         ibz_neg(y, y);
//     }
// }

// Modular inverse: a^(-1) mod m
int
ibz_invmod(ibz_t *inv, const ibz_t *a, const ibz_t *mod)
{
    if (ibz_is_zero(mod) || ibz_cmp_int32(mod, 0) <= 0) {
        return 0;
    }

    if (ibz_is_zero(a)) {
        return 0;
    }

    ibz_t gcd, x, y;
    ibz_t a_mod;

    ibz_mod(&a_mod, a, mod);

    if (ibz_is_zero(&a_mod)) {
        return 0;
    }

    ibz_gcdext(&gcd, &x, &y, &a_mod, mod);

    if (!ibz_is_one(&gcd)) {
        return 0;
    }

    ibz_mod(inv, &x, mod);

    return 1;
}

// Check divisibility
int
ibz_divides(const ibz_t *a, const ibz_t *b)
{
    if (ibz_is_zero(b))
        return 0;
    ibz_t r;
    ibz_mod(&r, a, b);
    return ibz_is_zero(&r);
}

// Integer square root
int
ibz_sqrt(ibz_t *sqrt, const ibz_t *a)
{
    if (ibz_is_negative(a)) {
        return 0;
    }

    if (ibz_is_zero(a)) {
        ibz_init(sqrt);
        return 1;
    }

    ibz_t x, temp;

    /* Start from a power-of-two upper bound on sqrt(a).  Starting from a
     * itself initially only halves the approximation, so a roughly 2050-bit
     * level-5 ideal index needs more than the 1000 iterations allowed below
     * and is falsely reported as a non-square. */
    const int bits = ibz_bitsize(a);
    ibz_set(&x, 1);
    ibz_mul_2exp(&x, &x, (bits + 1) / 2);

    // Newton-Raphson: x_new = (x + a/x) / 2
    for (int iter = 0; iter < 1000; iter++) {
        ibz_div(&temp, NULL, a, &x);
        ibz_add(&temp, &x, &temp);
        ibz_div_2exp(&temp, &temp, 1);

        if (ibz_cmp(&temp, &x) >= 0)
            break;

        ibz_copy(&x, &temp);
    }

    ibz_mul(&temp, &x, &x);
    if (ibz_cmp(&temp, a) == 0) {
        ibz_copy(sqrt, &x);
        return 1;
    }

    return 0;
}

// Square root floor
void
ibz_sqrt_floor(ibz_t *sqrt, const ibz_t *a)
{
    if (ibz_is_negative(a) || ibz_is_zero(a)) {
        ibz_init(sqrt);
        return;
    }

    ibz_t x, temp;

    int bits = ibz_bitsize(a);
    ibz_set(&x, 1);
    ibz_mul_2exp(&x, &x, (bits + 1) / 2);

    // Newton-Raphson with proper termination
    for (int iter = 0; iter < 1000; iter++) {
        ibz_div(&temp, NULL, a, &x);
        ibz_add(&temp, &x, &temp);
        ibz_div_2exp(&temp, &temp, 1);

        if (ibz_cmp(&temp, &x) >= 0) {
            break;
        }

        ibz_copy(&x, &temp);
    }

    ibz_copy(sqrt, &x);
}

// Legendre symbol
int
ibz_legendre(const ibz_t *a, const ibz_t *p)
{
    ibz_t a_mod, exp, result;

    if (ibz_cmp_int32(p, 2) <= 0 || ibz_is_even(p))
        return 0;

    ibz_mod(&a_mod, a, p);

    if (ibz_is_zero(&a_mod)) {
        return 0;
    }

    ibz_sub(&exp, p, (const ibz_t *)&ibz_const_one);
    ibz_div_2exp(&exp, &exp, 1);

    ibz_pow_mod(&result, &a_mod, &exp, p);

    if (ibz_is_one(&result)) {
        return 1;
    } else {
        return -1;
    }
}

// Modular square root
// int
// ibz_sqrt_mod_p(ibz_t *sqrt, const ibz_t *a, const ibz_t *p)
// {
//     ibz_t a_mod;
//     ibz_mod(&a_mod, a, p);

//     if (ibz_is_zero(&a_mod)) {
//         ibz_init(sqrt);
//         return 1;
//     }

//     if (ibz_legendre(&a_mod, p) != 1) {
//         return 0;
//     }

//     ibz_t p_mod_4, three, four;
//     ibz_set(&three, 3);
//     ibz_set(&four, 4);
//     ibz_mod(&p_mod_4, p, &four);

//     if (ibz_cmp(&p_mod_4, &three) == 0) {
//         ibz_t exp;
//         ibz_add(&exp, p, (const ibz_t *)&ibz_const_one);
//         ibz_div_2exp(&exp, &exp, 2);
//         ibz_pow_mod(sqrt, &a_mod, &exp, p);
//         return 1;
//     }

//     return 0;
// }

// Modular square root for odd prime p (Tonelli–Shanks)
// returns 1 if sqrt exists (and writes it), 0 otherwise.
static INTBIG_NOINLINE int
ibz_sqrt_mod_p_tonelli(ibz_t *root, const ibz_t *n, const ibz_t *p)
{
    ibz_t q, c, x, t;
    ibz_sub(&q, p, (const ibz_t *)&ibz_const_one);
    const int s = ibz_two_adic(&q);
    if (s <= 0)
        return 0;
    ibz_div_2exp(&q, &q, (uint32_t)s);

    ibz_set(&c, 2);
    while (ibz_legendre(&c, p) != -1)
        ibz_add(&c, &c, (const ibz_t *)&ibz_const_one);
    ibz_pow_mod(&c, &c, &q, p);

    ibz_pow_mod(&t, n, &q, p);
    ibz_add(&q, &q, (const ibz_t *)&ibz_const_one);
    ibz_div_2exp(&q, &q, 1);
    ibz_pow_mod(&x, n, &q, p);

    int m = s;
    while (!ibz_is_one(&t)) {
        ibz_copy(&q, &t);
        int i;
        for (i = 1; i < m; ++i) {
            ibz_mul_mod_positive(&q, &q, &q, p);
            if (ibz_is_one(&q))
                break;
        }
        if (i == m)
            return 0;

        ibz_copy(&q, &c);
        for (int j = 0; j < m - i - 1; ++j)
            ibz_mul_mod_positive(&q, &q, &q, p);
        ibz_mul_mod_positive(&x, &x, &q, p);
        ibz_mul_mod_positive(&q, &q, &q, p);
        ibz_mul_mod_positive(&t, &t, &q, p);
        ibz_copy(&c, &q);
        m = i;
    }

    ibz_copy(root, &x);
    return 1;
}

int
ibz_sqrt_mod_p(ibz_t *sqrt, const ibz_t *a, const ibz_t *p)
{
    /* Do not publish until every input read is complete: Cornacchia calls
     * this with sqrt == a, and public callers may also alias sqrt with p. */
    if (ibz_cmp_int32(p, 2) < 0 ||
        (ibz_cmp_int32(p, 2) > 0 && ibz_is_even(p))) {
        ibz_init(sqrt);
        return 0;
    }

    ibz_t n, result, scratch;
    ibz_mod(&n, a, p);
    if (ibz_cmp_int32(p, 2) == 0 || ibz_is_zero(&n)) {
        ibz_copy(sqrt, &n);
        return 1;
    }

    /* Keep only the three direct-formula objects in this frame; the general
     * Tonelli--Shanks scratch lives in a separate noinline helper. */
    const unsigned long modulus_mod_8 = ibz_mod_ui(p, 8);
    if ((modulus_mod_8 & 3UL) == 3UL) {
        ibz_div_2exp(&scratch, p, 2);
        ibz_add(&scratch, &scratch, (const ibz_t *)&ibz_const_one);
        ibz_pow_mod(&result, &n, &scratch, p);
    } else if (modulus_mod_8 == 5) {
        /* Atkin's p == 5 (mod 8) formula.  First distinguish the two residue
         * cases with d = n^((p-1)/4).  The final square check also rejects a
         * non-residue without a separate Legendre exponentiation. */
        ibz_div_2exp(&scratch, p, 2);
        ibz_pow_mod(&result, &n, &scratch, p);
        if (ibz_is_one(&result)) {
            ibz_div_2exp(&scratch, p, 3);
            ibz_add(&scratch,
                    &scratch,
                    (const ibz_t *)&ibz_const_one);
            ibz_pow_mod(&result, &n, &scratch, p);
        } else {
            ibz_div_2exp(&scratch, p, 3);
            ibz_mul_mod_positive(&result,
                                 &n,
                                 (const ibz_t *)&ibz_const_two,
                                 p);
            ibz_mul_mod_positive(&result,
                                 &result,
                                 (const ibz_t *)&ibz_const_two,
                                 p);
            ibz_pow_mod(&result, &result, &scratch, p);
            ibz_mul_mod_positive(&scratch,
                                 &n,
                                 (const ibz_t *)&ibz_const_two,
                                 p);
            ibz_mul_mod_positive(&result, &scratch, &result, p);
        }
    } else {
        if (ibz_legendre(&n, p) != 1 ||
            !ibz_sqrt_mod_p_tonelli(&result, &n, p)) {
            ibz_init(sqrt);
            return 0;
        }
    }

    ibz_mul_mod_positive(&scratch, &result, &result, p);
    if (ibz_cmp(&scratch, &n) != 0) {
        ibz_init(sqrt);
        return 0;
    }
    ibz_copy(sqrt, &result);
    return 1;
}
