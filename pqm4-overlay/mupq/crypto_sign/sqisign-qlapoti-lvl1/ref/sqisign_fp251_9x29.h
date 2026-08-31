#ifndef SQISIGN_FP251_9X29_H
#define SQISIGN_FP251_9X29_H

/*
 * Scalar ARM Cortex-M4 field backend for the SQIsign level-I prime
 *
 *     p = 5 * 2^248 - 1.
 *
 * The representation follows the unsaturated 9-limb layout studied in
 * "SQIsign on ARM": eight 29-bit limbs and one 19-bit limb.  Values held in
 * fp251_t are canonical Montgomery residues in [0,p), with Montgomery radix
 * R = 2^(9*29) = 2^261.  The implementation is freestanding C11 and requires
 * only 32x32->64 integer multiplication, which is available on Cortex-M4.
 *
 * This module is independently testable.  It is not, by itself, the complete
 * SQIsign protocol.
 */

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SQISIGN_FP251_LIMBS 9u
#define SQISIGN_FP251_RADIX_BITS 29u
#define SQISIGN_FP251_ENCODED_BYTES 32u
#define SQISIGN_FP251_MASK UINT32_C(0x1fffffff)

/* Canonical Montgomery-domain field element. */
typedef struct {
    uint32_t limb[SQISIGN_FP251_LIMBS];
} fp251_t;

/*
 * Non-canonical unsaturated value used for statically scheduled lazy
 * reduction.  The arithmetic layer, not the type, tracks the proven bound.
 * Safe inputs to the helpers below have lower limbs < 8*2^29 and top limb
 * < 32*(5*2^16); see the function contracts.
 */
typedef struct {
    uint32_t limb[SQISIGN_FP251_LIMBS];
} fp251_loose_t;

/* F_{p^2} = F_p[i]/(i^2+1). */
typedef struct {
    fp251_t re;
    fp251_t im;
} fp2_251_t;

/* Constants, all in the documented domain. */
extern const fp251_t fp251_zero;
extern const fp251_t fp251_one;       /* Montgomery representation of 1. */
extern const fp251_t fp251_modulus;   /* Standard integer representation. */
extern const fp251_t fp251_r2;        /* R^2 mod p, standard limbs. */

/* Base-field construction and conversion. */
void fp251_set_zero(fp251_t *r);
void fp251_set_one(fp251_t *r);
void fp251_copy(fp251_t *r, const fp251_t *a);
void fp251_from_u32(fp251_t *r, uint32_t x);

/*
 * Decode/encode canonical little-endian 32-byte integers.
 * fp251_decode() returns 1 on success and 0 if the input is >= p.
 * The in-memory fp251_t result is in Montgomery representation.
 */
int fp251_decode(fp251_t *r, const uint8_t in[SQISIGN_FP251_ENCODED_BYTES]);
void fp251_encode(uint8_t out[SQISIGN_FP251_ENCODED_BYTES], const fp251_t *a);


/*
 * Scalar form of the paper's static lazy-reduction schedule.
 *
 * - from_canonical: bound 1.
 * - add_bound4: each input limb is below 4 times its base; output is below
 *   8 times the base.
 * - sub_bound4: computes a + 4p - b component-wise under the same bound.
 * - incomplete_reduce: propagates lower-limb carries only.
 * - complete_reduce: two incomplete passes with a fixed Barrett fold of the
 *   top limb, then canonicalizes into [0,p).
 *
 * These functions have a fixed instruction schedule.  Call sites must place
 * reduction points using a static bound analysis; no runtime bound test is
 * performed.
 */
void fp251_loose_from_canonical(fp251_loose_t *r, const fp251_t *a);
void fp251_loose_add_bound4(fp251_loose_t *r, const fp251_loose_t *a,
                            const fp251_loose_t *b);
void fp251_loose_sub_bound4(fp251_loose_t *r, const fp251_loose_t *a,
                            const fp251_loose_t *b);
void fp251_loose_incomplete_reduce(fp251_loose_t *r,
                                   const fp251_loose_t *a);
void fp251_loose_complete_reduce(fp251_t *r, const fp251_loose_t *a);

/* Constant-schedule canonical field arithmetic. */
void fp251_add(fp251_t *r, const fp251_t *a, const fp251_t *b);
void fp251_sub(fp251_t *r, const fp251_t *a, const fp251_t *b);
void fp251_neg(fp251_t *r, const fp251_t *a);
void fp251_double(fp251_t *r, const fp251_t *a);
void fp251_mul(fp251_t *r, const fp251_t *a, const fp251_t *b);
void fp251_sqr(fp251_t *r, const fp251_t *a);

/* Returns 0 for input zero; otherwise computes a^(p-2). */
int fp251_inv(fp251_t *r, const fp251_t *a);

/*
 * Since p == 3 mod 4, sqrt(a)=a^((p+1)/4) when a is a square.
 * Returns 1 and a root on success, otherwise returns 0 and zeroes r.
 */
int fp251_sqrt(fp251_t *r, const fp251_t *a);

/* Returns 0 for zero, 1 for a nonzero square, and -1 for a nonsquare. */
int fp251_legendre(const fp251_t *a);

/* Comparison returns are 0 or 1. */
int fp251_iszero(const fp251_t *a);
int fp251_equal(const fp251_t *a, const fp251_t *b);

/* Explicitly wipe a field object. */
void fp251_secure_clear(fp251_t *a);

/* Extension-field arithmetic. */
void fp2_251_set_zero(fp2_251_t *r);
void fp2_251_set_one(fp2_251_t *r);
void fp2_251_copy(fp2_251_t *r, const fp2_251_t *a);
void fp2_251_add(fp2_251_t *r, const fp2_251_t *a, const fp2_251_t *b);
void fp2_251_sub(fp2_251_t *r, const fp2_251_t *a, const fp2_251_t *b);
void fp2_251_neg(fp2_251_t *r, const fp2_251_t *a);
void fp2_251_mul(fp2_251_t *r, const fp2_251_t *a, const fp2_251_t *b);
void fp2_251_sqr(fp2_251_t *r, const fp2_251_t *a);
int fp2_251_inv(fp2_251_t *r, const fp2_251_t *a);
int fp2_251_iszero(const fp2_251_t *a);
int fp2_251_equal(const fp2_251_t *a, const fp2_251_t *b);
void fp2_251_secure_clear(fp2_251_t *a);

/*
 * Two-way scalar batching for independent operations.  These routines do not
 * use NEON (which Cortex-M4 does not have); they expose independent dependency
 * chains so an M4 compiler/scheduler can interleave multiply-accumulates.
 */
void fp251_mul_2way(fp251_t out[2], const fp251_t a[2], const fp251_t b[2]);
void fp2_251_mul_2way(fp2_251_t out[2], const fp2_251_t a[2],
                      const fp2_251_t b[2]);

#ifdef __cplusplus
}
#endif

#endif /* SQISIGN_FP251_9X29_H */
