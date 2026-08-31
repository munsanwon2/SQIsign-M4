#ifndef SQISIGN_FP505_19X27_H
#define SQISIGN_FP505_19X27_H

/*
 * Scalar ARM Cortex-M4 field backend for the SQIsign level-5 prime
 *
 *     p = 27 * 2^500 - 1.
 *
 * The representation follows the unsaturated 19-limb layout studied in
 * "SQIsign on ARM": 18 27-bit limbs and one 19-bit limb.  Values held in
 * fp505_t are canonical Montgomery residues in [0,p), with Montgomery radix
 * R = 2^(19*27) = 2^513.  The implementation is freestanding C11 and requires
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

#define SQISIGN_FP505_LIMBS 19u
#define SQISIGN_FP505_RADIX_BITS 27u
#define SQISIGN_FP505_ENCODED_BYTES 64u
#define SQISIGN_FP505_MASK UINT32_C(0x07ffffff)

/* Canonical Montgomery-domain field element. */
typedef struct {
    uint32_t limb[SQISIGN_FP505_LIMBS];
} fp505_t;

/*
 * Non-canonical unsaturated value used for statically scheduled lazy
 * reduction.  The arithmetic layer, not the type, tracks the proven bound.
 * Safe inputs to the helpers below have lower limbs < 8*2^27 and top limb
 * < 32*(27*2^14); see the function contracts.
 */
typedef struct {
    uint32_t limb[SQISIGN_FP505_LIMBS];
} fp505_loose_t;

/* F_{p^2} = F_p[i]/(i^2+1). */
typedef struct {
    fp505_t re;
    fp505_t im;
} fp2_505_t;

/* Constants, all in the documented domain. */
extern const fp505_t fp505_zero;
extern const fp505_t fp505_one;       /* Montgomery representation of 1. */
extern const fp505_t fp505_modulus;   /* Standard integer representation. */
extern const fp505_t fp505_r2;        /* R^2 mod p, standard limbs. */

/* Base-field construction and conversion. */
void fp505_set_zero(fp505_t *r);
void fp505_set_one(fp505_t *r);
void fp505_copy(fp505_t *r, const fp505_t *a);
void fp505_from_u32(fp505_t *r, uint32_t x);

/*
 * Decode/encode canonical little-endian 64-byte integers.
 * fp505_decode() returns 1 on success and 0 if the input is >= p.
 * The in-memory fp505_t result is in Montgomery representation.
 */
int fp505_decode(fp505_t *r, const uint8_t in[SQISIGN_FP505_ENCODED_BYTES]);
void fp505_encode(uint8_t out[SQISIGN_FP505_ENCODED_BYTES], const fp505_t *a);


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
void fp505_loose_from_canonical(fp505_loose_t *r, const fp505_t *a);
void fp505_loose_add_bound4(fp505_loose_t *r, const fp505_loose_t *a,
                            const fp505_loose_t *b);
void fp505_loose_sub_bound4(fp505_loose_t *r, const fp505_loose_t *a,
                            const fp505_loose_t *b);
void fp505_loose_incomplete_reduce(fp505_loose_t *r,
                                   const fp505_loose_t *a);
void fp505_loose_complete_reduce(fp505_t *r, const fp505_loose_t *a);

/* Constant-schedule canonical field arithmetic. */
void fp505_add(fp505_t *r, const fp505_t *a, const fp505_t *b);
void fp505_sub(fp505_t *r, const fp505_t *a, const fp505_t *b);
void fp505_neg(fp505_t *r, const fp505_t *a);
void fp505_double(fp505_t *r, const fp505_t *a);
void fp505_mul(fp505_t *r, const fp505_t *a, const fp505_t *b);
void fp505_sqr(fp505_t *r, const fp505_t *a);

/* Returns 0 for input zero; otherwise computes a^(p-2). */
int fp505_inv(fp505_t *r, const fp505_t *a);

/*
 * Since p == 3 mod 4, sqrt(a)=a^((p+1)/4) when a is a square.
 * Returns 1 and a root on success, otherwise returns 0 and zeroes r.
 */
int fp505_sqrt(fp505_t *r, const fp505_t *a);

/* Returns 0 for zero, 1 for a nonzero square, and -1 for a nonsquare. */
int fp505_legendre(const fp505_t *a);

/* Comparison returns are 0 or 1. */
int fp505_iszero(const fp505_t *a);
int fp505_equal(const fp505_t *a, const fp505_t *b);

/* Explicitly wipe a field object. */
void fp505_secure_clear(fp505_t *a);

/* Extension-field arithmetic. */
void fp2_505_set_zero(fp2_505_t *r);
void fp2_505_set_one(fp2_505_t *r);
void fp2_505_copy(fp2_505_t *r, const fp2_505_t *a);
void fp2_505_add(fp2_505_t *r, const fp2_505_t *a, const fp2_505_t *b);
void fp2_505_sub(fp2_505_t *r, const fp2_505_t *a, const fp2_505_t *b);
void fp2_505_neg(fp2_505_t *r, const fp2_505_t *a);
void fp2_505_mul(fp2_505_t *r, const fp2_505_t *a, const fp2_505_t *b);
void fp2_505_sqr(fp2_505_t *r, const fp2_505_t *a);
int fp2_505_inv(fp2_505_t *r, const fp2_505_t *a);
int fp2_505_iszero(const fp2_505_t *a);
int fp2_505_equal(const fp2_505_t *a, const fp2_505_t *b);
void fp2_505_secure_clear(fp2_505_t *a);

/*
 * Two-way scalar batching for independent operations.  These routines do not
 * use NEON (which Cortex-M4 does not have); they expose independent dependency
 * chains so an M4 compiler/scheduler can interleave multiply-accumulates.
 */
void fp505_mul_2way(fp505_t out[2], const fp505_t a[2], const fp505_t b[2]);
void fp2_505_mul_2way(fp2_505_t out[2], const fp2_505_t a[2],
                      const fp2_505_t b[2]);

#ifdef __cplusplus
}
#endif

#endif /* SQISIGN_FP505_19X27_H */
