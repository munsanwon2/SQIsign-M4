#ifndef SQISIGN_FP383_14X28_H
#define SQISIGN_FP383_14X28_H

/*
 * Scalar ARM Cortex-M4 field backend for the SQIsign level-3 prime
 *
 *     p = 65 * 2^376 - 1.
 *
 * The representation follows the unsaturated 14-limb layout studied in
 * "SQIsign on ARM": 13 28-bit limbs and one 19-bit limb.  Values held in
 * fp383_t are canonical Montgomery residues in [0,p), with Montgomery radix
 * R = 2^(14*28) = 2^392.  The implementation is freestanding C11 and requires
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

#define SQISIGN_FP383_LIMBS 14u
#define SQISIGN_FP383_RADIX_BITS 28u
#define SQISIGN_FP383_ENCODED_BYTES 48u
#define SQISIGN_FP383_MASK UINT32_C(0x0fffffff)

/* Canonical Montgomery-domain field element. */
typedef struct {
    uint32_t limb[SQISIGN_FP383_LIMBS];
} fp383_t;

/*
 * Non-canonical unsaturated value used for statically scheduled lazy
 * reduction.  The arithmetic layer, not the type, tracks the proven bound.
 * Safe inputs to the helpers below have lower limbs < 8*2^28 and top limb
 * < 32*(65*2^12); see the function contracts.
 */
typedef struct {
    uint32_t limb[SQISIGN_FP383_LIMBS];
} fp383_loose_t;

/* F_{p^2} = F_p[i]/(i^2+1). */
typedef struct {
    fp383_t re;
    fp383_t im;
} fp2_383_t;

/* Constants, all in the documented domain. */
extern const fp383_t fp383_zero;
extern const fp383_t fp383_one;       /* Montgomery representation of 1. */
extern const fp383_t fp383_modulus;   /* Standard integer representation. */
extern const fp383_t fp383_r2;        /* R^2 mod p, standard limbs. */

/* Base-field construction and conversion. */
void fp383_set_zero(fp383_t *r);
void fp383_set_one(fp383_t *r);
void fp383_copy(fp383_t *r, const fp383_t *a);
void fp383_from_u32(fp383_t *r, uint32_t x);

/*
 * Decode/encode canonical little-endian 48-byte integers.
 * fp383_decode() returns 1 on success and 0 if the input is >= p.
 * The in-memory fp383_t result is in Montgomery representation.
 */
int fp383_decode(fp383_t *r, const uint8_t in[SQISIGN_FP383_ENCODED_BYTES]);
void fp383_encode(uint8_t out[SQISIGN_FP383_ENCODED_BYTES], const fp383_t *a);


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
void fp383_loose_from_canonical(fp383_loose_t *r, const fp383_t *a);
void fp383_loose_add_bound4(fp383_loose_t *r, const fp383_loose_t *a,
                            const fp383_loose_t *b);
void fp383_loose_sub_bound4(fp383_loose_t *r, const fp383_loose_t *a,
                            const fp383_loose_t *b);
void fp383_loose_incomplete_reduce(fp383_loose_t *r,
                                   const fp383_loose_t *a);
void fp383_loose_complete_reduce(fp383_t *r, const fp383_loose_t *a);

/* Constant-schedule canonical field arithmetic. */
void fp383_add(fp383_t *r, const fp383_t *a, const fp383_t *b);
void fp383_sub(fp383_t *r, const fp383_t *a, const fp383_t *b);
void fp383_neg(fp383_t *r, const fp383_t *a);
void fp383_double(fp383_t *r, const fp383_t *a);
void fp383_mul(fp383_t *r, const fp383_t *a, const fp383_t *b);
void fp383_sqr(fp383_t *r, const fp383_t *a);

/* Returns 0 for input zero; otherwise computes a^(p-2). */
int fp383_inv(fp383_t *r, const fp383_t *a);

/*
 * Since p == 3 mod 4, sqrt(a)=a^((p+1)/4) when a is a square.
 * Returns 1 and a root on success, otherwise returns 0 and zeroes r.
 */
int fp383_sqrt(fp383_t *r, const fp383_t *a);

/* Returns 0 for zero, 1 for a nonzero square, and -1 for a nonsquare. */
int fp383_legendre(const fp383_t *a);

/* Comparison returns are 0 or 1. */
int fp383_iszero(const fp383_t *a);
int fp383_equal(const fp383_t *a, const fp383_t *b);

/* Explicitly wipe a field object. */
void fp383_secure_clear(fp383_t *a);

/* Extension-field arithmetic. */
void fp2_383_set_zero(fp2_383_t *r);
void fp2_383_set_one(fp2_383_t *r);
void fp2_383_copy(fp2_383_t *r, const fp2_383_t *a);
void fp2_383_add(fp2_383_t *r, const fp2_383_t *a, const fp2_383_t *b);
void fp2_383_sub(fp2_383_t *r, const fp2_383_t *a, const fp2_383_t *b);
void fp2_383_neg(fp2_383_t *r, const fp2_383_t *a);
void fp2_383_mul(fp2_383_t *r, const fp2_383_t *a, const fp2_383_t *b);
void fp2_383_sqr(fp2_383_t *r, const fp2_383_t *a);
int fp2_383_inv(fp2_383_t *r, const fp2_383_t *a);
int fp2_383_iszero(const fp2_383_t *a);
int fp2_383_equal(const fp2_383_t *a, const fp2_383_t *b);
void fp2_383_secure_clear(fp2_383_t *a);

/*
 * Two-way scalar batching for independent operations.  These routines do not
 * use NEON (which Cortex-M4 does not have); they expose independent dependency
 * chains so an M4 compiler/scheduler can interleave multiply-accumulates.
 */
void fp383_mul_2way(fp383_t out[2], const fp383_t a[2], const fp383_t b[2]);
void fp2_383_mul_2way(fp2_383_t out[2], const fp2_383_t a[2],
                      const fp2_383_t b[2]);

#ifdef __cplusplus
}
#endif

#endif /* SQISIGN_FP383_14X28_H */
