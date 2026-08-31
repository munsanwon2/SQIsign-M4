#ifndef SQISIGN_THETA_FP2X4_H
#define SQISIGN_THETA_FP2X4_H

/*
 * Scalar four-coordinate theta kernels for ARM Cortex-M4.
 *
 * The application-class ARM implementation maps the four independent theta
 * coordinates to NEON lanes.  Cortex-M4 has no NEON, so this module preserves
 * the same four-coordinate data flow while executing each F_{p^2} operation
 * with the scalar 9x29 backend.  It is a reusable arithmetic primitive, not a
 * complete (2,2)-isogeny-chain implementation.
 */

#include "sqisign_fp251_9x29.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    fp2_251_t c[4];
} theta_fp2x4_t;

void theta_fp2x4_set_zero(theta_fp2x4_t *r);
void theta_fp2x4_copy(theta_fp2x4_t *r, const theta_fp2x4_t *a);

/*
 * H(x0,x1,x2,x3) =
 *   (x0+x1+x2+x3,
 *    x0-x1+x2-x3,
 *    x0+x1-x2-x3,
 *    x0-x1-x2+x3).
 *
 * Input/output aliasing is supported.
 */
void theta_fp2x4_hadamard(theta_fp2x4_t *r, const theta_fp2x4_t *a);

/* Coordinate-wise operations used by theta doubling/isogeny formulae. */
void theta_fp2x4_square(theta_fp2x4_t *r, const theta_fp2x4_t *a);
void theta_fp2x4_mul(theta_fp2x4_t *r,
                     const theta_fp2x4_t *a,
                     const theta_fp2x4_t *b);
void theta_fp2x4_hadamard_square(theta_fp2x4_t *r,
                                 const theta_fp2x4_t *a);

int theta_fp2x4_equal(const theta_fp2x4_t *a, const theta_fp2x4_t *b);
void theta_fp2x4_secure_clear(theta_fp2x4_t *a);

#ifdef __cplusplus
}
#endif

#endif /* SQISIGN_THETA_FP2X4_H */
