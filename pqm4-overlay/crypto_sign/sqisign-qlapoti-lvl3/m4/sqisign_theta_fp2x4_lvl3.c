#include "sqisign_theta_fp2x4_lvl3.h"

#include <stddef.h>

void theta_fp2x4_lvl3_set_zero(theta_fp2x4_lvl3_t *r)
{
    size_t i;
    for (i = 0; i < 4u; i++) {
        fp2_383_set_zero(&r->c[i]);
    }
}

void theta_fp2x4_lvl3_copy(theta_fp2x4_lvl3_t *r, const theta_fp2x4_lvl3_t *a)
{
    *r = *a;
}

void theta_fp2x4_lvl3_hadamard(theta_fp2x4_lvl3_t *r, const theta_fp2x4_lvl3_t *a)
{
    /*
     * Holding the four pairwise sums/differences is sufficient for exact
     * in-place operation.  Avoiding a second four-coordinate output object
     * saves 288 bytes of stack on the 14x28 Level-3 representation.
     */
    fp2_383_t s01, d01, s23, d23;

    fp2_383_add(&s01, &a->c[0], &a->c[1]);
    fp2_383_sub(&d01, &a->c[0], &a->c[1]);
    fp2_383_add(&s23, &a->c[2], &a->c[3]);
    fp2_383_sub(&d23, &a->c[2], &a->c[3]);

    fp2_383_add(&r->c[0], &s01, &s23);
    fp2_383_add(&r->c[1], &d01, &d23);
    fp2_383_sub(&r->c[2], &s01, &s23);
    fp2_383_sub(&r->c[3], &d01, &d23);

    fp2_383_secure_clear(&s01);
    fp2_383_secure_clear(&d01);
    fp2_383_secure_clear(&s23);
    fp2_383_secure_clear(&d23);
}

void theta_fp2x4_lvl3_square(theta_fp2x4_lvl3_t *r, const theta_fp2x4_lvl3_t *a)
{
    size_t i;
    /* Each output coordinate depends only on the matching input coordinate. */
    for (i = 0; i < 4u; i++) {
        fp2_383_sqr(&r->c[i], &a->c[i]);
    }
}

void theta_fp2x4_lvl3_mul(theta_fp2x4_lvl3_t *r,
                     const theta_fp2x4_lvl3_t *a,
                     const theta_fp2x4_lvl3_t *b)
{
    size_t i;
    /* Exact r==a and r==b aliasing are safe coordinate by coordinate. */
    for (i = 0; i < 4u; i++) {
        fp2_383_mul(&r->c[i], &a->c[i], &b->c[i]);
    }
}

void theta_fp2x4_lvl3_hadamard_square(theta_fp2x4_lvl3_t *r,
                                 const theta_fp2x4_lvl3_t *a)
{
    theta_fp2x4_lvl3_t t;
    theta_fp2x4_lvl3_hadamard(&t, a);
    theta_fp2x4_lvl3_square(r, &t);
    theta_fp2x4_lvl3_secure_clear(&t);
}

int theta_fp2x4_lvl3_equal(const theta_fp2x4_lvl3_t *a, const theta_fp2x4_lvl3_t *b)
{
    int equal = 1;
    size_t i;
    for (i = 0; i < 4u; i++) {
        equal &= fp2_383_equal(&a->c[i], &b->c[i]);
    }
    return equal;
}

void theta_fp2x4_lvl3_secure_clear(theta_fp2x4_lvl3_t *a)
{
    size_t i;
    for (i = 0; i < 4u; i++) {
        fp2_383_secure_clear(&a->c[i]);
    }
}
