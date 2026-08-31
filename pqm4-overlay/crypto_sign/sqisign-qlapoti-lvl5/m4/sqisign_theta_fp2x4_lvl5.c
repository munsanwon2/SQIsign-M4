#include "sqisign_theta_fp2x4_lvl5.h"

#include <stddef.h>

void theta_fp2x4_lvl5_set_zero(theta_fp2x4_lvl5_t *r)
{
    size_t i;
    for (i = 0; i < 4u; i++) {
        fp2_505_set_zero(&r->c[i]);
    }
}

void theta_fp2x4_lvl5_copy(theta_fp2x4_lvl5_t *r, const theta_fp2x4_lvl5_t *a)
{
    *r = *a;
}

void theta_fp2x4_lvl5_hadamard(theta_fp2x4_lvl5_t *r, const theta_fp2x4_lvl5_t *a)
{
    /*
     * Holding the four pairwise sums/differences is sufficient for exact
     * in-place operation.  Avoiding a second four-coordinate output object
     * saves 288 bytes of stack on the 19x27 Level-5 representation.
     */
    fp2_505_t s01, d01, s23, d23;

    fp2_505_add(&s01, &a->c[0], &a->c[1]);
    fp2_505_sub(&d01, &a->c[0], &a->c[1]);
    fp2_505_add(&s23, &a->c[2], &a->c[3]);
    fp2_505_sub(&d23, &a->c[2], &a->c[3]);

    fp2_505_add(&r->c[0], &s01, &s23);
    fp2_505_add(&r->c[1], &d01, &d23);
    fp2_505_sub(&r->c[2], &s01, &s23);
    fp2_505_sub(&r->c[3], &d01, &d23);

    fp2_505_secure_clear(&s01);
    fp2_505_secure_clear(&d01);
    fp2_505_secure_clear(&s23);
    fp2_505_secure_clear(&d23);
}

void theta_fp2x4_lvl5_square(theta_fp2x4_lvl5_t *r, const theta_fp2x4_lvl5_t *a)
{
    size_t i;
    /* Each output coordinate depends only on the matching input coordinate. */
    for (i = 0; i < 4u; i++) {
        fp2_505_sqr(&r->c[i], &a->c[i]);
    }
}

void theta_fp2x4_lvl5_mul(theta_fp2x4_lvl5_t *r,
                     const theta_fp2x4_lvl5_t *a,
                     const theta_fp2x4_lvl5_t *b)
{
    size_t i;
    /* Exact r==a and r==b aliasing are safe coordinate by coordinate. */
    for (i = 0; i < 4u; i++) {
        fp2_505_mul(&r->c[i], &a->c[i], &b->c[i]);
    }
}

void theta_fp2x4_lvl5_hadamard_square(theta_fp2x4_lvl5_t *r,
                                 const theta_fp2x4_lvl5_t *a)
{
    theta_fp2x4_lvl5_t t;
    theta_fp2x4_lvl5_hadamard(&t, a);
    theta_fp2x4_lvl5_square(r, &t);
    theta_fp2x4_lvl5_secure_clear(&t);
}

int theta_fp2x4_lvl5_equal(const theta_fp2x4_lvl5_t *a, const theta_fp2x4_lvl5_t *b)
{
    int equal = 1;
    size_t i;
    for (i = 0; i < 4u; i++) {
        equal &= fp2_505_equal(&a->c[i], &b->c[i]);
    }
    return equal;
}

void theta_fp2x4_lvl5_secure_clear(theta_fp2x4_lvl5_t *a)
{
    size_t i;
    for (i = 0; i < 4u; i++) {
        fp2_505_secure_clear(&a->c[i]);
    }
}
