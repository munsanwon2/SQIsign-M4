#include <quaternion.h>
#include "internal.h"

// internal helper functions
void
ibz_mat_4x4_mul(ibz_mat_4x4_t *res, const ibz_mat_4x4_t *a, const ibz_mat_4x4_t *b)
{
    ibz_mat_4x4_t mat;
    ibz_t prod;
    ibz_init(&prod);
    ibz_mat_4x4_init(&mat);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&(mat[i][j]), 0);
            for (int k = 0; k < 4; k++) {
                ibz_mul(&prod, &((*a)[i][k]), &((*b)[k][j]));
                ibz_add(&(mat[i][j]), &(mat[i][j]), &prod);
            }
        }
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_copy(&((*res)[i][j]), &(mat[i][j]));
        }
    }
    ibz_mat_4x4_finalize(&mat);
    ibz_finalize(&prod);
}

// helper functions for lattices
void
ibz_vec_4_set(ibz_vec_4_t *vec, int32_t coord0, int32_t coord1, int32_t coord2, int32_t coord3)
{
    ibz_set(&((*vec)[0]), coord0);
    ibz_set(&((*vec)[1]), coord1);
    ibz_set(&((*vec)[2]), coord2);
    ibz_set(&((*vec)[3]), coord3);
}

void
ibz_vec_4_copy(ibz_vec_4_t *new, const ibz_vec_4_t *vec)
{
    for (int i = 0; i < 4; i++) {
        ibz_copy(&((*new)[i]), &((*vec)[i]));
    }
}

void
ibz_vec_4_copy_ibz(ibz_vec_4_t *res, const ibz_t *coord0, const ibz_t *coord1, const ibz_t *coord2, const ibz_t *coord3)
{
    ibz_copy(&((*res)[0]), coord0);
    ibz_copy(&((*res)[1]), coord1);
    ibz_copy(&((*res)[2]), coord2);
    ibz_copy(&((*res)[3]), coord3);
}

void
ibz_vec_4_content(ibz_t *content, const ibz_vec_4_t *v)
{
    ibz_gcd(content, &((*v)[0]), &((*v)[1]));
    ibz_gcd(content, &((*v)[2]), content);
    ibz_gcd(content, &((*v)[3]), content);
}

void
ibz_vec_4_negate(ibz_vec_4_t *neg, const ibz_vec_4_t *vec)
{
    for (int i = 0; i < 4; i++) {
        ibz_neg(&((*neg)[i]), &((*vec)[i]));
    }
}

void
ibz_vec_4_add(ibz_vec_4_t *res, const ibz_vec_4_t *a, const ibz_vec_4_t *b)
{
    ibz_add(&((*res)[0]), &((*a)[0]), &((*b)[0]));
    ibz_add(&((*res)[1]), &((*a)[1]), &((*b)[1]));
    ibz_add(&((*res)[2]), &((*a)[2]), &((*b)[2]));
    ibz_add(&((*res)[3]), &((*a)[3]), &((*b)[3]));
}

void
ibz_vec_4_sub(ibz_vec_4_t *res, const ibz_vec_4_t *a, const ibz_vec_4_t *b)
{
    ibz_sub(&((*res)[0]), &((*a)[0]), &((*b)[0]));
    ibz_sub(&((*res)[1]), &((*a)[1]), &((*b)[1]));
    ibz_sub(&((*res)[2]), &((*a)[2]), &((*b)[2]));
    ibz_sub(&((*res)[3]), &((*a)[3]), &((*b)[3]));
}

int
ibz_vec_4_is_zero(const ibz_vec_4_t *x)
{
    int res = 1;
    for (int i = 0; i < 4; i++) {
        res &= ibz_is_zero(&((*x)[i]));
    }
    return (res);
}

void
ibz_vec_4_linear_combination(ibz_vec_4_t *lc,
                             const ibz_t *coeff_a,
                             const ibz_vec_4_t *vec_a,
                             const ibz_t *coeff_b,
                             const ibz_vec_4_t *vec_b)
{
    ibz_t prod;
    ibz_vec_4_t sums;
    ibz_vec_4_init(&sums);
    ibz_init(&prod);
    for (int i = 0; i < 4; i++) {
        ibz_mul(&(sums[i]), coeff_a, &((*vec_a)[i]));
        ibz_mul(&prod, coeff_b, &((*vec_b)[i]));
        ibz_add(&(sums[i]), &(sums[i]), &prod);
    }
    for (int i = 0; i < 4; i++) {
        ibz_copy(&((*lc)[i]), &(sums[i]));
    }
    ibz_finalize(&prod);
    ibz_vec_4_finalize(&sums);
}

void
ibz_vec_4_scalar_mul(ibz_vec_4_t *prod, const ibz_t *scalar, const ibz_vec_4_t *vec)
{
    for (int i = 0; i < 4; i++) {
        ibz_mul(&((*prod)[i]), &((*vec)[i]), scalar);
    }
}

int
ibz_vec_4_scalar_div(ibz_vec_4_t *quot, const ibz_t *scalar, const ibz_vec_4_t *vec)
{
    int res = 1;
    ibz_t r;
    ibz_init(&r);
    for (int i = 0; i < 4; i++) {
        ibz_div(&((*quot)[i]), &r, &((*vec)[i]), scalar);
        res = res && ibz_is_zero(&r);
    }
    ibz_finalize(&r);
    return (res);
}

void
ibz_mat_4x4_copy(ibz_mat_4x4_t *new, const ibz_mat_4x4_t *mat)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_copy(&((*new)[i][j]), &((*mat)[i][j]));
        }
    }
}

void
ibz_mat_4x4_negate(ibz_mat_4x4_t *neg, const ibz_mat_4x4_t *mat)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_neg(&((*neg)[i][j]), &((*mat)[i][j]));
        }
    }
}

void
ibz_mat_4x4_transpose(ibz_mat_4x4_t *transposed, const ibz_mat_4x4_t *mat)
{
    ibz_mat_4x4_t work;
    ibz_mat_4x4_init(&work);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_copy(&(work[i][j]), &((*mat)[j][i]));
        }
    }
    ibz_mat_4x4_copy(transposed, &work);
    ibz_mat_4x4_finalize(&work);
}

void
ibz_mat_4x4_zero(ibz_mat_4x4_t *zero)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&((*zero)[i][j]), 0);
        }
    }
}

void
ibz_mat_4x4_identity(ibz_mat_4x4_t *id)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_set(&((*id)[i][j]), 0);
        }
        ibz_set(&((*id)[i][i]), 1);
    }
}

int
ibz_mat_4x4_is_identity(const ibz_mat_4x4_t *mat)
{
    int res = 1;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            res = res && ibz_is_one(&((*mat)[i][j])) == (i == j);
        }
    }
    return (res);
}

int
ibz_mat_4x4_equal(const ibz_mat_4x4_t *mat1, const ibz_mat_4x4_t *mat2)
{
    int res = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            res = res | ibz_cmp(&((*mat1)[i][j]), &((*mat2)[i][j]));
        }
    }
    return (!res);
}

void
ibz_mat_4x4_scalar_mul(ibz_mat_4x4_t *prod, const ibz_t *scalar, const ibz_mat_4x4_t *mat)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_mul(&((*prod)[i][j]), &((*mat)[i][j]), scalar);
        }
    }
}

void
ibz_mat_4x4_gcd(ibz_t *gcd, const ibz_mat_4x4_t *mat)
{
    ibz_t d;
    ibz_init(&d);
    ibz_copy(&d, &((*mat)[0][0]));
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_gcd(&d, &d, &((*mat)[i][j]));
        }
    }
    ibz_copy(gcd, &d);
    ibz_finalize(&d);
}

int
ibz_mat_4x4_scalar_div(ibz_mat_4x4_t *quot, const ibz_t *scalar, const ibz_mat_4x4_t *mat)
{
    int res = 1;
    ibz_t r;
    ibz_init(&r);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_div(&((*quot)[i][j]), &r, &((*mat)[i][j]), scalar);
            res = res && ibz_is_zero(&r);
        }
    }
    ibz_finalize(&r);
    return (res);
}

// 4x4 inversion helper functions
void
ibz_inv_dim4_make_coeff_pmp(ibz_t *coeff,
                            const ibz_t *a1,
                            const ibz_t *a2,
                            const ibz_t *b1,
                            const ibz_t *b2,
                            const ibz_t *c1,
                            const ibz_t *c2)
{
    ibz_t prod, sum;
    ibz_init(&prod);
    ibz_init(&sum);
    ibz_mul(&sum, a1, a2);
    ibz_mul(&prod, b1, b2);
    ibz_sub(&sum, &sum, &prod);
    ibz_mul(&prod, c1, c2);
    ibz_add(coeff, &sum, &prod);
    ibz_finalize(&prod);
    ibz_finalize(&sum);
}

void
ibz_inv_dim4_make_coeff_mpm(ibz_t *coeff,
                            const ibz_t *a1,
                            const ibz_t *a2,
                            const ibz_t *b1,
                            const ibz_t *b2,
                            const ibz_t *c1,
                            const ibz_t *c2)
{
    ibz_t prod, sum;
    ibz_init(&prod);
    ibz_init(&sum);
    ibz_mul(&sum, b1, b2);
    ibz_mul(&prod, a1, a2);
    ibz_sub(&sum, &sum, &prod);
    ibz_mul(&prod, c1, c2);
    ibz_sub(coeff, &sum, &prod);
    ibz_finalize(&prod);
    ibz_finalize(&sum);
}

static int
minor_2x2(ibz_t *det, const ibz_t *a, const ibz_t *b, const ibz_t *c, const ibz_t *d)
{
    ibz_t ad, bc;
    return ibz_mul_checked(&ad, a, d) && ibz_mul_checked(&bc, b, c) &&
           ibz_sub_checked(det, &ad, &bc);
}

/* Laplace expansion along the first two rows. Compute the adjugate only
 * when the caller needs it. None of the signed minors may wrap. */
int
ibz_mat_4x4_inv_with_det_as_denom(ibz_mat_4x4_t *inv, ibz_t *det, const ibz_mat_4x4_t *mat)
{
    ibz_t top[6], bottom[6], value, product;
    int index = 0;
    for (int a = 0; a < 4; a++) {
        for (int b = a + 1; b < 4; b++, index++) {
            if (!minor_2x2(&top[index], &(*mat)[0][a], &(*mat)[0][b], &(*mat)[1][a], &(*mat)[1][b]) ||
                !minor_2x2(&bottom[index], &(*mat)[2][a], &(*mat)[2][b], &(*mat)[3][a], &(*mat)[3][b]))
                return 0;
        }
    }
    ibz_set(&value, 0);
    for (int i = 0; i < 6; i++) {
        if (!ibz_mul_checked(&product, &top[i], &bottom[5 - i]))
            return 0;
        if (!(i == 1 || i == 4 ? ibz_sub_checked(&value, &value, &product)
                               : ibz_add_checked(&value, &value, &product)))
            return 0;
    }

    if (inv != NULL) {
        ibz_mat_4x4_t adjugate = { 0 };
        if (!ibz_is_zero(&value)) {
            for (int row = 0; row < 4; row++) {
                const int paired_row = row ^ 1;
                const ibz_t *minors = row < 2 ? bottom : top;
                for (int col = 0; col < 4; col++) {
                    const int cols[3] = {col != 0 ? 0 : 1, col < 2 ? 2 : 1, col != 3 ? 3 : 2};
                    const int minor_index[3] = {6 - col - (col == 0), 4 - col - (col == 1),
                                                3 - col - (col == 1) - (col == 2)};
                    ibz_t *entry = &adjugate[col][row];
                    for (int i = 0; i < 3; i++) {
                        if (!ibz_mul_checked(&product, &(*mat)[paired_row][cols[i]], &minors[minor_index[i]]) ||
                            !((row + col + i) % 2 ? ibz_sub_checked(entry, entry, &product)
                                                  : ibz_add_checked(entry, entry, &product)))
                            return 0;
                    }
                }
            }
        }
        ibz_mat_4x4_copy(inv, &adjugate);
    }
    if (det != NULL)
        ibz_copy(det, &value);
    return !ibz_is_zero(&value);
}

int
ibz_mat_4x4_eval_checked(ibz_vec_4_t *res, const ibz_mat_4x4_t *mat, const ibz_vec_4_t *vec)
{
    ibz_vec_4_t candidate = { 0 };
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_t product;
            if (!ibz_mul_checked(&product, &(*mat)[i][j], &(*vec)[j]) ||
                !ibz_add_checked(&candidate[i], &candidate[i], &product))
                return 0;
        }
    }
    ibz_vec_4_copy(res, &candidate);
    return 1;
}

// matrix evaluation

void
ibz_mat_4x4_eval(ibz_vec_4_t *res, const ibz_mat_4x4_t *mat, const ibz_vec_4_t *vec)
{
    ibz_vec_4_t sum;
    ibz_t prod;
    ibz_init(&prod);
    ibz_vec_4_init(&sum);
    // assume initialization to 0
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_mul(&prod, &(*mat)[i][j], &(*vec)[j]);
            ibz_add(&(sum[i]), &(sum[i]), &prod);
        }
    }
    ibz_vec_4_copy(res, &sum);
    ibz_finalize(&prod);
    ibz_vec_4_finalize(&sum);
}

void
ibz_mat_4x4_eval_t(ibz_vec_4_t *res, const ibz_vec_4_t *vec, const ibz_mat_4x4_t *mat)
{
    ibz_vec_4_t sum;
    ibz_t prod;
    ibz_init(&prod);
    ibz_vec_4_init(&sum);
    // assume initialization to 0
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            ibz_mul(&prod, &(*mat)[j][i], &(*vec)[j]);
            ibz_add(&(sum[i]), &(sum[i]), &prod);
        }
    }
    ibz_vec_4_copy(res, &sum);
    ibz_finalize(&prod);
    ibz_vec_4_finalize(&sum);
}

// quadratic forms

void
quat_qf_eval(ibz_t *res, const ibz_mat_4x4_t *qf, const ibz_vec_4_t *coord)
{
    ibz_vec_4_t sum;
    ibz_t prod;
    ibz_init(&prod);
    ibz_vec_4_init(&sum);
    ibz_mat_4x4_eval(&sum, qf, coord);
    for (int i = 0; i < 4; i++) {
        ibz_mul(&prod, &(sum[i]), &(*coord)[i]);
        if (i > 0) {
            ibz_add(&(sum[0]), &(sum[0]), &prod);
        } else {
            ibz_copy(&sum[0], &prod);
        }
    }
    ibz_copy(res, &sum[0]);
    ibz_finalize(&prod);
    ibz_vec_4_finalize(&sum);
}
