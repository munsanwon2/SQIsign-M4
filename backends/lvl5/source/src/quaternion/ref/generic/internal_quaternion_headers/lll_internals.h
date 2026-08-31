#ifndef LLL_INTERNALS_H
#define LLL_INTERNALS_H

/** @file
 *
 * @authors Sina Schaeffler
 *
 * @brief Declarations of functions only used for the LLL tets
 */

#include <quaternion.h>
#include "internal.h"

/** @internal
 * @ingroup quat_helpers
 * @defgroup lll_internal Functions only used for LLL or its tests
 * @{
 */

/** @internal
 * @ingroup lll_internal
 * @defgroup lll_params Parameters used by the L2 implementation (floats) and its tests (ints)
 * @{
 */

#define DELTABAR 0.995
#define DELTA_NUM 99
#define DELTA_DENOM 100

#define ETABAR 0.505
#define EPSILON_NUM 1
#define EPSILON_DENOM 100

#define PREC 64
/**
 * @}
 */

/**  @internal
 * @ingroup lll_internal
 * @defgroup quat_lll_verify_helpers Helper functions for lll verification in dimension 4
 * @{
 */

/** @brief Set ibq to parameters delta and eta = 1/2 + epsilon using L2 constants
 */
void quat_lll_set_ibq_parameters(ibq_t *delta, ibq_t *eta);

/** @brief Set an ibq vector to 4 given integer coefficients
 */
void ibq_vec_4_copy_ibz(ibq_vec_4_t *vec,
                        const ibz_t *coeff0,
                        const ibz_t *coeff1,
                        const ibz_t *coeff2,
                        const ibz_t *coeff3); // dim4, test/dim4

/** @brief Bilinear form vec00*vec10+vec01*vec11+q*vec02*vec12+q*vec03*vec13 for ibz_q
 */
void quat_lll_bilinear(ibq_t *b, const ibq_vec_4_t *vec0, const ibq_vec_4_t *vec1,
                       const ibz_t *q); // dim4, test/dim4

/** @brief Outputs the transposition of the orthogonalised matrix of mat (as fractions)
 *
 * For the bilinear form vec00*vec10+vec01*vec11+q*vec02*vec12+q*vec03*vec13
 */
void quat_lll_gram_schmidt_transposed_with_ibq(ibq_mat_4x4_t *orthogonalised_transposed,
                                               const ibz_mat_4x4_t *mat,
                                               const ibz_t *q); // dim4

/** @brief Verifies if mat is lll-reduced for parameter coeff and norm defined by q
 *
 * For the bilinear form vec00*vec10+vec01*vec11+q*vec02*vec12+q*vec03*vec13
 */
int quat_lll_verify(const ibz_mat_4x4_t *mat,
                    const ibq_t *delta,
                    const ibq_t *eta,
                    const quat_alg_t *alg); // test/lattice, test/dim4
                                            /** @}
                                             */

/** @internal
 * @ingroup lll_internal
 * @defgroup lll_internal_gram Internal LLL function
 * @{
 */

/** @brief In-place L2 reduction core function
 *
 * Given a lattice basis represented by the columns of a 4x4 matrix
 * and the Gram matrix of its bilinear form, L2-reduces the basis
 * in-place and updates the Gram matrix accordingly.
 *
 * Implements the L2 Algorithm of Nguyen-Stehlé, also known as fplll:
 * https://iacr.org/archive/eurocrypt2005/34940217/34940217.pdf
 *
 * Parameters are in lll/lll_internals.h
 *
 * @param G In/Output: Gram matrix of the lattice basis
 * @param basis In/Output: lattice basis
 */
/** Status-returning L2 core.  Integer updates are capacity checked and a zero
 * return means neither the mutated working Gram matrix nor basis may be used. */
int quat_lll_core_checked(ibz_mat_4x4_t *G, ibz_mat_4x4_t *basis);

/** Compatibility wrapper for test/benchmark code that supplies known-valid
 * bounded inputs.  Production callers must use quat_lll_core_checked. */
void quat_lll_core(ibz_mat_4x4_t *G, ibz_mat_4x4_t *basis);

/**
 * @brief LLL reduction on 4-dimensional lattice
 *
 * Implements the L2 Algorithm of Nguyen-Stehlé, also known as fplll:
 * https://iacr.org/archive/eurocrypt2005/34940217/34940217.pdf
 *
 * Parameters are in lll/lll_internals.h
 *
 * @param red Output: LLL reduced basis
 * @param lattice In/Output:  lattice with 4-dimensional basis
 * @param alg The quaternion algebra
 */
int quat_lattice_lll(ibz_mat_4x4_t *red, const quat_lattice_t *lattice, const quat_alg_t *alg);

/**
 * @}
 */

// end of lll_internal
/** @}
 */

/** ML2 — Modified L2 algorithm (NS09 Figure 9), generic d-rank input.
 *
 * Reduces d input generators (4-dimensional integer vectors). Writes up to
 * `out_capacity` reduced basis vectors into `output`, returns the rank.
 *
 * Used to replace `ibz_mat_4xn_hnf_mod_core` and `quat_lattice_hnf` in the
 * lattice operations (paper Issue 8: "replace only the HNF portion with ML2").
 * When `alg` is non-NULL, reduction uses the quaternion reduced-norm Gram
 * form diag(1,1,p,p), as required by compact MLLL.  A NULL `alg` selects a
 * Euclidean fallback only for algebra-independent span computations. */
int quat_ml2(ibz_vec_4_t *output,
             int out_capacity,
             const ibz_vec_4_t *input,
             int d,
             const quat_alg_t *alg);

/** Reducer signature shared by ML2 retry and fault-injection tests. */
typedef int (*quat_ml2_reducer_t)(ibz_vec_4_t *output,
                                  int out_capacity,
                                  const ibz_vec_4_t *input,
                                  int d,
                                  const quat_alg_t *alg);

/** ML2 could not be run without exceeding the active signed integer width. */
#define QUAT_ML2_ERR_PRECISION (-2)

/** Maximum number of ML2 attempts, including the unpermuted first attempt. */
#define QUAT_ML2_RETRY_MAX_ATTEMPTS 4

/**
 * Retry a full-rank ML2 reduction using deterministic generator permutations.
 *
 * The original order is tried first. If it returns rank four, its output is
 * published unchanged. Otherwise, the fixed schedule is one-step rotation,
 * reversal, then even-indexed generators followed by odd-indexed generators.
 * Only a rank-four result is published. `out_capacity` must be at least four;
 * output is unchanged for an invalid capacity or if every attempt fails.  The
 * native `quat_ml2_retry` entry point permits `output` to overlap `input`: it
 * copies the complete selected ordering before publishing a successful basis.
 */
int quat_ml2_retry_with_reducer(ibz_vec_4_t *output,
                                int out_capacity,
                                const ibz_vec_4_t *input,
                                int d,
                                const quat_alg_t *alg,
                                quat_ml2_reducer_t reducer);

/** Generic MLLL retry driver.  A valid lower-rank first result is published
 * without retries; negative failures try the remaining permutations.  This
 * entry point shares the same atomic profiling path as full-rank retry. */
int quat_ml2_mlll_with_reducer(ibz_vec_4_t *output,
                               int out_capacity,
                               const ibz_vec_4_t *input,
                               int d,
                               const quat_alg_t *alg,
                               quat_ml2_reducer_t reducer);

int quat_ml2_retry(ibz_vec_4_t *output,
                   int out_capacity,
                   const ibz_vec_4_t *input,
                   int d,
                   const quat_alg_t *alg);

/** Native generic-MLLL retry path.  A valid lower-rank first result is
 * published directly; negative failures use the deterministic permutations. */
int quat_ml2_mlll(ibz_vec_4_t *output,
                  int out_capacity,
                  const ibz_vec_4_t *input,
                  int d,
                  const quat_alg_t *alg);

#if defined(SQISIGN_ML2_PROFILE)
typedef struct {
    uint64_t inputs;
    uint64_t precision_rejected;
    uint64_t first_attempt_failures;
    uint64_t recovered[QUAT_ML2_RETRY_MAX_ATTEMPTS - 1];
    uint64_t exhausted;
    uint64_t underlying_attempts;
} quat_ml2_profile_dimension_t;

typedef struct {
    quat_ml2_profile_dimension_t d4;
    quat_ml2_profile_dimension_t d8;
    quat_ml2_profile_dimension_t d16;
} quat_ml2_profile_t;

void quat_ml2_profile_reset(void);
void quat_ml2_profile_get(quat_ml2_profile_t *profile);
#endif

#endif
