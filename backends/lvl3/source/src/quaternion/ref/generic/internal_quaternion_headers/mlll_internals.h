#ifndef MLLL_INTERNALS_H
#define MLLL_INTERNALS_H

/** @file
 *
 * @brief Declarations for Modified LLL (Pohst 1987) algorithm
 *
 * MLLL takes a generating set of a lattice (possibly linearly dependent)
 * and outputs an LLL-reduced basis. Used to replace HNF in ideal operations.
 */

#include <quaternion.h>
#include "lll_internals.h"

/** @brief Maximum number of generators MLLL can accept */
#define MLLL_MAX_GENERATORS 16

/** @brief Lattice dimension (quaternion algebra = 4) */
#define MLLL_DIM 4

/**
 * @brief Compact modified LLL algorithm backed by NS09 ML2
 *
 * Given g generating vectors of a lattice in Z^4, compute an LLL-reduced
 * basis of the lattice they span. Handles linear dependencies automatically.
 *
 * @param basis Output: 4x4 matrix whose columns are the LLL-reduced basis.
 *              Only the first `rank` columns are meaningful.
 * @param rank Output: rank of the lattice (number of basis vectors found)
 * @param generators Array of g vectors in Z^4 (each is ibz_vec_4_t)
 * @param g Number of generators (must be <= MLLL_MAX_GENERATORS)
 * @param alg The quaternion algebra; selects the reduced-norm Gram form
 */
void quat_mlll(ibz_mat_4x4_t *basis,
               int *rank,
               const ibz_vec_4_t *generators,
               int g,
               const quat_alg_t *alg);

/**
 * @brief Lattice multiplication using MLLL instead of HNF
 *
 * Replaces quat_lattice_mul with MLLL-based approach.
 * Computes res = lat1 * lat2 using MLLL on 16 product generators.
 *
 * @param res Output: product lattice
 * @param lat1 First lattice
 * @param lat2 Second lattice
 * @param alg The quaternion algebra
 * @returns 1 on success, 0 with res unchanged on invalid input, precision
 *          exhaustion, or incomplete reduction.  Aliasing res with either
 *          input is supported.
 */
int quat_lattice_mul_mlll(quat_lattice_t *res,
                           const quat_lattice_t *lat1,
                           const quat_lattice_t *lat2,
                           const quat_alg_t *alg);

/**
 * @brief Lattice addition using ML2 instead of HNF
 *
 * Replaces quat_lattice_add with MLLL-based approach.
 *
 * @param res Output: sum lattice
 * @param lat1 First lattice
 * @param lat2 Second lattice
 * @param alg The quaternion algebra
 * @param reducer Reducer implementation; production passes quat_ml2_retry
 * @returns 1 on success, 0 if the reducer did not return a rank-four basis.
 *          On failure, res is left unchanged.
 */
int quat_lattice_add_mlll_with_reducer(quat_lattice_t *res,
                                       const quat_lattice_t *lat1,
                                       const quat_lattice_t *lat2,
                                       const quat_alg_t *alg,
                                       quat_ml2_reducer_t reducer);

int quat_lattice_add_mlll(quat_lattice_t *res,
                          const quat_lattice_t *lat1,
                          const quat_lattice_t *lat2,
                          const quat_alg_t *alg);

/** @brief Paper Algorithm 3, CompactIdealIntersection. */
int quat_lattice_intersect_mlll(quat_lattice_t *res,
                                 ibz_t *intersection_norm,
                                 const quat_lattice_t *lat1,
                                 const ibz_t *norm1,
                                 const quat_lattice_t *lat2,
                                 const ibz_t *norm2,
                                 const quat_alg_t *alg);

#endif
