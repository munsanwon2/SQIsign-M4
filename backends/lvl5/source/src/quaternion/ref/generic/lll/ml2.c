/* Section 4 ML2 (NS09 Figure 9).
 * Vectors and the Gram matrix are exact; Cholesky and Babai decisions use
 * a 53-bit significand with an independent exponent. A zero prefix accounts
 * for dependent generators without changing the Z-span.
 */

#include <fenv.h>
#include <quaternion.h>
#include "lll_internals.h"
#include "internal.h"
#include "dpe.h"

#if defined(TARGET_ARM)
#include <sqisign_m4_workspace.h>
#endif

#if defined(SQISIGN_ML2_PROFILE)
#include <stdatomic.h>
#endif

/* Two rank-four bases give 16 ordered product generators. */
#define ML2_MAX_D 16

#define ML2_MU_COUNT ((ML2_MAX_D * (ML2_MAX_D - 1)) / 2)
#define ML2_MU(ST, I, J) \
    ((ST)->mu[((I) * ((I) - 1)) / 2 + (J)])

/* ========== state ========== */

typedef struct {
    int d;
    int zeta;
    int failure;
    /* NULL selects the Euclidean fallback used by algebra-independent
     * span-only callers.  Compact quaternion algorithms always pass alg and
     * therefore use the reduced-norm form diag(1, 1, p, p). */
    const ibz_t *p;
    ibz_t coefficients[ML2_MAX_D - 1];
    /* Only the current CFA row and the persistent diagonal of r are ever
     * read.  Historical strict-lower mu rows are packed densely. */
    dpe_t r_diag[ML2_MAX_D];
    dpe_t r_row[ML2_MAX_D];
    dpe_t mu[ML2_MU_COUNT];
    dpe_t s[ML2_MAX_D + 1];
    dpe_t delta_bar;
#if defined(TARGET_ARM)
    ibz_vec_4_t b[];
#else
    ibz_vec_4_t b[ML2_MAX_D];
#endif
} ml2_state_t;

#if defined(TARGET_ARM)
/* CCM scratch includes the complete pre-update Babai coefficient vector. */
_Static_assert(sizeof(ml2_state_t) + ML2_MAX_D * sizeof(ibz_vec_4_t) <= 46344u,
               "ML2 state exceeds the Cortex-M4 CCM scratch budget");
#endif

static void
ml2_state_init(ml2_state_t *st, int d, const quat_alg_t *alg)
{
    st->d = d;
    st->zeta = 0;
    st->failure = 0;
    st->p = alg == NULL ? NULL : &alg->p;
    for (int i = 0; i < d; i++) {
        ibz_vec_4_init(&st->b[i]);
        dpe_init(st->r_diag[i]);
        dpe_init(st->r_row[i]);
    }
    for (int i = 0; i < d * (d - 1) / 2; i++)
        dpe_init(st->mu[i]);
    for (int i = 0; i <= d; i++)
        dpe_init(st->s[i]);
    dpe_init(st->delta_bar);
    dpe_set_d(st->delta_bar, DELTABAR);
}

static void
ml2_state_finalize(ml2_state_t *st)
{
    int d = st->d;
    for (int i = 0; i < d; i++) {
        ibz_vec_4_finalize(&st->b[i]);
        dpe_clear(st->r_diag[i]);
        dpe_clear(st->r_row[i]);
    }
    for (int i = 0; i < d * (d - 1) / 2; i++)
        dpe_clear(st->mu[i]);
    for (int i = 0; i <= d; i++)
        dpe_clear(st->s[i]);
    dpe_clear(st->delta_bar);
}

/* ML2's active Gram--Schmidt block must start with a positive diagonal.
 * Move exact zero generators into the skipped prefix before any floating
 * division is attempted.  The stable rotation keeps the relative order of
 * nonzero generators, which also keeps retry permutations deterministic. */
static void
ml2_partition_initial_zeros(ml2_state_t *st)
{
    int zeta = 0;

    for (int i = 0; i < st->d; i++) {
        if (!ibz_vec_4_is_zero(&st->b[i]))
            continue;
        for (int j = i; j > zeta; j--)
            for (int coordinate = 0; coordinate < 4; coordinate++)
                ibz_swap(&st->b[j][coordinate],
                         &st->b[j - 1][coordinate]);
        zeta++;
    }
    st->zeta = zeta;
}

/* Half the trace pairing is diag(1,1,p,p) in ambient quaternion
 * coordinates. With no algebra the span-only metric is Euclidean. */
static int
ml2_gram_entry(ibz_t *value, const ml2_state_t *st, int i, int j)
{
    ibz_t term;
    ibz_set(value, 0);
    for (int c = 0; c < 4; c++) {
        if (!ibz_mul_checked(&term, &st->b[i][c], &st->b[j][c]) ||
            (c >= 2 && st->p != NULL && !ibz_mul_checked(&term, &term, st->p)) ||
            !ibz_add_checked(value, value, &term))
            return 0;
    }
    return 1;
}

static int
ml2_cfa_kappa(ml2_state_t *st, int kappa)
{
    dpe_t tmp;
    ibz_t gram;
    const int z = st->zeta;
    for (int j = z; j < kappa; j++) {
        if (!ml2_gram_entry(&gram, st, kappa, j))
            goto precision;
        dpe_set_z(st->r_row[j], &gram);
        for (int k = z; k < j; k++) {
            dpe_mul(tmp, st->r_row[k], ML2_MU(st, j, k));
            dpe_sub(st->r_row[j], st->r_row[j], tmp);
        }
        if (dpe_cmp_d(st->r_diag[j], 0.0) <= 0)
            return 0;
        dpe_div(ML2_MU(st, kappa, j), st->r_row[j], st->r_diag[j]);
    }

    if (!ml2_gram_entry(&gram, st, kappa, kappa))
        goto precision;
    dpe_set_z(st->s[z], &gram);
    for (int j = z; j < kappa; j++) {
        dpe_mul(tmp, ML2_MU(st, kappa, j), st->r_row[j]);
        dpe_sub(st->s[j + 1], st->s[j], tmp);
    }
    dpe_set(st->r_diag[kappa], st->s[kappa]);
    return 1;

precision:
    st->failure = QUAT_ML2_ERR_PRECISION;
    return 0;
}

/* For the old basis and the complete Babai vector X, form the exact
 * b'[k] = b[k] - sum_i X[i] b[i]. No CFA observes a partially updated row.
 * In the next CFA, bilinearity makes the recomputed Gram entry
 * <b'[k], b[h]> equal to G[k,h] - sum_i X[i] G[i,h] for h != k.
 * The same expansion of <b'[k], b'[k]> gives the snapshot diagonal. */
static int
ml2_apply_size_reduction(ml2_state_t *st, int kappa)
{
    ibz_vec_4_t candidate;
    ibz_t correction, term;
    const int z = st->zeta;
    for (int c = 0; c < 4; c++) {
        ibz_set(&correction, 0);
        for (int i = z; i < kappa; i++) {
            if (ibz_is_zero(&st->coefficients[i]))
                continue;
            if (!ibz_mul_checked(&term, &st->coefficients[i], &st->b[i][c]) ||
                !ibz_add_checked(&correction, &correction, &term))
                return 0;
        }
        if (!ibz_sub_checked(&candidate[c], &st->b[kappa][c], &correction))
            return 0;
    }
    ibz_vec_4_copy(&st->b[kappa], &candidate);
    return 1;
}

static void
ml2_lazy_size_reduce(ml2_state_t *st, int kappa)
{
    dpe_t tmp, rounded;
    const int z = st->zeta;
    for (int iter = 0; iter < 65536; iter++) {
        if (!ml2_cfa_kappa(st, kappa)) {
            if (!st->failure)
                st->failure = -1;
            return;
        }
        int need_reduce = 0;
        for (int j = z; j < kappa; j++) {
            dpe_abs(tmp, ML2_MU(st, kappa, j));
            if (dpe_cmp_d(tmp, ETABAR) > 0) {
                need_reduce = 1;
                break;
            }
        }
        if (!need_reduce)
            return;

        for (int i = kappa - 1; i >= z; i--) {
            dpe_round(rounded, ML2_MU(st, kappa, i));
            if (!dpe_get_z(&st->coefficients[i], rounded)) {
                st->failure = QUAT_ML2_ERR_PRECISION;
                return;
            }
            for (int j = z; j < i; j++) {
                dpe_mul(tmp, rounded, ML2_MU(st, i, j));
                dpe_sub(ML2_MU(st, kappa, j), ML2_MU(st, kappa, j), tmp);
            }
        }
        if (!ml2_apply_size_reduction(st, kappa)) {
            st->failure = QUAT_ML2_ERR_PRECISION;
            return;
        }
    }
    st->failure = -1;
}

static void
ml2_insert_before(ml2_state_t *st, int from, int to)
{
    for (int i = from; i > to; i--)
        for (int c = 0; c < 4; c++)
            ibz_swap(&st->b[i][c], &st->b[i - 1][c]);
}

/* ========== ML2 main loop (paper Algorithm 13, 0-indexed) ==========
 * Init after exact-zero partitioning:
 *   r[zeta][zeta] = G[zeta][zeta]; kappa = zeta + 1.
 * Loop: while kappa < d
 *   LazySizeReduce(kappa)
 *   kappa_prime = kappa
 *   while kappa >= zeta + 1 && delta_bar * r[kappa-1][kappa-1] >= s[kappa-1] : kappa--
 *   copy mu, r row from kappa_prime into kappa, r[kappa][kappa] = s[kappa]
 *   insert b[kappa_prime] before b[kappa]
 *   if b[kappa] is zero: zeta++
 *   kappa = max(zeta + 1, kappa + 1)
 */
static void
ml2_main_loop(ml2_state_t *st)
{
    dpe_t delta_r;
    dpe_init(delta_r);

    if (st->zeta == st->d) {
        dpe_clear(delta_r);
        return;
    }

    if (!ml2_cfa_kappa(st, st->zeta)) {
        if (!st->failure)
            st->failure = -1;
        return;
    }
    int kappa = st->zeta + 1;
    int outer_iter = 0;
    /* paper Issue 13 (author reply): "size-reduce keeps entry count unchanged and size monotonically non-increasing".
     * Paper ML2 (NS09 Fig 9) does not have abort/oscillation detection — convergence
     * guaranteed by Lemma lazy-size-reduction-bound + size-monotonic invariant.
     * Set OUTER_MAX large enough that termination is dominated by paper's convergence
     * theorem, not by our cap. No same-state early-abort. */
    const int OUTER_MAX = 65536;

    while (kappa < st->d) {
        outer_iter++;
        if (outer_iter > OUTER_MAX) {
            st->failure = -1;
            break;
        }

        ml2_lazy_size_reduce(st, kappa);
        if (st->failure)
            break;

        int kappa_prime = kappa;
        while (kappa >= st->zeta + 1) {
            dpe_mul(delta_r, st->delta_bar, st->r_diag[kappa - 1]);
            if (dpe_cmp(delta_r, st->s[kappa - 1]) < 0)
                break;
            kappa--;
        }

        if (kappa != kappa_prime) {
            for (int i = st->zeta; i < kappa; i++)
                dpe_set(ML2_MU(st, kappa, i),
                        ML2_MU(st, kappa_prime, i));
            dpe_set(st->r_diag[kappa], st->s[kappa]);
            ml2_insert_before(st, kappa_prime, kappa);
        }

        if (ibz_vec_4_is_zero(&st->b[kappa])) {
            /* paper invariant: zero vectors accumulate at positions 0..zeta-1.
             * Output is b[zeta..d-1]. Swap the just-detected zero into b[zeta]
             * so CFA's j range (zeta+1..kappa-1) skips it. */
            if (kappa != st->zeta) {
                for (int c = 0; c < 4; c++)
                    ibz_swap(&st->b[kappa][c], &st->b[st->zeta][c]);
            }
            st->zeta++;
            /* Insertion of a dependent generator can expose a new first
             * active vector.  Restart its factorization from an exact,
             * strictly positive diagonal before processing another row;
             * otherwise the stale zero diagonal again becomes a divisor. */
            if (st->zeta < st->d && !ml2_cfa_kappa(st, st->zeta)) {
                if (!st->failure)
                    st->failure = -1;
                break;
            }
            kappa = st->zeta + 1;
            continue;
        }

        /* Insertion invalidates the later Gram--Schmidt rows.  Recompute
         * each lazily in LazySizeReduce, after its predecessors have been
         * reduced and any zero vectors removed.  Eagerly rebuilding all
         * rows here can divide by a zero diagonal from a dependent row. */

        int next = kappa + 1;
        if (st->zeta + 1 > next)
            next = st->zeta + 1;
        kappa = next;
    }

    dpe_clear(delta_r);
}

/* If nonzero positive integers x and y satisfy x < 2^x_bits and
 * y < 2^y_bits, then x+y < 2^max(x_bits,y_bits)+1.  Zero has bound zero. */
static int
ml2_positive_sum_bound(int x_bits, int y_bits)
{
    if (x_bits == 0)
        return y_bits;
    if (y_bits == 0)
        return x_bits;
    return (x_bits > y_bits ? x_bits : y_bits) + 1;
}

/* Section 4 reserves 512 times the largest initial squared norm for
 * exact ML2 intermediates. Bound each weighted coordinate before forming the
 * initial Gram matrix. Near the boundary use the exact diagonal alone. */
static int
ml2_input_fits_precision(const ibz_vec_4_t *input,
                         int d,
                         const quat_alg_t *alg)
{
    int p_bits = alg == NULL ? 0 : ibz_bitsize(&alg->p);

    for (int i = 0; i < d; i++) {
        int norm_bound = 0;
        for (int coordinate = 0; coordinate < 4; coordinate++) {
            int bits = ibz_bitsize(&input[i][coordinate]);
            int term_bound = bits == 0 ? 0 : 2 * bits;
            if (term_bound != 0 && coordinate >= 2 && alg != NULL &&
                !ibz_is_one(&alg->p))
                term_bound += p_bits;
            norm_bound = ml2_positive_sum_bound(norm_bound, term_bound);
        }
        if (norm_bound <= IBZ_BITS - 1 - 9)
            continue;

        ibz_t norm, term;
        ibz_set(&norm, 0);
        for (int c = 0; c < 4; c++) {
            const ibz_t *value = &input[i][c];
            if (!ibz_mul_checked(&term, value, value) ||
                (c >= 2 && alg != NULL && !ibz_mul_checked(&term, &term, &alg->p)) ||
                !ibz_add_checked(&norm, &norm, &term) ||
                ibz_bitsize(&norm) > IBZ_BITS - 1 - 9)
                return 0;
        }
    }
    return 1;
}

static int ml2_retry_source_index(int attempt, int position, int d);

typedef struct {
    const ibz_mat_4x4_t *left;
    const ibz_mat_4x4_t *right;
} ml2_product_bases_t;

static int
ml2_load_product(ibz_vec_4_t *out, const ml2_product_bases_t *bases, int index, const quat_alg_t *alg)
{
    ibz_vec_4_t left;
    for (int c = 0; c < 4; c++) {
        ibz_copy(&left[c], &(*bases->left)[c][index / 4]);
        ibz_copy(&(*out)[c], &(*bases->right)[c][index % 4]);
    }
    if (!quat_alg_coord_mul_fits(&left, out, alg))
        return 0;
    quat_alg_coord_mul(out, &left, out, alg);
    return 1;
}

/* ========== entry point ==========
 * Reduce d input generators. On return, the LLL-reduced basis is in
 * st.b[zeta..d-1] (rank = d - zeta). Copies up to `out_capacity` of those
 * into `output`. `publish_lower_rank` lets the full-rank retry path leave its
 * caller's output untouched after an unsuccessful lower-rank attempt.
 */
static int
ml2_reduce_ordered(ibz_vec_4_t *output,
                   int out_capacity,
                   const ibz_vec_4_t *input,
                   int d,
                   const quat_alg_t *alg,
                   int order_attempt,
                   int publish_lower_rank,
                   const ml2_product_bases_t *products)
{
    if ((input == NULL && products == NULL) || (out_capacity > 0 && output == NULL) ||
        d < 1 || d > ML2_MAX_D || out_capacity < 0 ||
        order_attempt < 0 ||
        order_attempt >= QUAT_ML2_RETRY_MAX_ATTEMPTS)
        return -1;
    if (alg != NULL && ibz_cmp(&alg->p, &ibz_const_zero) <= 0)
        return -1;

    if (input != NULL && !ml2_input_fits_precision(input, d, alg))
        return QUAT_ML2_ERR_PRECISION;

#if defined(TARGET_ARM)
    ml2_state_t *st = sqisign_m4_workspace_acquire(sizeof(*st) + d * sizeof(ibz_vec_4_t));
    if (st == NULL)
        return -1;
#else
    ml2_state_t automatic_state;
    ml2_state_t *st = &automatic_state;
#endif
    ml2_state_init(st, d, alg);
    int rho = -1;
    for (int i = 0; i < d; i++) {
        int source = ml2_retry_source_index(order_attempt, i, d);
        if (products == NULL) {
            ibz_vec_4_copy(&st->b[i], &input[source]);
        } else if (!ml2_load_product(&st->b[i], products, source, alg)) {
            rho = QUAT_ML2_ERR_PRECISION;
            goto cleanup;
        }
    }
    if (products != NULL && !ml2_input_fits_precision(st->b, d, alg)) {
        rho = QUAT_ML2_ERR_PRECISION;
        goto cleanup;
    }
    ml2_partition_initial_zeros(st);

    const int rounding = fegetround();
    if (rounding < 0 || (rounding != FE_TONEAREST && fesetround(FE_TONEAREST)))
        goto cleanup;
    ml2_main_loop(st);
    if (rounding != FE_TONEAREST && fesetround(rounding))
        goto cleanup;
    if (st->failure) {
        rho = st->failure;
        goto cleanup;
    }
    rho = st->d - st->zeta;
    if (rho > 4) {
        rho = -1;
        goto cleanup;
    }
    if (rho == 4 || publish_lower_rank) {
        const int n = rho < out_capacity ? rho : out_capacity;
        for (int i = 0; i < n; i++)
            ibz_vec_4_copy(&output[i], &st->b[st->zeta + i]);
    }

cleanup:
    ml2_state_finalize(st);
#if defined(TARGET_ARM)
    sqisign_m4_workspace_release(st);
#endif
    return rho;
}

int
quat_ml2(ibz_vec_4_t *output,
         int out_capacity,
         const ibz_vec_4_t *input,
         int d,
         const quat_alg_t *alg)
{
    return ml2_reduce_ordered(
        output, out_capacity, input, d, alg, 0, 1, NULL);
}

/* Return the source position for one of the retry permutations. Generator
 * coordinates are copied verbatim: there are no sign changes or arithmetic
 * operations, so retries cannot increase the fixed-precision bit budget. */
static int
ml2_retry_source_index(int attempt, int position, int d)
{
    switch (attempt) {
    case 1: /* one-step cyclic rotation */
        return (position + 1) % d;
    case 2: /* reverse order */
        return d - 1 - position;
    case 3: { /* even positions, followed by odd positions */
        int even_count = (d + 1) / 2;
        if (position < even_count)
            return 2 * position;
        return 2 * (position - even_count) + 1;
    }
    default:
        return position;
    }
}

#if defined(SQISIGN_ML2_PROFILE)
typedef struct {
    atomic_uint_fast64_t inputs;
    atomic_uint_fast64_t precision_rejected;
    atomic_uint_fast64_t first_attempt_failures;
    atomic_uint_fast64_t recovered[QUAT_ML2_RETRY_MAX_ATTEMPTS - 1];
    atomic_uint_fast64_t exhausted;
    atomic_uint_fast64_t underlying_attempts;
} ml2_atomic_profile_dimension_t;

static struct {
    ml2_atomic_profile_dimension_t d4;
    ml2_atomic_profile_dimension_t d8;
    ml2_atomic_profile_dimension_t d16;
} ml2_profile;

static ml2_atomic_profile_dimension_t *
ml2_profile_for_dimension(int d)
{
    if (d == 4)
        return &ml2_profile.d4;
    if (d == 8)
        return &ml2_profile.d8;
    if (d == 16)
        return &ml2_profile.d16;
    return NULL;
}

static void
ml2_profile_dimension_reset(ml2_atomic_profile_dimension_t *dimension)
{
    atomic_store(&dimension->inputs, 0);
    atomic_store(&dimension->precision_rejected, 0);
    atomic_store(&dimension->first_attempt_failures, 0);
    for (int i = 0; i < QUAT_ML2_RETRY_MAX_ATTEMPTS - 1; i++)
        atomic_store(&dimension->recovered[i], 0);
    atomic_store(&dimension->exhausted, 0);
    atomic_store(&dimension->underlying_attempts, 0);
}

static void
ml2_profile_dimension_get(quat_ml2_profile_dimension_t *output,
                          const ml2_atomic_profile_dimension_t *dimension)
{
    output->inputs = atomic_load(&dimension->inputs);
    output->precision_rejected = atomic_load(&dimension->precision_rejected);
    output->first_attempt_failures = atomic_load(&dimension->first_attempt_failures);
    for (int i = 0; i < QUAT_ML2_RETRY_MAX_ATTEMPTS - 1; i++)
        output->recovered[i] = atomic_load(&dimension->recovered[i]);
    output->exhausted = atomic_load(&dimension->exhausted);
    output->underlying_attempts = atomic_load(&dimension->underlying_attempts);
}

void
quat_ml2_profile_reset(void)
{
    ml2_profile_dimension_reset(&ml2_profile.d4);
    ml2_profile_dimension_reset(&ml2_profile.d8);
    ml2_profile_dimension_reset(&ml2_profile.d16);
}

void
quat_ml2_profile_get(quat_ml2_profile_t *output)
{
    if (output == NULL)
        return;
    ml2_profile_dimension_get(&output->d4, &ml2_profile.d4);
    ml2_profile_dimension_get(&output->d8, &ml2_profile.d8);
    ml2_profile_dimension_get(&output->d16, &ml2_profile.d16);
}
#endif

/* Native retry driver for ML2.  The requested permutation is folded into
 * ML2's initial exact input copy, so retries need neither a 16-generator
 * permutation buffer nor a second four-vector publication buffer.  ML2
 * publishes rank four directly; a lower-rank first result is published only
 * for generic MLLL. */
static int
ml2_retry_driver_native(ibz_vec_4_t *output,
                        int out_capacity,
                        const ibz_vec_4_t *input,
                        int d,
                        const quat_alg_t *alg,
                        int retain_first_lower_rank,
                        const ml2_product_bases_t *products)
{
    int first_rho = -1;
    int result = -1;
#if defined(SQISIGN_ML2_PROFILE)
    ml2_atomic_profile_dimension_t *profile_dimension;
#endif

    if (output == NULL || (input == NULL && products == NULL) ||
        d < 1 || d > ML2_MAX_D || out_capacity < 4)
        return -1;
    if (alg != NULL && ibz_cmp(&alg->p, &ibz_const_zero) <= 0)
        return -1;

#if defined(SQISIGN_ML2_PROFILE)
    profile_dimension = ml2_profile_for_dimension(d);
    if (profile_dimension != NULL)
        atomic_fetch_add(&profile_dimension->inputs, 1);
#endif

    if (input != NULL && !ml2_input_fits_precision(input, d, alg)) {
#if defined(SQISIGN_ML2_PROFILE)
        if (profile_dimension != NULL)
            atomic_fetch_add(&profile_dimension->precision_rejected, 1);
#endif
        return QUAT_ML2_ERR_PRECISION;
    }

    for (int attempt = 0; attempt < QUAT_ML2_RETRY_MAX_ATTEMPTS; attempt++) {
        int rho;
#if defined(SQISIGN_ML2_PROFILE)
        if (profile_dimension != NULL)
            atomic_fetch_add(&profile_dimension->underlying_attempts, 1);
#endif
        rho = ml2_reduce_ordered(output,
                                 out_capacity,
                                 input,
                                 d,
                                 alg,
                                 attempt,
                                 retain_first_lower_rank && attempt == 0,
                                 products);
        if (attempt == 0)
            first_rho = rho;

        if (rho == 4) {
            result = 4;
#if defined(SQISIGN_ML2_PROFILE)
            if (profile_dimension != NULL && attempt > 0)
                atomic_fetch_add(&profile_dimension->recovered[attempt - 1], 1);
#endif
            break;
        }
#if defined(SQISIGN_ML2_PROFILE)
        if (profile_dimension != NULL && attempt == 0)
            atomic_fetch_add(&profile_dimension->first_attempt_failures, 1);
#endif
        if (attempt == 0 && retain_first_lower_rank && rho >= 0 && rho <= 4) {
            result = rho;
            break;
        }
    }

    if (result < 0 || (!retain_first_lower_rank && result != 4)) {
        result = first_rho;
#if defined(SQISIGN_ML2_PROFILE)
        if (profile_dimension != NULL)
            atomic_fetch_add(&profile_dimension->exhausted, 1);
#endif
    }

    return result;
}

int
quat_ml2_retry_products(ibz_mat_4x4_t *out,
                        const ibz_mat_4x4_t *left,
                        const ibz_mat_4x4_t *right,
                        const quat_alg_t *alg)
{
    if (out == NULL || left == NULL || right == NULL || alg == NULL)
        return -1;
    const ml2_product_bases_t products = { left, right };
    return ml2_retry_driver_native(*out, 4, NULL, 16, alg, 0, &products);
}

/* Callback retry driver retained for fault-injection tests and alternate
 * reducers.  `retain_first_lower_rank` is used by generic MLLL:
 * a nonnegative rank below four is a valid result and is published without
 * pointless permutations.  Full-rank lattice operations leave it false and
 * publish only a rank-four attempt.  Keeping attempt 0 inside this driver is
 * also essential for profiling: one logical input and every actually
 * executed reduction attempt are accounted for in exactly one place. */
#if !defined(TARGET_ARM)
static int
ml2_retry_driver_callback(ibz_vec_4_t *output,
                          int out_capacity,
                          const ibz_vec_4_t *input,
                          int d,
                          const quat_alg_t *alg,
                          quat_ml2_reducer_t reducer,
                          int retain_first_lower_rank)
{
    ibz_vec_4_t permuted[ML2_MAX_D];
    ibz_vec_4_t scratch[4];
    int first_rho = -1;
    int result = -1;
#if defined(SQISIGN_ML2_PROFILE)
    ml2_atomic_profile_dimension_t *profile_dimension;
#endif

    if (output == NULL || input == NULL || reducer == NULL ||
        d < 1 || d > ML2_MAX_D || out_capacity < 4)
        return -1;
    if (alg != NULL && ibz_cmp(&alg->p, &ibz_const_zero) <= 0)
        return -1;

#if defined(SQISIGN_ML2_PROFILE)
    profile_dimension = ml2_profile_for_dimension(d);
    if (profile_dimension != NULL)
        atomic_fetch_add(&profile_dimension->inputs, 1);
#endif

    if (!ml2_input_fits_precision(input, d, alg)) {
#if defined(SQISIGN_ML2_PROFILE)
        if (profile_dimension != NULL)
            atomic_fetch_add(&profile_dimension->precision_rejected, 1);
#endif
        return QUAT_ML2_ERR_PRECISION;
    }

    for (int i = 0; i < d; i++)
        ibz_vec_4_init(&permuted[i]);
    for (int i = 0; i < 4; i++)
        ibz_vec_4_init(&scratch[i]);

    for (int attempt = 0; attempt < QUAT_ML2_RETRY_MAX_ATTEMPTS; attempt++) {
        const ibz_vec_4_t *attempt_input = input;

        if (attempt != 0) {
            for (int i = 0; i < d; i++) {
                int source = ml2_retry_source_index(attempt, i, d);
                for (int c = 0; c < 4; c++)
                    ibz_copy(&permuted[i][c], &input[source][c]);
            }
            attempt_input = permuted;
        }

        int rho;
#if defined(SQISIGN_ML2_PROFILE)
        if (profile_dimension != NULL)
            atomic_fetch_add(&profile_dimension->underlying_attempts, 1);
#endif
        rho = reducer(scratch, 4, attempt_input, d, alg);
        if (attempt == 0)
            first_rho = rho;

        if (rho == 4) {
            for (int i = 0; i < 4; i++)
                for (int c = 0; c < 4; c++)
                    ibz_copy(&output[i][c], &scratch[i][c]);
            result = 4;
#if defined(SQISIGN_ML2_PROFILE)
            if (profile_dimension != NULL && attempt > 0)
                atomic_fetch_add(&profile_dimension->recovered[attempt - 1], 1);
#endif
            break;
        }
#if defined(SQISIGN_ML2_PROFILE)
        if (profile_dimension != NULL && attempt == 0)
            atomic_fetch_add(&profile_dimension->first_attempt_failures, 1);
#endif
        if (attempt == 0 && retain_first_lower_rank && rho >= 0 && rho <= 4) {
            for (int i = 0; i < rho; i++)
                for (int c = 0; c < 4; c++)
                    ibz_copy(&output[i][c], &scratch[i][c]);
            result = rho;
            break;
        }
    }

    if (result < 0 || (!retain_first_lower_rank && result != 4)) {
        result = first_rho;
#if defined(SQISIGN_ML2_PROFILE)
        if (profile_dimension != NULL)
            atomic_fetch_add(&profile_dimension->exhausted, 1);
#endif
    }

    for (int i = 0; i < 4; i++)
        ibz_vec_4_finalize(&scratch[i]);
    for (int i = 0; i < d; i++)
        ibz_vec_4_finalize(&permuted[i]);
    return result;
}
#endif

static int
ml2_retry_driver(ibz_vec_4_t *output,
                 int out_capacity,
                 const ibz_vec_4_t *input,
                 int d,
                 const quat_alg_t *alg,
                 quat_ml2_reducer_t reducer,
                 int retain_first_lower_rank)
{
    if (reducer == quat_ml2)
        return ml2_retry_driver_native(output,
                                       out_capacity,
                                       input,
                                       d,
                                       alg,
                                       retain_first_lower_rank,
                                       NULL);
#if defined(TARGET_ARM)
    /* Alternate reducer callbacks are test instrumentation.  Excluding their
     * permutation workspace from the embedded build keeps an unreachable
     * 16-vector frame out of both the binary and conservative stack audits. */
    return -1;
#else
    return ml2_retry_driver_callback(output,
                                     out_capacity,
                                     input,
                                     d,
                                     alg,
                                     reducer,
                                     retain_first_lower_rank);
#endif
}

int
quat_ml2_retry_with_reducer(ibz_vec_4_t *output,
                            int out_capacity,
                            const ibz_vec_4_t *input,
                            int d,
                            const quat_alg_t *alg,
                            quat_ml2_reducer_t reducer)
{
    return ml2_retry_driver(
        output, out_capacity, input, d, alg, reducer, 0);
}

int
quat_ml2_mlll_with_reducer(ibz_vec_4_t *output,
                            int out_capacity,
                            const ibz_vec_4_t *input,
                            int d,
                            const quat_alg_t *alg,
                            quat_ml2_reducer_t reducer)
{
    return ml2_retry_driver(
        output, out_capacity, input, d, alg, reducer, 1);
}

int
quat_ml2_retry(ibz_vec_4_t *output,
               int out_capacity,
               const ibz_vec_4_t *input,
               int d,
               const quat_alg_t *alg)
{
    return ml2_retry_driver_native(
        output, out_capacity, input, d, alg, 0, NULL);
}

int
quat_ml2_mlll(ibz_vec_4_t *output,
              int out_capacity,
              const ibz_vec_4_t *input,
              int d,
              const quat_alg_t *alg)
{
    return ml2_retry_driver_native(
        output, out_capacity, input, d, alg, 1, NULL);
}
