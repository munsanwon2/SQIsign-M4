/* ML2 — Modified L2 algorithm (NS09 Figure 9), generic d-rank input.
 *
 * Implements:
 *   - Algorithm \ref{alg:ML2}                (03Ideal.tex)
 *   - Algorithm \ref{alg:lazy-size-reduction} (91proof.tex)
 *
 * Model (paper Issue 13):
 *   - exact ibz_t for vectors b_i and on-demand integer Gram entries
 *   - 53-bit float (dpe_t) for Cholesky r, GSO mu, intermediate s
 *
 * Indexing: paper is 1-indexed, this file is 0-indexed.
 *   - paper b_1..b_d        ->  b[0..d-1]
 *   - paper r,mu             ->  r[i][j], mu[i][j] for j <= i, both 0-indexed
 *   - paper zeta in [0..d-1] ->  zeta in [0..d-1]  (same semantics: skip count)
 *   - paper kappa in [1..d]  ->  kappa in [0..d-1]
 *
 * delta_bar/eta_bar use repo defaults DELTABAR=0.995, ETABAR=0.505 (Schnorr-Euchner).
 * Paper uses 0.875/0.525; the repo defaults LLL-reduce a strict superset of bases
 * the paper does and converge faster.
 */

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

#ifndef ML2_MAX_D
/* d=16 needed for paper Issue 11/12 spec-faithful quat_lideal_lideal_mul_reduced:
 * the bar(J_t)·I product yields 16 column generators that we LLL directly
 * (avoids hnf_mod_core's m^2 transient blow-up at 28-limb). */
#define ML2_MAX_D 16
#endif

#define ML2_MU_COUNT ((ML2_MAX_D * (ML2_MAX_D - 1)) / 2)
#define ML2_MU(ST, I, J) \
    ((ST)->mu[((I) * ((I) - 1)) / 2 + (J)])

/* ========== state ========== */

typedef struct {
    int d;
    int zeta;
    int aborted; /* 1 if main_loop hit oscillation/iter cap */
    /* NULL selects the Euclidean fallback used by algebra-independent
     * span-only callers.  Compact quaternion algorithms always pass alg and
     * therefore use the reduced-norm form diag(1, 1, p, p). */
    const ibz_t *p;
    ibz_vec_4_t b[ML2_MAX_D];
    /* Only the current CFA row and the persistent diagonal of r are ever
     * read.  Historical strict-lower mu rows are packed densely. */
    dpe_t r_diag[ML2_MAX_D];
    dpe_t r_row[ML2_MAX_D];
    dpe_t mu[ML2_MU_COUNT];
    dpe_t s[ML2_MAX_D + 1];
    dpe_t delta_bar;
} ml2_state_t;

#if defined(TARGET_ARM)
/* Reserve room in the 60-KiB shared CCM arena for allocator metadata and
 * later phase overlays.  The largest current state is 29,360 bytes at lvl-V. */
_Static_assert(sizeof(ml2_state_t) <= 29360u,
               "ML2 state exceeds the Cortex-M4 CCM scratch budget");
#endif

static void
ml2_state_init(ml2_state_t *st, int d, const quat_alg_t *alg)
{
    st->d = d;
    st->zeta = 0;
    st->aborted = 0;
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

/* Compute one exact integer Gram entry on demand.  In the quaternion
 * algebra (-1,-p), half the trace pairing (equivalently, the reduced-norm
 * bilinear form) is
 *
 *   <a,b> = a0*b0 + a1*b1 + p*(a2*b2 + a3*b3).
 *
 * The omitted common factor 2 has no effect on ML2's decisions.  A NULL p is
 * the documented Euclidean fallback for algebra-independent span reduction.
 *
 * Keeping the full d-by-d fixed-width Gram matrix would dominate transient
 * memory (more than 100 KiB at level V for d=16).  Every use of a Gram entry
 * is followed by conversion to dpe_t, and only b[kappa] changes during lazy
 * size reduction, so recomputing the requested entry from the current exact
 * basis is equivalent to maintaining that matrix.  `scratch` must not alias
 * `value`. */
static void
ml2_gram_entry(ibz_t *value,
               ibz_t *scratch,
               const ml2_state_t *st,
               int i,
               int j)
{
    ibz_set(value, 0);
    for (int k = 0; k < 4; k++) {
        ibz_mul(scratch, &st->b[i][k], &st->b[j][k]);
        if (k >= 2 && st->p != NULL)
            ibz_mul(scratch, scratch, st->p);
        ibz_add(value, value, scratch);
    }
}

/* ========== CFA (steps 2-7 of NS09 Figure 2) ==========
 * Compute r[kappa][j], mu[kappa][j] for zeta < j < kappa, and
 * s[zeta..kappa] (partial sums of <b_kappa, b_kappa> minus mu^2*r diagonals).
 *
 * r[i][j] = <b_i, b_j*> ; mu[i][j] = r[i][j] / r[j][j]
 * s[k] = G[kappa][kappa] - sum_{j=zeta..k-1} mu[kappa][j]^2 * r[j][j]
 */
static void
ml2_cfa_kappa(ml2_state_t *st, int kappa)
{
    dpe_t tmp;
    ibz_t gram, gram_scratch;
    dpe_init(tmp);
    ibz_init(&gram);
    ibz_init(&gram_scratch);

    int z = st->zeta;
    for (int j = z; j < kappa; j++) {
        ml2_gram_entry(&gram, &gram_scratch, st, kappa, j);
        dpe_set_z(st->r_row[j], &gram);
        for (int k = z; k < j; k++) {
            dpe_mul(tmp, st->r_row[k], ML2_MU(st, j, k));
            dpe_sub(st->r_row[j], st->r_row[j], tmp);
        }
        dpe_div(ML2_MU(st, kappa, j), st->r_row[j], st->r_diag[j]);
    }

    /* s[z] = G[kappa][kappa]; s[k+1] = s[k] - mu[kappa][k] * r[kappa][k] */
    ml2_gram_entry(&gram, &gram_scratch, st, kappa, kappa);
    dpe_set_z(st->s[z], &gram);
    for (int k = z; k < kappa; k++) {
        dpe_mul(tmp, ML2_MU(st, kappa, k), st->r_row[k]);
        dpe_sub(st->s[k + 1], st->s[k], tmp);
    }
    /* r[kappa][kappa] := s[kappa] */
    dpe_set(st->r_diag[kappa], st->s[kappa]);

    ibz_finalize(&gram_scratch);
    ibz_finalize(&gram);
    dpe_clear(tmp);
}

/* ========== LazySizeReduce (paper Algorithm 14) ==========
 * Repeatedly: CFA -> if all |mu[kappa][j]| <= eta_bar, return;
 * else for i = kappa-1 down to zeta: X_i = round(mu[kappa][i]); update mu;
 * then b[kappa] -= sum X_i * b[i]; loop.  The next CFA recomputes the
 * affected exact Gram row directly from the updated basis.
 */
static void
ml2_lazy_size_reduce(ml2_state_t *st, int kappa)
{
    dpe_t tmp_f, X_f;
    ibz_t X_z, tmp_z;
    dpe_init(tmp_f);
    dpe_init(X_f);
    ibz_init(&X_z);
    ibz_init(&tmp_z);

    int z = st->zeta;

    /* paper Algorithm LazySizeReduce: loop until max |mu_kappa,j| <= eta_bar.
     * Paper guarantees size-monotonic convergence (Issue 13 author reply).
     * Large cap to let paper's monotonic argument dominate, not artificial abort. */
    for (int iter = 0; iter < 65536; iter++) {
        ml2_cfa_kappa(st, kappa);

        /* check if all |mu[kappa][j]| <= eta_bar */
        int need_reduce = 0;
        for (int j = z; j < kappa; j++) {
            if (dpe_cmp_d(ML2_MU(st, kappa, j), ETABAR) > 0 ||
                dpe_cmp_d(ML2_MU(st, kappa, j), -ETABAR) < 0) {
                need_reduce = 1;
                break;
            }
        }
        if (!need_reduce)
            goto done;

        /* for i = kappa-1 down to z: X_i = round(mu[kappa][i]); update mu[kappa][j]
         * for j < i; then b[kappa] -= X_i * b[i]. */
        for (int i = kappa - 1; i >= z; i--) {
            dpe_set(X_f, ML2_MU(st, kappa, i));
            dpe_round(X_f, X_f);
            dpe_get_z(&X_z, X_f);
            if (ibz_is_zero(&X_z))
                continue;

            /* mu[kappa][j] -= X_i * mu[i][j] for j = z..i-1 (Fp arithmetic) */
            for (int j = z; j < i; j++) {
                dpe_mul(tmp_f, X_f, ML2_MU(st, i, j));
                dpe_sub(ML2_MU(st, kappa, j),
                        ML2_MU(st, kappa, j),
                        tmp_f);
            }

            /* b[kappa] -= X_i * b[i]  (exact integer vector update) */
            for (int c = 0; c < 4; c++) {
                ibz_mul(&tmp_z, &X_z, &st->b[i][c]);
                ibz_sub(&st->b[kappa][c], &st->b[kappa][c], &tmp_z);
            }

        }
    }

    /* The paper algorithm has no cap, but a fixed-precision implementation
     * must not publish a partially size-reduced basis if its safeguard is
     * reached.  Mark the attempt as failed so the wrapper can retry ML2 with
     * another span-preserving generator order. */
    st->aborted = 1;

done:
    dpe_clear(tmp_f);
    dpe_clear(X_f);
    ibz_finalize(&X_z);
    ibz_finalize(&tmp_z);
}

/* Insert b[from] right before b[to], 0 <= to <= from < d.
 * After: b'[to] = b[from], b'[to+1..from] = old b[to..from-1].
 * Exact Gram entries are computed on demand, so no matrix permutation or
 * rebuild is required here.
 */
static void
ml2_insert_before(ml2_state_t *st, int from, int to)
{
    if (from == to)
        return;
    /* rotate b */
    ibz_vec_4_t tmp;
    ibz_vec_4_init(&tmp);
    for (int c = 0; c < 4; c++)
        ibz_copy(&tmp[c], &st->b[from][c]);
    for (int i = from; i > to; i--)
        for (int c = 0; c < 4; c++)
            ibz_copy(&st->b[i][c], &st->b[i - 1][c]);
    for (int c = 0; c < 4; c++)
        ibz_copy(&st->b[to][c], &tmp[c]);
    ibz_vec_4_finalize(&tmp);
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

    ml2_cfa_kappa(st, st->zeta);
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
            st->aborted = 1;
            break;
        }

        ml2_lazy_size_reduce(st, kappa);
        if (st->aborted)
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
            if (st->zeta < st->d)
                ml2_cfa_kappa(st, st->zeta);
            kappa = st->zeta + 1;
            continue;
        }

        if (kappa != kappa_prime) {
            /* Insertion reorders b[zeta..kappa_prime], so r[j][*] and
             * mu[j][*] for j in (kappa,kappa_prime] are stale.  Recompute
             * those rows only after ruling out a zero at the insertion
             * point, since a zero row has no valid Gram--Schmidt divisor. */
            for (int j = kappa + 1; j <= kappa_prime; j++)
                ml2_cfa_kappa(st, j);
        }

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

/* Lemma 8 bounds every exact ML2 integer by max_i ||b_i||_p^2.  Bound each
 * initial diagonal Gram entry coordinate-by-coordinate, so a large p is
 * charged only to coordinates 2 and 3:
 *
 *   b0^2 + b1^2 + p*(b2^2 + b3^2).
 *
 * This is both rigorous before any fixed-width multiplication takes place
 * and substantially tighter than adding bits(p) to twice the largest of all
 * four coordinates.  Permuting generators cannot change the bound. */
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
        if (norm_bound > IBZ_BITS - 1)
            return 0;
    }
    return 1;
}

static int ml2_retry_source_index(int attempt, int position, int d);

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
                   int publish_lower_rank)
{
    if (input == NULL || (out_capacity > 0 && output == NULL) ||
        d < 1 || d > ML2_MAX_D || out_capacity < 0 ||
        order_attempt < 0 ||
        order_attempt >= QUAT_ML2_RETRY_MAX_ATTEMPTS)
        return -1;
    if (alg != NULL && ibz_cmp(&alg->p, &ibz_const_zero) <= 0)
        return -1;

    if (!ml2_input_fits_precision(input, d, alg))
        return QUAT_ML2_ERR_PRECISION;

#if defined(TARGET_ARM)
    ml2_state_t *st =
        (ml2_state_t *)sqisign_m4_workspace_acquire(sizeof(ml2_state_t));
    if (st == NULL)
        return -1;
#else
    ml2_state_t automatic_state;
    ml2_state_t *st = &automatic_state;
#endif
    ml2_state_init(st, d, alg);
    for (int i = 0; i < d; i++) {
        int source = ml2_retry_source_index(order_attempt, i, d);
        for (int c = 0; c < 4; c++) {
            ibz_copy(&st->b[i][c], &input[source][c]);
        }
    }
    ml2_partition_initial_zeros(st);

    ml2_main_loop(st);

    int rho;
    int n;
    if (st->aborted) {
        /* Signal to caller: oscillation/iter-cap. Return -1 so caller can
         * fallback to HNF (or another rank-preserving reducer). The partial
         * LLL state has no guarantee that the first 4 vectors are linearly
         * independent, and the input generators are only a sublattice basis
         * (not the union lattice), so neither is safe to return blindly. */
        rho = -1;
        goto cleanup;
    }
    /* Normal path: in-place zero handling kept b[zeta..d-1] nonzero. */
    int rho_nz = 0;
    for (int i = st->zeta; i < d; i++) {
        if (!ibz_vec_4_is_zero(&st->b[i])) {
            if (rho_nz != i - (int)st->zeta) {
                for (int c = 0; c < 4; c++)
                    ibz_swap(&st->b[st->zeta + rho_nz][c], &st->b[i][c]);
            }
            rho_nz++;
        }
    }
    rho = rho_nz;
    n = rho < out_capacity ? rho : out_capacity;
    if (rho == 4 || publish_lower_rank) {
        for (int i = 0; i < n; i++)
            for (int c = 0; c < 4; c++)
                ibz_copy(&output[i][c], &st->b[st->zeta + i][c]);
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
        output, out_capacity, input, d, alg, 0, 1);
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
                        int retain_first_lower_rank)
{
    int first_rho = -1;
    int result = -1;
#if defined(SQISIGN_ML2_PROFILE)
    ml2_atomic_profile_dimension_t *profile_dimension;
#endif

    if (output == NULL || input == NULL ||
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
                                 retain_first_lower_rank && attempt == 0);
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
                                       retain_first_lower_rank);
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
        output, out_capacity, input, d, alg, 0);
}

int
quat_ml2_mlll(ibz_vec_4_t *output,
              int out_capacity,
              const ibz_vec_4_t *input,
              int d,
              const quat_alg_t *alg)
{
    return ml2_retry_driver_native(
        output, out_capacity, input, d, alg, 1);
}
