#include "generate_cff.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/**
 * Computes the binomial coefficient C(n,k) via Manolopoulos' iterative method
 * (SIGCSE Bulletin 34(4), 2002 — Comb6).
 *
 *   C(n,k) = prod_{j=1..k} (n - j + 1) / j,   with k := min(k, n - k).
 *
 * Alternating one multiply and one divide keeps the running value equal to
 * C(n,j) at each step — always an integer — which avoids overflow on the
 * intermediate products.
 *
 * @param n Size of the set.
 * @param k Size of the subset.
 * @return C(n,k), or 0 when k < 0, n < 0 or k > n.
 */
long long binom(int n, int k)
{
    if (k < 0 || n < 0 || k > n) return 0;
    if (k > n - k) k = n - k;
    long long t = 1;
    for (int j = 1; j <= k; j++)
        t = t * (n - j + 1) / j;
    return t;
}

/**
 * Computes the colex rank of a k-subset of {0,...,S-1}.
 *
 *   rank({a_0 < a_1 < ... < a_{k-1}}) = sum_{i=0}^{k-1} C(a_i, i + 1).
 *
 * @param a Subset elements in strictly increasing order.
 * @param k Subset size.
 * @return Colex rank of the subset.
 */
long long colex_rank(const int *a, int k)
{
    long long r = 0;
    for (int i = 0; i < k; i++)
        r += binom(a[i], i + 1);
    return r;
}

/**
 * Recovers the k-subset of {0,...,S-1} with a given colex rank (inverse of
 * colex_rank).
 *
 * From the highest position down, picks the largest a_i such that
 * C(a_i, i + 1) <= remaining rank, then subtracts it off.
 *
 * @param r Colex rank to unrank.
 * @param k Subset size.
 * @param S Ground-set size.
 * @param a Output: subset elements in strictly increasing order.
 */
void colex_unrank(long long r, int k, int S, int *a)
{
    for (int i = k - 1; i >= 0; i--) {
        int c = i;
        while (c + 1 < S && binom(c + 1, i + 1) <= r)
            c++;
        a[i] = c;
        r -= binom(c, i + 1);
    }
}

/**
 * Maps an output column index of an embedded M(d,k,S) to the standard colex
 * rank of its k-subset.
 *
 * At the innermost chain level (k_base = ck[clen-1], S_base = cs[clen-1]):
 *   E     = C(S_base, k_base)                            inherited columns
 *   DELTA = sum_{i=0}^{m-1} C(S_base + i, k_base + i + 1),  m = k - k_base
 *           (colex-earlier k-subsets that do not contain the added elements)
 *
 * The output index j maps as:
 *   j < E            -> old rank = DELTA + chain_col(j, ..., k_base)  (recurse)
 *   E <= j < E+DELTA -> old rank = j - E   (shifted block)
 *   j >= E + DELTA   -> old rank = j       (identity)
 *
 * @param j Output column index.
 * @param ck Chain k values (oldest first).
 * @param cs Chain S values (oldest first).
 * @param clen Chain length (0 means standard colex ordering).
 * @param k Current column subset size.
 * @return Standard colex rank of column j's k-subset.
 */
static long long chain_col(long long j,
                            const int *ck, const int *cs, int clen,
                            int k)
{
    if (clen == 0) return j;   /* no chain: standard colex */

    int k_base = ck[clen - 1];
    int S_base = cs[clen - 1];
    long long E = binom(S_base, k_base);
    int m = k - k_base;
    long long DELTA = 0;
    for (int i = 0; i < m; i++)
        DELTA += binom(S_base + i, k_base + i + 1);

    if (j < E) {
        /* Column j comes from the predecessor; recurse to get its colex rank
           at the previous level, then shift it into the current level's range. */
        long long ci_base = chain_col(j, ck, cs, clen - 1, k_base);
        return DELTA + ci_base;
    }
    if (j < E + DELTA) return j - E;
    return j;
}

/**
 * Initialises the state for enumerating t-combinations of {0,...,n-1}
 * (Knuth, TAOCP 7.2.1.3, Algorithm L).
 *
 * c[0..t-1] holds the current combination; c[t] = n and c[t+1] = 0 are
 * sentinels used by next_comb to terminate the gap-search without an explicit
 * bounds check. The caller must allocate t+2 ints.
 *
 * @param c State array of at least t+2 ints.
 * @param t Combination size.
 * @param n Ground-set size.
 */
static void next_comb_init(int *c, int t, int n)
{
    for (int j = 0; j < t; j++) c[j] = j;   /* L1: c[j] = j (0-indexed) */
    c[t]   = n;
    c[t+1] = 0;
}

/**
 * Advances c[0..t-1] to the next combination in Algorithm L order.
 *
 * @param c State array previously set up by next_comb_init.
 * @param t Combination size.
 * @return 1 if a successor exists, 0 once all C(n,t) combinations are done.
 */
static int next_comb(int *c, int t)
{
    int j = 0;
    while (c[j] + 1 == c[j+1]) {            /* L3: find smallest j with a gap */
        c[j] = j;
        j++;
    }
    if (j >= t) return 0;                   /* L4: done */
    c[j]++;                                  /* L5: increase c[j] */
    return 1;
}

/**
 * Builds M(d,k,S): a d-CFF with t = C(S,d) rows and n = C(S,k) columns.
 *
 * Column j is the k-subset of colex rank chain_col(j, ...); its d-subsets index
 * the rows set to 1 (entry 1 iff the d-subset is contained in the k-subset).
 * When chain_len > 0, chain_col reorders the columns so each predecessor lands
 * in the upper-left block. Columns are processed in 64-wide bit-packed blocks,
 * one per OpenMP thread.
 *
 * @param d Strength (row subset size).
 * @param k Column subset size.
 * @param S Ground-set size.
 * @param chain_len Number of predecessors in the embedding chain.
 * @param chain_k Predecessor k values (oldest first), or NULL for a fresh matrix.
 * @param chain_S Predecessor S values (oldest first), or NULL for a fresh matrix.
 * @return Newly allocated CFF, or NULL on allocation failure.
 */
CFF *generate_cff(int d, int k, int S,
                  int chain_len, const int *chain_k, const int *chain_S)
{
    CFF *cff = malloc(sizeof(CFF));
    if (!cff) return NULL;

    cff->d = d;  cff->k = k;  cff->S = S;
    cff->t = binom(S, d);
    cff->n = binom(S, k);
    cff->words_per_row = (cff->n + 63) / 64;

    /* Store chain */
    cff->chain_len = (chain_len < MAX_CHAIN) ? chain_len : MAX_CHAIN;
    if (chain_len > 0 && chain_k && chain_S) {
        memcpy(cff->chain_k, chain_k, (size_t)cff->chain_len * sizeof(int));
        memcpy(cff->chain_S, chain_S, (size_t)cff->chain_len * sizeof(int));
    }

    cff->rows = calloc((size_t)cff->t, sizeof(uint64_t *));
    if (!cff->rows) { free(cff); return NULL; }

    for (long long i = 0; i < cff->t; i++) {
        cff->rows[i] = calloc((size_t)cff->words_per_row, sizeof(uint64_t));
        if (!cff->rows[i]) {
            for (long long p = 0; p < i; p++) free(cff->rows[p]);
            free(cff->rows);
            free(cff);
            return NULL;
        }
    }

    /* Local copies of chain arrays for use inside the parallel region */
    int local_ck[MAX_CHAIN], local_cs[MAX_CHAIN];
    memcpy(local_ck, cff->chain_k, (size_t)cff->chain_len * sizeof(int));
    memcpy(local_cs, cff->chain_S, (size_t)cff->chain_len * sizeof(int));
    int local_clen = cff->chain_len;

    long long n_blocks = (cff->n + 63) / 64;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (long long b = 0; b < n_blocks; b++) {
        long long j_start = b * 64;
        long long j_end   = j_start + 64;
        if (j_end > cff->n) j_end = cff->n;

        int *K   = malloc((size_t)k * sizeof(int));
        int *idx = malloc((size_t)(d + 2) * sizeof(int));   /* +2 sentinels for Algorithm L */
        int *D   = malloc((size_t)d * sizeof(int));
        if (!K || !idx || !D) { free(K); free(idx); free(D); continue; }

        for (long long j = j_start; j < j_end; j++) {
            long long ci = chain_col(j, local_ck, local_cs, local_clen, k);
            colex_unrank(ci, k, S, K);

            next_comb_init(idx, d, k);
            do {
                for (int i = 0; i < d; i++) D[i] = K[idx[i]];

                long long row = colex_rank(D, d);
                cff->rows[row][b] |= 1ULL << (j - j_start);
            } while (next_comb(idx, d));
        }

        free(K); free(idx); free(D);
    }

    return cff;
}

/**
 * Frees a CFF and all of its rows.
 *
 * @param cff CFF to free (may be NULL).
 */
void free_cff(CFF *cff)
{
    if (!cff) return;
    for (long long i = 0; i < cff->t; i++)
        free(cff->rows[i]);
    free(cff->rows);
    free(cff);
}
