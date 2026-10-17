/**
 * @file bench_subset.c
 * @brief Subset-inclusion growth benchmark.
 *
 * The Subset embedding reorders columns so the predecessor occupies the
 * upper-left block; the new columns are therefore exactly the indices
 * [n_old, n_new). Because the algorithm is column-parallel, a column is the
 * natural atomic unit of work, and the growth cost is dominated by materialising
 * these new columns.
 *
 * This benchmark times only the new columns, using the same inner loop as
 * generate_cff (chain_col -> colex_unrank -> enumerate d-subsets -> colex_rank
 * -> set bit). It applies to every regime of the plan:
 *   - monotone (k' = k): old columns hold only old cells, so the new columns are
 *     exactly the new cells;
 *   - nested (k = floor(S/2)) and embedding (k grows freely): old columns are
 *     deliberately not processed (per the convention that old columns/rows are
 *     excluded), and each new column is materialised in full.
 *
 * Allocating the destination rows is done outside the timed region.
 */

#include "bench.h"
#include "generate_cff.h"     /* CFF, binom, colex_rank/unrank, chain_col,
                                 next_comb_init, next_comb */

#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

double bench_subset_step(const SubsetStep *s)
{
    const int d = s->d, k = s->k_new, S = s->S_new;

    long long t = binom(S, d);
    long long n = binom(S, k);
    long long words_per_row = (n + 63) / 64;
    long long n_old = s->n_old;

    /* Destination rows (untimed setup). */
    uint64_t **rows = calloc((size_t)t, sizeof(uint64_t *));
    if (!rows) return -1.0;
    for (long long i = 0; i < t; i++) {
        rows[i] = calloc((size_t)words_per_row, sizeof(uint64_t));
        if (!rows[i]) {
            for (long long p = 0; p < i; p++) free(rows[p]);
            free(rows);
            return -1.0;
        }
    }

    /* Local copies of the chain for the parallel region. */
    int ck[BENCH_MAX_STEPS], cs[BENCH_MAX_STEPS];
    for (int i = 0; i < s->chain_len; i++) { ck[i] = s->chain_k[i]; cs[i] = s->chain_S[i]; }
    int clen = s->chain_len;

    long long first_block = n_old / 64;         /* block containing column n_old */
    long long n_blocks    = (n + 63) / 64;

    double t0 = now_sec();

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
#endif
    for (long long b = first_block; b < n_blocks; b++) {
        long long j_start = b * 64;
        long long j_end   = j_start + 64;
        if (j_end > n) j_end = n;

        int *K   = malloc((size_t)k * sizeof(int));
        int *idx = malloc((size_t)(d + 2) * sizeof(int));   /* +2 sentinels */
        int *D   = malloc((size_t)d * sizeof(int));
        if (!K || !idx || !D) { free(K); free(idx); free(D); continue; }

        for (long long j = j_start; j < j_end; j++) {
            if (j < n_old) continue;            /* skip old columns in the boundary block */

            long long ci = chain_col(j, ck, cs, clen, k);
            colex_unrank(ci, k, S, K);

            next_comb_init(idx, d, k);
            do {
                for (int i = 0; i < d; i++) D[i] = K[idx[i]];
                long long row = colex_rank(D, d);
                rows[row][b] |= 1ULL << (j - j_start);
            } while (next_comb(idx, d));
        }

        free(K); free(idx); free(D);
    }

    double elapsed = now_sec() - t0;

    for (long long i = 0; i < t; i++) free(rows[i]);
    free(rows);
    return elapsed;
}
