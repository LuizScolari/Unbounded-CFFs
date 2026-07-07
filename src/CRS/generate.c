#include "generate.h"
#include "cff.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>

/**
 * Allocates and validates a CFF for the given prime powers.
 *
 * Checks the covering condition prod(q[j]) >= n^d (compared in log space to
 * avoid overflow), then allocates the bit-packed matrix with t = sum(q[j])
 * rows and n columns.
 *
 * @param k Number of prime powers.
 * @param q Array of prime powers q[0..k-1].
 * @param n Number of columns (items).
 * @param d Target CFF strength.
 * @return Newly allocated CFF, or NULL on invalid parameters or allocation failure.
 */
CFF *cff_alloc(int k, const int *q, int n, int d)
{
    /* Validate: prod(q[j]) >= n^d  using logarithms */
    double log_prod = 0.0;
    int    t        = 0;
    for (int j = 0; j < k; j++) {
        log_prod += log((double)q[j]);
        t        += q[j];
    }
    double log_nd = (double)d * log((double)n);
    if (log_prod < log_nd - 1e-9) {
        fprintf(stderr, "Error: product of prime powers (%.2f) < n^d (%.2f)\n",
                exp(log_prod), exp(log_nd));
        return NULL;
    }

    CFF *cff = calloc(1, sizeof(CFF));
    if (!cff) return NULL;

    cff->k         = k;
    cff->n         = n;
    cff->d         = d;
    cff->t         = t;
    cff->row_words = ((size_t)n + 63) / 64;

    cff->q = malloc((size_t)k * sizeof(int));
    if (!cff->q) { free(cff); return NULL; }
    memcpy(cff->q, q, (size_t)k * sizeof(int));

    size_t total = (size_t)cff->t * cff->row_words;
    cff->data = calloc(total, sizeof(uint64_t));
    if (!cff->data) { free(cff->q); free(cff); return NULL; }

    return cff;
}

/**
 * Frees a CFF and its owned buffers.
 *
 * @param cff CFF to free (may be NULL).
 */
void cff_free(CFF *cff)
{
    if (!cff) return;
    free(cff->q);
    free(cff->data);
    free(cff);
}

/**
 * Fills one prime-power block of the matrix (thread body).
 *
 * For prime power q at row_offset, row y (0..q-1) gets a 1 in every column c
 * with c ≡ y (mod q). Rows are filled one at a time so each row's contiguous
 * words stay cache-local.
 *
 * @param arg Pointer to a BlockArg describing the block to fill.
 * @return Always NULL.
 */
static void *generate_block(void *arg)
{
    BlockArg *a = (BlockArg *)arg;

    const int       q   = a->q;
    const size_t    rw  = a->row_words;
    const int       cs  = a->col_start;
    const int       ce  = a->col_end;
    uint64_t *const base = a->data + (size_t)a->row_offset * rw;

    for (int y = 0; y < q; y++) {
        uint64_t *row = base + (size_t)y * rw;
        /* first column congruent to y in [cs, ce) */
        int first = cs + ((y - cs % q) + q) % q;
        for (int c = first; c < ce; c += q) {
            row[c >> 6] |= 1ULL << (c & 63);
        }
    }
    return NULL;
}

/**
 * Generates the full CFF matrix in parallel, one thread per prime power.
 *
 * @param cff CFF to fill (already allocated).
 * @param nthreads Requested thread count (clamped to at least 1).
 */
void cff_generate(CFF *cff, int nthreads)
{
    if (nthreads < 1) nthreads = 1;

    int k = cff->k;

    pthread_t *thr  = malloc((size_t)k * sizeof(pthread_t));
    BlockArg  *args = malloc((size_t)k * sizeof(BlockArg));

    int offset = 0;
    for (int j = 0; j < k; j++) {
        args[j] = (BlockArg){
            .data       = cff->data,
            .row_words  = cff->row_words,
            .row_offset = offset,
            .q          = cff->q[j],
            .col_start  = 0,
            .col_end    = cff->n
        };
        offset += cff->q[j];
        pthread_create(&thr[j], NULL, generate_block, &args[j]);
    }

    for (int j = 0; j < k; j++)
        pthread_join(thr[j], NULL);

    free(thr);
    free(args);
}
