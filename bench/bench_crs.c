/**
 * @file bench_crs.c
 * @brief CRS (Chinese Remainder Sieve) growth benchmark.
 *
 * The CRS embedding never reads the old matrix: it regenerates the whole matrix
 * from the merged parameters, which re-sets the 1s of the old top-left block as
 * well. This benchmark times only the 1s of the NEW positions:
 *
 *   - old prime-power blocks (j < old_k): only the new columns [old_n, new_n);
 *   - new prime-power blocks (j >= old_k): all columns [0, new_n).
 *
 * It reuses generate_block -- the exact per-block inner loop of cff_generate --
 * with a per-block col_start, so what is timed is identical to the real fill,
 * just restricted to the new positions. Allocation (which zero-fills the whole
 * matrix, old region included) is done outside the timed region.
 */

#include "bench.h"
#include "cff.h"
#include "generate.h"

#include <stdlib.h>
#include <string.h>
#include <pthread.h>

double bench_crs_step(const CRSStep *s)
{
    int new_k = s->oldk + s->addk;

    int *q = malloc((size_t)new_k * sizeof(int));
    if (!q) return -1.0;
    memcpy(q,           s->oldq, (size_t)s->oldk * sizeof(int));
    memcpy(q + s->oldk, s->addq, (size_t)s->addk * sizeof(int));

    /* Full grown matrix; the calloc zero-fills old+new region (untimed setup). */
    CFF *cff = cff_alloc(new_k, q, (int)s->new_n, s->d);
    free(q);
    if (!cff) return -1.0;

    pthread_t *thr  = malloc((size_t)new_k * sizeof(pthread_t));
    BlockArg  *args = malloc((size_t)new_k * sizeof(BlockArg));
    if (!thr || !args) { free(thr); free(args); cff_free(cff); return -1.0; }

    int offset = 0;
    double t0 = now_sec();
    for (int j = 0; j < new_k; j++) {
        /* Old prime powers only fill the newly added columns; new prime powers
           fill every column. */
        int col_start = (j < s->oldk) ? (int)s->old_n : 0;
        args[j] = (BlockArg){
            .data       = cff->data,
            .row_words  = cff->row_words,
            .row_offset = offset,
            .q          = cff->q[j],
            .col_start  = col_start,
            .col_end    = cff->n
        };
        offset += cff->q[j];
        pthread_create(&thr[j], NULL, generate_block, &args[j]);
    }
    for (int j = 0; j < new_k; j++)
        pthread_join(thr[j], NULL);
    double t = now_sec() - t0;

    free(thr);
    free(args);
    cff_free(cff);
    return t;
}
