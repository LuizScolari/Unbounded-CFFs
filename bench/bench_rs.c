/**
 * @file bench_rs.c
 * @brief RS (polynomial / monotone) growth benchmark.
 *
 * Measures one embedding step by (1) generating the three new blocks with the
 * library's own generate_new_cff_blocks_timed -- which is where all the
 * finite-field work lives -- and (2) placing (concatenating) only those new
 * blocks into the final matrix. The copy of the old matrix and every file
 * operation are left out of the timed region, matching "time of generation and
 * concatenation of the new parts".
 *
 * Two times are reported, mirroring the two ways the matrix can be filled:
 *
 *   Inv-CFF = setup + inverted evaluation index + fill + concatenation
 *   CFF     =                                     fill + concatenation
 *
 * i.e. the second one is the cost of the same growth when the inverted
 * evaluation index is already available; the difference between them is what
 * building that index costs.
 *
 * The old-block geometry is derived from the generated blocks themselves
 * (old_rows == rows of the old_new block, old_cols == columns of the new_old
 * block), exactly as embed_cff lays them out, so the placement is always
 * self-consistent and never depends on reading the old file.
 */

#include "bench.h"
#include "cff_builder.h"     /* generated_cffs, SET_BIT/GET_BIT, WORDS_FOR_BITS,
                                generate_new_cff_blocks_timed, rs_phase_times,
                                free_matrix */

#include <stdlib.h>

int bench_rs_step(const RSStep *s, double *t_inv_out, double *t_cff_out,
                  double *t_concat_out, long long *rows_out, long long *cols_out)
{
    long Fq[BENCH_MAX_STEPS], K[BENCH_MAX_STEPS];
    for (int i = 0; i < s->steps; i++) { Fq[i] = s->q[i]; K[i] = s->k[i]; }

    /* (1) Generation of the new blocks -- the dominant, genuinely-new work,
           split into setup / inverted index / fill by the library. */
    rs_phase_times ph = {0};
    generated_cffs nb = generate_new_cff_blocks_timed(s->construction, s->block,
                                                      s->d, Fq, K, s->steps, &ph);

    if (nb.cff_new == NULL && nb.cff_new_old == NULL && nb.cff_old_new == NULL)
        return -1;                   /* generation failed (bad parameters) */

    /* Geometry, exactly as embed_cff assembles it. The old block shares its row
       range with old_new and its column range with new_old. */
    long old_rows = nb.rows_old_new;
    long old_cols = nb.cols_new_old;

    long new_total_rows = old_rows + nb.rows_new_old;
    long width_top      = old_cols + nb.cols_old_new;
    long width_bottom   = nb.cols_new_old + nb.cols_new;
    long new_total_cols = (width_top > width_bottom) ? width_top : width_bottom;
    long words_per_row  = WORDS_FOR_BITS(new_total_cols);

    /* Allocation of the destination matrix is setup, not "new work": kept out
       of the timed region (it also zero-fills the old region). */
    uint64_t **final_cff = malloc((size_t)new_total_rows * sizeof(uint64_t *));
    if (!final_cff) {
        free_matrix(nb.cff_old_new, nb.rows_old_new);
        free_matrix(nb.cff_new_old, nb.rows_new_old);
        free_matrix(nb.cff_new,     nb.rows_new);
        return -1;
    }
    for (long i = 0; i < new_total_rows; i++)
        final_cff[i] = calloc((size_t)words_per_row, sizeof(uint64_t));

    /* (2) Concatenation of the NEW parts only (the old-matrix copy is skipped). */
    double t1 = now_sec();

    for (long i = 0; i < nb.rows_old_new; i++)
        for (long j = 0; j < nb.cols_old_new; j++)
            if (GET_BIT(nb.cff_old_new[i], j))
                SET_BIT(final_cff[i], old_cols + j);

    for (long i = 0; i < nb.rows_new_old; i++)
        for (long j = 0; j < nb.cols_new_old; j++)
            if (GET_BIT(nb.cff_new_old[i], j))
                SET_BIT(final_cff[old_rows + i], j);

    for (long i = 0; i < nb.rows_new; i++)
        for (long j = 0; j < nb.cols_new; j++)
            if (GET_BIT(nb.cff_new[i], j))
                SET_BIT(final_cff[old_rows + i], nb.cols_new_old + j);

    double t_concat = now_sec() - t1;

    for (long i = 0; i < new_total_rows; i++) free(final_cff[i]);
    free(final_cff);
    free_matrix(nb.cff_old_new, nb.rows_old_new);
    free_matrix(nb.cff_new_old, nb.rows_new_old);
    free_matrix(nb.cff_new,     nb.rows_new);

    if (t_inv_out)    *t_inv_out    = ph.setup + ph.index + ph.fill + t_concat;
    if (t_cff_out)    *t_cff_out    = ph.fill + t_concat;
    if (t_concat_out) *t_concat_out = t_concat;
    if (rows_out)     *rows_out     = new_total_rows;
    if (cols_out)     *cols_out     = new_total_cols;
    return 0;
}
