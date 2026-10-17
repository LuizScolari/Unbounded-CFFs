/**
 * @file bench_main.c
 * @brief Driver for the CFF growth benchmark.
 *
 * Runs every chain of the benchmark plan (bench_plan.c) and reports, for each
 * CFF of a chain, the wall-clock time to produce ONLY the new parts of that
 * matrix: the whole matrix for the first CFF of a chain (generated from
 * scratch) and just the new blocks/columns for each embedding after it. File
 * I/O is never involved: every step is described by parameters, and the timed
 * region excludes reading/writing.
 *
 * The output is one table per construction/regime, each one a list of blocks
 * (chains) whose first row is an initial CFF and whose remaining rows embed the
 * row above:
 *
 *   rs-embedding      RS polynomial, full q x q rows, d = floor((q-1)/k).
 *   rs-monotone       RS monotone: d, k fixed, only q grows.
 *   subset-embedding  Subset: d fixed, k and S grow.
 *   subset-nested     Subset nested (Sperner): d = 1, k = floor(S/2).
 *   subset-monotone   Subset monotone: d, k fixed, S grows.
 *   crs-embedding     CRS: one prime power added per step, n maximal.
 *
 * RS reports two times: Inv-CFF(s), the total including the construction of the
 * inverted evaluation index, and CFF(s), the same growth with that index
 * already computed. The other constructions fill the matrix directly and report
 * a single CFF(s).
 *
 * Usage: ./cff_bench [--reps N] [--warmup W] [--budget S] [--only ID]
 *                    [--csv FILE] [--max-bytes B] [--list]
 */

#include "bench.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int       ran;              /* 1 = measured, 0 = skipped/failed         */
    Stats     inv;              /* RS only: index + fill + concatenation    */
    Stats     cff;              /* fill + concatenation (RS) / new parts    */
    double    concat;           /* RS only: mean concatenation time         */
    long long got_t, got_n;     /* RS only: geometry the library produced   */
    int       geom_mismatch;    /* RS only: got_t/got_n differ from the plan */
} Result;

typedef enum { RUN_OK, RUN_SKIP, RUN_FAIL } RunStatus;

static BenchRow g_plan[BENCH_MAX_ROWS];
static Result   g_res[BENCH_MAX_ROWS];
static int      g_n;

/* ------------------------------------------------------------------ helpers */

static long long parse_bytes(const char *s)
{
    char *end;
    double v = strtod(s, &end);
    if      (*end == 'K' || *end == 'k') v *= 1024.0;
    else if (*end == 'M' || *end == 'm') v *= 1024.0 * 1024.0;
    else if (*end == 'G' || *end == 'g') v *= 1024.0 * 1024.0 * 1024.0;
    return (long long)v;
}

/* Does this table match the --only selector? Accepts a table name, a
   construction group (rs / subset / crs) or "all". */
static int table_selected(TableId id, const char *only)
{
    if (!only || strcmp(only, "all") == 0) return 1;
    const char *name = bench_table_name(id);
    if (strcmp(only, name) == 0) return 1;
    if (strcmp(only, "rs")     == 0 && strncmp(name, "rs-",     3) == 0) return 1;
    if (strcmp(only, "subset") == 0 && strncmp(name, "subset-", 7) == 0) return 1;
    if (strcmp(only, "crs")    == 0 && strncmp(name, "crs-",    4) == 0) return 1;
    return 0;
}

/* "1048576" -> "1,048,576" (thousands separators, as in the paper tables). */
static void group_ll(char *out, size_t cap, long long v)
{
    char raw[32];
    snprintf(raw, sizeof(raw), "%lld", v);
    int len = (int)strlen(raw), o = 0;
    for (int i = 0; i < len && o < (int)cap - 2; i++) {
        if (i > 0 && (len - i) % 3 == 0) out[o++] = ',';
        out[o++] = raw[i];
    }
    out[o] = '\0';
}

/* ------------------------------------------------------------- measurement */

/**
 * Runs one plan row: warm-up runs (discarded, to absorb one-time costs such as
 * FLINT/GLib context init, thread-pool spin-up and first-touch page faults),
 * then repetitions until either @p reps samples are collected or @p budget
 * seconds of wall time have been spent (always at least one sample).
 */
static RunStatus run_row(const BenchRow *r, Result *out, int reps, int warmup,
                         double budget, long long max_bytes)
{
    Kind kind = bench_table_kind(r->table);
    double est_peak = bench_matrix_bytes(r->t, r->n);
    if (kind == K_RS) {
        /* For RS the matrix is the small part: the inverted evaluation index
           stores one polynomial index per (evaluation point, polynomial) pair
           inside GLib arrays and hash tables, which measures at roughly 70
           bytes per entry. That, plus the new blocks held alongside the
           assembled matrix, is what actually decides whether a row fits. */
        long q_last = r->rs.q[r->rs.steps - 1];
        est_peak = est_peak * 2.2 + 70.0 * (double)q_last * (double)r->n;
    }

    memset(out, 0, sizeof(*out));

    if (est_peak > (double)max_bytes) return RUN_SKIP;

    /* Warm-up runs (discarded): they absorb one-time costs -- FLINT/GLib
       context init, OpenMP/pthread pool spin-up, first-touch page faults --
       which otherwise dominate the sub-millisecond rows. They are capped at a
       tenth of the budget so that an expensive row still gets exactly one. */
    double warm_cap = budget / 10.0, warm_spent = 0.0;
    for (int w = 0; w < warmup; w++) {
        double a = 0.0, b = 0.0, c = 0.0, t0 = now_sec();
        int bad = 0;
        switch (kind) {
            case K_RS:  bad = bench_rs_step(&r->rs, &a, &b, &c, NULL, NULL) < 0; break;
            case K_CRS: bad = bench_crs_step(&r->crs) < 0.0;                     break;
            case K_SUB: bad = bench_subset_step(&r->sub) < 0.0;                  break;
        }
        if (bad) return RUN_FAIL;
        warm_spent += now_sec() - t0;
        if (warm_spent >= warm_cap) break;
    }

    double *s_inv = malloc((size_t)reps * sizeof(double));
    double *s_cff = malloc((size_t)reps * sizeof(double));
    if (!s_inv || !s_cff) { free(s_inv); free(s_cff); return RUN_FAIL; }

    double concat_sum = 0.0, spent = 0.0;
    int    good = 0, failed = 0;

    for (int i = 0; i < reps; i++) {
        double t_inv = 0.0, t_cff = 0.0, t_cat = 0.0;
        double t0 = now_sec();
        switch (kind) {
            case K_RS:
                if (bench_rs_step(&r->rs, &t_inv, &t_cff, &t_cat,
                                  &out->got_t, &out->got_n) < 0) failed = 1;
                break;
            case K_CRS:
                t_cff = bench_crs_step(&r->crs);
                if (t_cff < 0.0) failed = 1;
                break;
            case K_SUB:
                t_cff = bench_subset_step(&r->sub);
                if (t_cff < 0.0) failed = 1;
                break;
        }
        if (failed) break;

        s_inv[good] = t_inv;
        s_cff[good] = t_cff;
        good++;
        concat_sum += t_cat;

        spent += now_sec() - t0;
        if (spent >= budget) break;          /* keep the whole run affordable */
    }

    if (failed || good == 0) {
        free(s_inv); free(s_cff);
        return RUN_FAIL;
    }

    out->cff    = stats_compute(s_cff, good);
    out->inv    = stats_compute(s_inv, good);
    out->concat = concat_sum / good;
    out->ran    = 1;
    free(s_inv); free(s_cff);

    /* The RS geometry is assembled by the library, so it doubles as a check of
       the (t,n) written in the table. */
    out->geom_mismatch = (kind == K_RS && (out->got_t != r->t || out->got_n != r->n));
    return RUN_OK;
}

/* ------------------------------------------------------------------ output */

static const char *table_title(TableId id)
{
    switch (id) {
        case TBL_RS_EMBEDDING:     return "RS embedding -- polynomial, full q x q rows, d = floor((q-1)/k)";
        case TBL_RS_MONOTONE:      return "RS monotone -- d, k fixed, only q grows (rows (dk+1)q)";
        case TBL_SUBSET_EMBEDDING: return "Subset embedding -- d fixed, k and S grow";
        case TBL_SUBSET_NESTED:    return "Subset nested (Sperner) -- d = 1, k = floor(S/2)";
        case TBL_SUBSET_MONOTONE:  return "Subset monotone -- d, k fixed, S grows";
        case TBL_CRS_EMBEDDING:    return "CRS embedding -- one prime power added per step, n maximal";
        default:                   return "";
    }
}

static void print_table(TableId id)
{
    int any = 0;
    for (int i = 0; i < g_n; i++)
        if (g_plan[i].table == id && g_res[i].ran) { any = 1; break; }
    if (!any) return;

    printf("\n%s\n%s\n", bench_table_name(id), table_title(id));

    Kind kind = bench_table_kind(id);
    if (kind == K_RS)
        printf("%8s %4s %4s %12s %14s %6s %13s %13s\n",
               "q", "k", "d", "t", "n", "reps", "Inv-CFF(s)", "CFF(s)");
    else if (kind == K_SUB)
        printf("%8s %4s %4s %12s %14s %6s %13s\n",
               "S", "k", "d", "t", "n", "reps", "CFF(s)");
    else
        printf("%16s %4s %12s %14s %6s %13s\n",
               "q", "d", "t", "n", "reps", "CFF(s)");

    int prev_block = -1;
    for (int i = 0; i < g_n; i++) {
        const BenchRow *r = &g_plan[i];
        const Result   *x = &g_res[i];
        if (r->table != id || !x->ran) continue;
        if (prev_block != -1 && r->block != prev_block) printf("\n");
        prev_block = r->block;

        char t_s[32], n_s[32];
        group_ll(t_s, sizeof(t_s), r->t);
        group_ll(n_s, sizeof(n_s), r->n);

        if (kind == K_RS)
            printf("%8s %4d %4d %12s %14s %6d %13.6f %13.6f\n",
                   r->q, r->k, r->d, t_s, n_s, x->cff.n, x->inv.mean, x->cff.mean);
        else if (kind == K_SUB)
            printf("%8d %4d %4d %12s %14s %6d %13.6f\n",
                   r->S, r->k, r->d, t_s, n_s, x->cff.n, x->cff.mean);
        else
            printf("%16s %4d %12s %14s %6d %13.6f\n",
                   r->q, r->d, t_s, n_s, x->cff.n, x->cff.mean);
    }
}

static void write_csv(const char *path)
{
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); return; }

    fprintf(f, "table,construction,block,step,q,k,d,S,t,n,t_prev,n_prev,reps,ran,"
               "inv_mean_s,inv_min_s,inv_std_s,inv_median_s,"
               "cff_mean_s,cff_min_s,cff_std_s,cff_median_s,"
               "concat_mean_s,new_cells,matrix_MB\n");

    for (int i = 0; i < g_n; i++) {
        const BenchRow *r = &g_plan[i];
        const Result   *x = &g_res[i];
        Kind kind = bench_table_kind(r->table);

        fprintf(f, "%s,%s,%d,%d,\"%s\",", bench_table_name(r->table),
                bench_table_construction(r->table), r->block, r->step, r->q);
        /* k / d / S are only meaningful for some constructions */
        if (kind == K_CRS) fprintf(f, ",%d,,", r->d);
        else if (kind == K_SUB) fprintf(f, "%d,%d,%d,", r->k, r->d, r->S);
        else fprintf(f, "%d,%d,,", r->k, r->d);

        fprintf(f, "%lld,%lld,%lld,%lld,%d,%d,",
                r->t, r->n, r->t_prev, r->n_prev, x->cff.n, x->ran);

        if (kind == K_RS)
            fprintf(f, "%.9f,%.9f,%.9f,%.9f,",
                    x->inv.mean, x->inv.min, x->inv.std, x->inv.median);
        else
            fprintf(f, ",,,,");

        fprintf(f, "%.9f,%.9f,%.9f,%.9f,%.9f,%lld,%.3f\n",
                x->cff.mean, x->cff.min, x->cff.std, x->cff.median, x->concat,
                r->t * r->n - r->t_prev * r->n_prev,
                bench_matrix_bytes(r->t, r->n) / (1024.0 * 1024.0));
    }
    fclose(f);
    printf("\nCSV written to %s\n", path);
}

static void list_plan(void)
{
    TableId prev = TBL_COUNT;
    for (int i = 0; i < g_n; i++) {
        const BenchRow *r = &g_plan[i];
        if (r->table != prev) {
            printf("\n%s (%s)\n", bench_table_name(r->table), table_title(r->table));
            prev = r->table;
        }
        printf("  block %d step %d: %d-CFF(%lld,%lld)  ~%.1f MB%s\n",
               r->block, r->step, r->d, r->t, r->n,
               bench_matrix_bytes(r->t, r->n) / (1024.0 * 1024.0),
               r->step == 0 ? "   [initial]" : "");
    }
    printf("\n%d rows in the plan.\n", g_n);
}

/* ------------------------------------------------------------------- main */

static void usage(const char *prog)
{
    printf("Usage: %s [--reps N] [--warmup W] [--budget S] [--only ID] [--csv FILE]\n"
           "          [--max-bytes B] [--list]\n"
           "  --reps N       maximum repetitions per row (default 100)\n"
           "  --warmup W     discarded warm-up runs per row (default 5; they stop\n"
           "                 early after budget/10 seconds, so heavy rows get one)\n"
           "  --budget S     stop repeating a row after S seconds (default 20; 0 = no limit)\n"
           "  --only ID      rs-embedding | rs-monotone | rs | subset-embedding |\n"
           "                 subset-nested | subset-monotone | subset | crs-embedding |\n"
           "                 crs | all (default all)\n"
           "  --csv FILE     write the results as CSV (for bench/gen_tables.py)\n"
           "  --max-bytes B  skip rows whose estimated peak memory exceeds B\n"
           "                 (K/M/G ok; default 4G)\n"
           "  --list         list the plan without running it\n", prog);
}

int main(int argc, char **argv)
{
    int         reps      = 100;
    int         warmup    = 5;
    double      budget    = 20.0;
    const char *only      = "all";
    const char *csv       = NULL;
    long long   max_bytes = 4LL * 1024 * 1024 * 1024;
    int         list_only = 0;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--reps")      && i + 1 < argc) reps   = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--warmup")    && i + 1 < argc) warmup = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--budget")    && i + 1 < argc) budget = atof(argv[++i]);
        else if (!strcmp(argv[i], "--only")      && i + 1 < argc) only   = argv[++i];
        else if (!strcmp(argv[i], "--csv")       && i + 1 < argc) csv    = argv[++i];
        else if (!strcmp(argv[i], "--max-bytes") && i + 1 < argc) max_bytes = parse_bytes(argv[++i]);
        else if (!strcmp(argv[i], "--list"))                      list_only = 1;
        else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) { usage(argv[0]); return 0; }
        else { fprintf(stderr, "Unknown option: %s\n", argv[i]); usage(argv[0]); return 1; }
    }
    if (reps   < 1)   reps = 1;
    if (warmup < 0)   warmup = 0;
    if (budget <= 0.0) budget = 1e18;

    /* Build the plan, then drop the tables --only excludes. */
    int all = bench_plan_build(g_plan, BENCH_MAX_ROWS);
    g_n = 0;
    for (int i = 0; i < all; i++)
        if (table_selected(g_plan[i].table, only)) g_plan[g_n++] = g_plan[i];

    if (g_n == 0) { fprintf(stderr, "No table matches --only %s\n", only); return 1; }

    if (list_only) { list_plan(); return 0; }

    printf("CFF growth benchmark  (reps<=%d, warmup=%d, budget=%.0fs/row, cap=%.0f MB, only=%s)\n",
           reps, warmup, budget, (double)max_bytes / (1024.0 * 1024.0), only);

    TableId prev = TBL_COUNT;
    for (int i = 0; i < g_n; i++) {
        if (g_plan[i].table != prev) {
            prev = g_plan[i].table;
            printf("\n== %s ==\n", bench_table_name(prev));
        }
        printf("  run   block %d step %d: %d-CFF(%lld,%lld) ...",
               g_plan[i].block, g_plan[i].step, g_plan[i].d, g_plan[i].t, g_plan[i].n);
        fflush(stdout);

        switch (run_row(&g_plan[i], &g_res[i], reps, warmup, budget, max_bytes)) {
            case RUN_OK:
                printf(" %d rep%s, %.6f s\n", g_res[i].cff.n,
                       g_res[i].cff.n == 1 ? "" : "s", g_res[i].cff.mean);
                if (g_res[i].geom_mismatch)
                    printf("  WARN  the plan says %lld x %lld, the library built %lld x %lld\n",
                           g_plan[i].t, g_plan[i].n, g_res[i].got_t, g_res[i].got_n);
                break;
            case RUN_SKIP:
                printf(" SKIPPED (~%.0f MB over the cap)\n",
                       bench_matrix_bytes(g_plan[i].t, g_plan[i].n) / (1024.0 * 1024.0));
                break;
            case RUN_FAIL:
                printf(" FAILED\n");
                break;
        }
    }

    printf("\n==================================================================================\n");
    printf(" GROWTH BENCHMARK -- time to build ONLY the new parts (file I/O excluded)\n");
    printf("==================================================================================\n");
    for (TableId id = 0; id < TBL_COUNT; id++) print_table(id);

    if (csv) write_csv(csv);
    return 0;
}
