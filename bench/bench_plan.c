/**
 * @file bench_plan.c
 * @brief The benchmark plan: which CFF chains each table is built from.
 *
 * Every table is a list of blocks (chains). A chain is a sequence of CFFs where
 * the first one is generated from scratch and each of the following is an
 * embedding of the one before it -- exactly one block of a benchmark table.
 *
 * Tables
 * ------
 *  rs-embedding      RS polynomial, full q x q rows, d = floor((q-1)/k) (the
 *                    largest strength the full block supports). A step raises
 *                    the degree k, jumps to a larger sub-field q, or both.
 *  rs-monotone       RS monotone: d and k are fixed, only q grows; rows (dk+1)q.
 *  subset-embedding  Subset M(d,k,S): d fixed per chain, k and S both grow with
 *                    no further constraint.
 *  subset-nested     Subset nested (Sperner) family: d = 1 and k = floor(S/2)
 *                    at every step, i.e. the maximum number of columns C(S,k).
 *  subset-monotone   Subset monotone family: d and k fixed, only S grows.
 *  crs-embedding     Chinese Remainder Sieve: one chain per strength d, adding
 *                    one prime power per step with n always maximal.
 *
 * Editing the arrays below is all it takes to benchmark other families.
 */

#include "bench.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

/* binom() lives in the Subset module; declared here instead of including its
   header, which defines a CFF type that clashes with the CRS one. */
extern long long binom(int n, int k);

/* ------------------------------------------------------------ table metadata */

static const struct {
    const char *name;
    const char *construction;
    Kind        kind;
} TABLE_META[TBL_COUNT] = {
#define X(id, name, cons, kind) { name, cons, kind },
    BENCH_TABLES
#undef X
};

const char *bench_table_name(TableId id)         { return TABLE_META[id].name; }
const char *bench_table_construction(TableId id) { return TABLE_META[id].construction; }
Kind        bench_table_kind(TableId id)         { return TABLE_META[id].kind; }

/* ------------------------------------------------------------------ helpers */

double bench_matrix_bytes(long long t, long long n)
{
    long long words = (n + 63) / 64;
    return (double)t * (double)words * 8.0;
}

static long long ipow_ll(long base, int exp)
{
    long long r = 1;
    for (int i = 0; i < exp; i++) r *= base;
    return r;
}

/* Rows of an RS CFF over F_q with degree k and strength d, using the same
   minimum/full test as the library. */
static long long rs_rows(char block, int d, long q, long k)
{
    int full = !(d < (q - 1) / k);
    if (block == 'f' && full) return (long long)q * q;
    return (long long)(d * k + 1) * q;
}

/* Columns of an RS CFF: one per polynomial of degree at most k over F_q. */
static long long rs_cols(long q, long k) { return ipow_ll(q, (int)k + 1); }

/* Largest n with (prod q)^(1/d) >= n, matching cff_alloc's log-space test. */
static long crs_max_n(const int *q, int k, int d)
{
    long double logprod = 0.0L;
    for (int i = 0; i < k; i++) logprod += logl((long double)q[i]);
    long double target = logprod / (long double)d;
    long n = (long)floorl(expl(target) + 0.5L);
    while (n > 1 && (long double)d * logl((long double)n) > logprod + 1e-9L) n--;
    while ((long double)d * logl((long double)(n + 1)) <= logprod + 1e-9L)   n++;
    return n;
}

/* ------------------------------------------------------------- plan building */

typedef struct { BenchRow *rows; int n, max; } Plan;

static BenchRow *plan_add(Plan *p, TableId table, int block, int step)
{
    if (p->n >= p->max) return NULL;
    BenchRow *r = &p->rows[p->n++];
    memset(r, 0, sizeof(*r));
    r->table = table;
    r->block = block;
    r->step  = step;
    return r;
}

/**
 * One RS chain. qs/ks hold the parameters of every CFF of the chain; ds holds
 * the strength of each one. The first row is generated from scratch by the
 * polynomial builder (the monotone construction can only embed, never start a
 * chain), and every later row embeds the whole prefix of the chain.
 */
static void add_rs_chain(Plan *p, TableId table, int block, char construction,
                         char blk, const long *qs, const long *ks, const int *ds, int m)
{
    for (int i = 0; i < m; i++) {
        BenchRow *r = plan_add(p, table, block, i);
        if (!r) return;

        snprintf(r->q, sizeof(r->q), "%ld", qs[i]);
        r->k = (int)ks[i];
        r->d = ds[i];
        r->t = rs_rows(blk, ds[i], qs[i], ks[i]);
        r->n = rs_cols(qs[i], ks[i]);
        if (i > 0) {
            r->t_prev = rs_rows(blk, ds[i - 1], qs[i - 1], ks[i - 1]);
            r->n_prev = rs_cols(qs[i - 1], ks[i - 1]);
        }

        r->rs.construction = (i == 0) ? 'p' : construction;
        r->rs.block        = blk;
        r->rs.d            = ds[i];
        r->rs.steps        = i + 1;
        for (int j = 0; j <= i; j++) { r->rs.q[j] = qs[j]; r->rs.k[j] = ks[j]; }
    }
}

/** One Subset chain M(d, ks[i], Ss[i]) with d fixed along the chain. */
static void add_subset_chain(Plan *p, TableId table, int block, int d,
                             const int *ks, const int *Ss, int m)
{
    for (int i = 0; i < m; i++) {
        BenchRow *r = plan_add(p, table, block, i);
        if (!r) return;

        r->d = d;
        r->k = ks[i];
        r->S = Ss[i];
        r->t = binom(Ss[i], d);
        r->n = binom(Ss[i], ks[i]);
        if (i > 0) {
            r->t_prev = binom(Ss[i - 1], d);
            r->n_prev = binom(Ss[i - 1], ks[i - 1]);
        }

        r->sub.d         = d;
        r->sub.k_new     = ks[i];
        r->sub.S_new     = Ss[i];
        r->sub.chain_len = i;                    /* i == 0: initial, no chain */
        for (int j = 0; j < i; j++) { r->sub.chain_k[j] = ks[j]; r->sub.chain_S[j] = Ss[j]; }
        r->sub.n_old = (i > 0) ? binom(Ss[i - 1], ks[i - 1]) : 0;
    }
}

/**
 * One CRS chain: the initial CFF uses every prime power in base, then each step
 * adds one prime power from add. n is always the maximum the moduli allow.
 */
static void add_crs_chain(Plan *p, TableId table, int block, int d,
                          const int *base, int nbase, const int *add, int nadd)
{
    int       cur[BENCH_MAX_STEPS];
    int       ncur  = 0;
    long long t_cur = 0;

    BenchRow *r = plan_add(p, table, block, 0);
    if (!r) return;
    for (int j = 0; j < nbase; j++) { cur[ncur++] = base[j]; t_cur += base[j]; }

    int off = 0;
    for (int j = 0; j < nbase && off < (int)sizeof(r->q) - 8; j++)
        off += snprintf(r->q + off, sizeof(r->q) - off, "%s%d", j ? "," : "", base[j]);
    r->d = d;
    r->t = t_cur;
    r->n = crs_max_n(cur, ncur, d);

    r->crs.d    = d;
    r->crs.oldk = 0;                             /* initial: every column new */
    r->crs.addk = nbase;
    for (int j = 0; j < nbase; j++) r->crs.addq[j] = base[j];
    r->crs.old_n = 0;
    r->crs.new_n = r->n;

    for (int i = 0; i < nadd; i++) {
        r = plan_add(p, table, block, i + 1);
        if (!r) return;

        int all[BENCH_MAX_STEPS];
        for (int j = 0; j < ncur; j++) all[j] = cur[j];
        all[ncur] = add[i];

        snprintf(r->q, sizeof(r->q), "%d", add[i]);
        r->d      = d;
        r->t_prev = t_cur;
        r->n_prev = crs_max_n(cur, ncur, d);
        r->t      = t_cur + add[i];
        r->n      = crs_max_n(all, ncur + 1, d);

        r->crs.d    = d;
        r->crs.oldk = ncur;
        for (int j = 0; j < ncur; j++) r->crs.oldq[j] = cur[j];
        r->crs.addk  = 1;
        r->crs.addq[0] = add[i];
        r->crs.old_n = r->n_prev;
        r->crs.new_n = r->n;

        cur[ncur++] = add[i];
        t_cur      += add[i];
    }
}

/* ------------------------------------------------------------------- tables */

/* RS embedding: full q x q rows, d = floor((q-1)/k) (the maximum strength the
   full block supports). Each chain stays inside one characteristic, so every q
   is a sub-field of the next. */
static void build_rs_embedding(Plan *p)
{
    static const long q1[] = {2, 4, 4, 16, 16, 16}, k1[] = {1, 1, 2, 2, 3, 4};
    static const long q2[] = {3, 3, 9, 9, 9, 9},    k2[] = {1, 2, 2, 3, 4, 5};
    static const long q3[] = {5, 5, 25, 25},        k3[] = {1, 2, 2, 3};

    static const struct { const long *q, *k; int m; } chains[] = {
        { q1, k1, 6 }, { q2, k2, 6 }, { q3, k3, 4 },
    };

    for (int c = 0; c < (int)(sizeof(chains) / sizeof(chains[0])); c++) {
        int ds[BENCH_MAX_STEPS];
        for (int i = 0; i < chains[c].m; i++)
            ds[i] = (int)((chains[c].q[i] - 1) / chains[c].k[i]);
        add_rs_chain(p, TBL_RS_EMBEDDING, c + 1, 'p', 'f',
                     chains[c].q, chains[c].k, ds, chains[c].m);
    }
}

/* RS monotone: d and k fixed, only q grows (rows (dk+1)q, linear in q). */
static void build_rs_monotone(Plan *p)
{
    static const struct { int d; long k; long qs[BENCH_MAX_STEPS]; int m; } chains[] = {
        { 1, 1, {2, 4, 16, 256}, 4 },
        { 2, 1, {3, 9, 81},      3 },
        { 1, 2, {3, 27},         2 },
        { 2, 2, {5, 25},         2 },
    };

    for (int c = 0; c < (int)(sizeof(chains) / sizeof(chains[0])); c++) {
        long ks[BENCH_MAX_STEPS];
        int  ds[BENCH_MAX_STEPS];
        for (int i = 0; i < chains[c].m; i++) { ks[i] = chains[c].k; ds[i] = chains[c].d; }
        add_rs_chain(p, TBL_RS_MONOTONE, c + 1, 'm', 'm',
                     chains[c].qs, ks, ds, chains[c].m);
    }
}

/* Subset embedding: d fixed per chain, k and S both grow. */
static void build_subset_embedding(Plan *p)
{
    static const struct { int d; int ks[BENCH_MAX_STEPS], Ss[BENCH_MAX_STEPS]; int m; } chains[] = {
        { 1, {5,  8, 11, 13}, {10, 16, 22, 26}, 4 },
        { 2, {4,  7, 10},     {10, 16, 22},     3 },
        { 3, {5,  8, 10},     {12, 18, 22},     3 },
    };

    for (int c = 0; c < (int)(sizeof(chains) / sizeof(chains[0])); c++)
        add_subset_chain(p, TBL_SUBSET_EMBEDDING, c + 1, chains[c].d,
                         chains[c].ks, chains[c].Ss, chains[c].m);
}

/* Subset nested (Sperner): d = 1 and k = floor(S/2) at every step, which is the
   maximum number of columns C(S,k) for that S. */
static void build_subset_nested(Plan *p)
{
    static const int S1[] = {12, 16, 20, 24, 28};
    static const int S2[] = {11, 15, 19, 23, 27};

    static const struct { const int *Ss; int m; } chains[] = {
        { S1, 5 }, { S2, 5 },
    };

    for (int c = 0; c < (int)(sizeof(chains) / sizeof(chains[0])); c++) {
        int ks[BENCH_MAX_STEPS];
        for (int i = 0; i < chains[c].m; i++) ks[i] = chains[c].Ss[i] / 2;
        add_subset_chain(p, TBL_SUBSET_NESTED, c + 1, 1,
                         ks, chains[c].Ss, chains[c].m);
    }
}

/* Subset monotone: d and k fixed, only S grows. */
static void build_subset_monotone(Plan *p)
{
    static const struct { int d, k; int Ss[BENCH_MAX_STEPS]; int m; } chains[] = {
        { 1, 10, {20, 24, 28, 32}, 4 },
        { 2,  8, {16, 20, 24, 28}, 4 },
        { 3,  6, {12, 16, 20, 24}, 4 },
    };

    for (int c = 0; c < (int)(sizeof(chains) / sizeof(chains[0])); c++) {
        int ks[BENCH_MAX_STEPS];
        for (int i = 0; i < chains[c].m; i++) ks[i] = chains[c].k;
        add_subset_chain(p, TBL_SUBSET_MONOTONE, c + 1, chains[c].d,
                         ks, chains[c].Ss, chains[c].m);
    }
}

/* CRS: one chain per strength d; each step adds one prime power, n maximal. */
static void build_crs_embedding(Plan *p)
{
    static const struct {
        int d;
        int base[8]; int nbase;
        int add[16]; int nadd;
    } chains[] = {
        { 1, {2, 3},          2, {5, 7, 11, 13, 17, 19},              6 },
        { 2, {2, 3, 5, 7},    4, {11, 13, 17, 19, 23, 29, 31, 37},    8 },
        { 3, {2, 3, 5, 7, 11}, 5, {13, 17, 19, 23, 29, 31},           6 },
    };

    for (int c = 0; c < (int)(sizeof(chains) / sizeof(chains[0])); c++)
        add_crs_chain(p, TBL_CRS_EMBEDDING, c + 1, chains[c].d,
                      chains[c].base, chains[c].nbase,
                      chains[c].add,  chains[c].nadd);
}

int bench_plan_build(BenchRow *rows, int max)
{
    Plan p = { rows, 0, max };
    build_rs_embedding(&p);
    build_rs_monotone(&p);
    build_subset_embedding(&p);
    build_subset_nested(&p);
    build_subset_monotone(&p);
    build_crs_embedding(&p);
    return p.n;
}
