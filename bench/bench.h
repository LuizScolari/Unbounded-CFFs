#ifndef BENCH_H
#define BENCH_H

/**
 * @file bench.h
 * @brief Growth (embedding) benchmark harness for the three CFF constructions.
 *
 * The benchmark walks a set of CFF chains. Each chain starts with an initial
 * CFF generated from scratch and every following row is an embedding of the row
 * above it, so a chain is exactly one block of a benchmark table. For each row
 * the harness measures the wall-clock time to produce ONLY the NEW parts of the
 * grown matrix, using each construction's real inner loop. It never touches
 * files: every step is described purely by parameters, so no file read/write
 * time is ever counted.
 *
 * What "new parts" means per construction (see bench/README.md for the full
 * rationale):
 *   - RS  : the three new blocks (old_new, new_old, new_new) plus placing
 *           (concatenating) them into the final matrix. The copy of the old
 *           matrix and all file I/O are excluded. Reported twice: with and
 *           without the construction of the inverted evaluation index.
 *   - CRS : setting the 1s only in the new positions -- for the old prime
 *           powers, the new columns [n_old, n_new); for the new prime powers,
 *           all columns. The 1s of the old top-left block are not re-set here.
 *   - Sub : generating only the new columns [n_old, n_new). Old columns (and
 *           hence the old rows) are never processed.
 */

#include <stdint.h>

/* ------------------------------------------------------------------ timing */

/** Monotonic wall-clock time in seconds (CLOCK_MONOTONIC). */
double now_sec(void);

/* --------------------------------------------------------------- statistics */

typedef struct {
    double mean;
    double min;
    double max;
    double std;      /* sample standard deviation (0 when n < 2)             */
    double median;
    int    n;
} Stats;

/** Summary statistics of an array of timing samples. */
Stats stats_compute(const double *x, int n);

/* ----------------------------------------------------------- step descriptors */

#define BENCH_MAX_STEPS 32

/**
 * RS growth step. The family is described by the full chain of prime-power
 * sub-fields q[0..steps-1] and matching degrees k[0..steps-1]; the growth being
 * measured is the addition of the LAST entry, i.e. the "old" CFF spans
 * q[0..steps-2]. A step may raise q, k, or both. construction is 'p'
 * (polynomial) or 'm' (monotone); block is 'f' (full, q*q rows) or 'm'
 * (minimum). The old/new dimensions used for placement are derived from the
 * generated blocks themselves, so no dimension needs to be supplied here.
 */
typedef struct {
    char construction;               /* 'p' polynomial | 'm' monotone         */
    char block;                      /* 'f' full | 'm' minimum                */
    int  d;
    long q[BENCH_MAX_STEPS];
    long k[BENCH_MAX_STEPS];
    int  steps;
} RSStep;

/**
 * CRS growth step. Old prime powers oldq[0..oldk-1] carry old_n columns; addq
 * are the prime powers added; the grown matrix has new_n columns and strength d.
 * An initial CFF has oldk == 0 (every column of every block is new).
 */
typedef struct {
    int  d;
    int  oldq[BENCH_MAX_STEPS]; int oldk;
    int  addq[BENCH_MAX_STEPS]; int addk;
    long old_n;
    long new_n;
} CRSStep;

/**
 * Subset growth step. chain_k/chain_S are the predecessors (oldest first,
 * including the immediate predecessor whose (k,S) fixes n_old); the grown matrix
 * is M(d, k_new, S_new). An initial CFF has chain_len == 0.
 */
typedef struct {
    int       d;
    int       chain_k[BENCH_MAX_STEPS];
    int       chain_S[BENCH_MAX_STEPS];
    int       chain_len;
    int       k_new, S_new;
    long long n_old;                 /* = C(S_prev, k_prev)                   */
} SubsetStep;

/* --------------------------------------------------------------- per-step timers */

/**
 * Times one RS step.
 *
 * @param t_inv    Out: setup + inverted index + fill + concatenation (seconds).
 * @param t_cff    Out: fill + concatenation, i.e. the same work with the
 *                 inverted evaluation index already computed.
 * @param t_concat Out (optional): the concatenation part alone.
 * @param rows     Out (optional): rows of the assembled grown matrix.
 * @param cols     Out (optional): columns of the assembled grown matrix.
 * @return 0 on success, -1 on failure (bad parameters or allocation).
 */
int bench_rs_step(const RSStep *s, double *t_inv, double *t_cff, double *t_concat,
                  long long *rows, long long *cols);

/** Seconds to set the 1s of the new positions, or a negative value on failure. */
double bench_crs_step(const CRSStep *s);

/** Seconds to generate the new columns, or a negative value on failure. */
double bench_subset_step(const SubsetStep *s);

/* ------------------------------------------------------------- benchmark plan */

typedef enum { K_RS, K_CRS, K_SUB } Kind;

/** Table identifiers, in output order. */
#define BENCH_TABLES                                                          \
    X(RS_EMBEDDING,     "rs-embedding",     "RS",     K_RS)                   \
    X(RS_MONOTONE,      "rs-monotone",      "RS",     K_RS)                   \
    X(SUBSET_EMBEDDING, "subset-embedding", "Subset", K_SUB)                  \
    X(SUBSET_NESTED,    "subset-nested",    "Subset", K_SUB)                  \
    X(SUBSET_MONOTONE,  "subset-monotone",  "Subset", K_SUB)                  \
    X(CRS_EMBEDDING,    "crs-embedding",    "CRS",    K_CRS)

typedef enum {
#define X(id, name, cons, kind) TBL_##id,
    BENCH_TABLES
#undef X
    TBL_COUNT
} TableId;

/** Table id -> stable name ("rs-embedding", ...). */
const char *bench_table_name(TableId id);
/** Table id -> construction name ("RS" / "Subset" / "CRS"). */
const char *bench_table_construction(TableId id);
/** Table id -> which step timer applies. */
Kind        bench_table_kind(TableId id);

/**
 * One row of one benchmark table: the parameters shown in the table plus the
 * step descriptor that produces it. step == 0 is the initial CFF of the block
 * (generated from scratch); step > 0 is an embedding of the row above.
 */
typedef struct {
    TableId     table;
    int         block;               /* 1-based chain index inside the table  */
    int         step;                /* 0-based row index inside the chain    */

    /* Parameters displayed in the table (unused fields are 0 / empty). */
    char        q[64];               /* RS: field size; CRS: prime powers     */
    int         k;                   /* RS: degree;     Subset: subset size   */
    int         d;
    int         S;                   /* Subset only                           */
    long long   t;                   /* rows of this CFF                      */
    long long   n;                   /* columns of this CFF                   */
    long long   t_prev, n_prev;      /* predecessor dimensions (0 when initial)*/

    /* Step descriptor (only the member matching the table's kind is used). */
    RSStep      rs;
    CRSStep     crs;
    SubsetStep  sub;
} BenchRow;

#define BENCH_MAX_ROWS 256

/**
 * Fills @p rows with the full benchmark plan (all tables, all blocks, all
 * steps), in table order.
 *
 * @return Number of rows written, at most @p max.
 */
int bench_plan_build(BenchRow *rows, int max);

/** Bytes of a bit-packed t x n matrix. */
double bench_matrix_bytes(long long t, long long n);

#endif /* BENCH_H */
