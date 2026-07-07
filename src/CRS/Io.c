#include "cff.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#define SAVE_CHUNK 2048   /* rows per batch */

typedef struct {
    const uint64_t *data;
    size_t  row_words;
    int     n;
    char   *buf;       /* output text buffer for this slice */
    int     r_start;
    int     r_end;
} TextArg;

/**
 * Converts a slice of bit-packed rows into '0'/'1' text (thread body).
 *
 * Each row is written as n characters followed by '\n' into the shared buffer,
 * processing 64 columns per machine word.
 *
 * @param arg Pointer to a TextArg describing the row slice.
 * @return Always NULL.
 */
static void *rows_to_text(void *arg)
{
    TextArg *a = (TextArg *)arg;
    const int     n   = a->n;
    const size_t  rw  = a->row_words;
    const size_t  stride = (size_t)n + 1;   /* n chars + '\n' */

    for (int r = a->r_start; r < a->r_end; r++) {
        const uint64_t *row = a->data + (size_t)r * rw;
        char *out = a->buf + (size_t)(r - a->r_start) * stride;

        /* Process 64 columns at a time */
        int full_words = n >> 6;
        for (int w = 0; w < full_words; w++) {
            uint64_t word = row[w];
            char *p = out + (w << 6);
            for (int b = 0; b < 64; b++)
                p[b] = ((word >> b) & 1) + '0';
        }
        /* Remaining columns */
        int rem = n & 63;
        if (rem) {
            uint64_t word = row[full_words];
            char *p = out + (full_words << 6);
            for (int b = 0; b < rem; b++)
                p[b] = ((word >> b) & 1) + '0';
        }
        out[n] = '\n';
    }
    return NULL;
}

/**
 * Saves a CFF to a text file.
 *
 * The first line is the header "d n k q1 q2 ... qk"; the following t lines are
 * the matrix rows, each n characters of '0'/'1'. Rows are converted to text in
 * parallel chunks and written sequentially.
 *
 * @param cff CFF to save.
 * @param filename Output file path.
 * @return 0 on success, -1 on I/O or allocation error.
 */
int cff_save(const CFF *cff, const char *filename)
{
    FILE *f = fopen(filename, "w");
    if (!f) { perror(filename); return -1; }

    /* Header: d n k q1 q2 ... qk */
    fprintf(f, "%d %d %d", cff->d, cff->n, cff->k);
    for (int j = 0; j < cff->k; j++)
        fprintf(f, " %d", cff->q[j]);
    fprintf(f, "\n");

    /* Body: convert in chunks, write sequentially */
    int nthreads = 4;
    size_t stride = (size_t)cff->n + 1;

    for (int r = 0; r < cff->t; r += SAVE_CHUNK) {
        int chunk = cff->t - r;
        if (chunk > SAVE_CHUNK) chunk = SAVE_CHUNK;

        char *buf = malloc((size_t)chunk * stride);
        if (!buf) { fclose(f); return -1; }

        /* Split chunk among threads */
        int nt = nthreads;
        if (nt > chunk) nt = chunk;

        pthread_t *thr  = malloc((size_t)nt * sizeof(pthread_t));
        TextArg   *args = malloc((size_t)nt * sizeof(TextArg));
        int rows_per = chunk / nt;

        for (int i = 0; i < nt; i++) {
            int rs = r + i * rows_per;
            int re = (i == nt - 1) ? r + chunk : rs + rows_per;
            args[i] = (TextArg){
                .data      = cff->data,
                .row_words = cff->row_words,
                .n         = cff->n,
                .buf       = buf + (size_t)(rs - r) * stride,
                .r_start   = rs,
                .r_end     = re
            };
            pthread_create(&thr[i], NULL, rows_to_text, &args[i]);
        }
        for (int i = 0; i < nt; i++)
            pthread_join(thr[i], NULL);

        fwrite(buf, 1, (size_t)chunk * stride, f);

        free(args);
        free(thr);
        free(buf);
    }

    fclose(f);
    return 0;
}

/**
 * Reads only the header line of a CFF file to recover its parameters.
 *
 * Parses "d n k q1 q2 ... qk" and allocates the prime-power array; the caller
 * must free(*q).
 *
 * @param filename Path of the file to read.
 * @param k Output: number of prime powers.
 * @param q Output: newly allocated array of prime powers (caller frees).
 * @param n Output: number of columns.
 * @param d Output: CFF strength.
 * @return 0 on success, -1 on I/O or allocation error.
 */
int cff_load_header(const char *filename,
                    int *k, int **q, int *n, int *d)
{
    FILE *f = fopen(filename, "r");
    if (!f) { perror(filename); return -1; }

    /* Read first line */
    char line[8192];
    if (!fgets(line, sizeof(line), f)) {
        fclose(f);
        return -1;
    }
    fclose(f);

    /* Header: d n k q1 q2 ... qk */
    char *p   = line;
    char *end;
    *d = (int)strtol(p, &end, 10); p = end;
    *n = (int)strtol(p, &end, 10); p = end;
    *k = (int)strtol(p, &end, 10); p = end;

    *q = malloc((size_t)(*k) * sizeof(int));
    if (!*q) return -1;
    for (int j = 0; j < *k; j++) {
        (*q)[j] = (int)strtol(p, &end, 10);
        p = end;
    }

    return 0;
}
