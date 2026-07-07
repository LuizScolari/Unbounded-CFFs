#include "cff_io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define IO_BUF_SIZE  (64 * 1024)   /* 64 KB I/O buffer        */
#define HEADER_MAX   2048          /* max header line length  */

/**
 * Writes a CFF to a text file.
 *
 * The first line is the header "d=<d> k=<k> S=<S> t=<t> n=<n>", optionally
 * followed by " chain=k0:S0,k1:S1,..." when the CFF was embedded. The matrix
 * follows as t lines of n '0'/'1' characters, emitted through a buffer.
 *
 * @param cff CFF to write.
 * @param filename Output file path.
 * @return 0 on success, -1 on I/O error.
 */
int write_cff(const CFF *cff, const char *filename)
{
    FILE *f = fopen(filename, "w");
    if (!f) { perror(filename); return -1; }

    /* Header line: mandatory fields + optional chain */
    fprintf(f, "d=%d k=%d S=%d t=%lld n=%lld",
            cff->d, cff->k, cff->S, cff->t, cff->n);
    if (cff->chain_len > 0) {
        fprintf(f, " chain=");
        for (int i = 0; i < cff->chain_len; i++) {
            if (i > 0) fputc(',', f);
            fprintf(f, "%d:%d", cff->chain_k[i], cff->chain_S[i]);
        }
    }
    fputc('\n', f);

    /* Matrix rows: bit-packed -> '0'/'1' characters, buffered */
    char buf[IO_BUF_SIZE];
    int  pos = 0;

    for (long long i = 0; i < cff->t; i++) {
        for (long long j = 0; j < cff->n; j++) {
            buf[pos++] = ((cff->rows[i][j / 64] >> (j % 64)) & 1) ? '1' : '0';
            if (pos == IO_BUF_SIZE) {
                if (fwrite(buf, 1, IO_BUF_SIZE, f) != IO_BUF_SIZE)
                    goto write_err;
                pos = 0;
            }
        }
        buf[pos++] = '\n';
        if (pos == IO_BUF_SIZE) {
            if (fwrite(buf, 1, IO_BUF_SIZE, f) != IO_BUF_SIZE)
                goto write_err;
            pos = 0;
        }
    }
    if (pos > 0 && fwrite(buf, 1, pos, f) != (size_t)pos)
        goto write_err;

    fclose(f);
    return 0;

write_err:
    perror(filename);
    fclose(f);
    return -1;
}

/**
 * Reads a CFF from a text file written by write_cff.
 *
 * Parses the header (mandatory d/k/S/t/n plus the optional chain), validates
 * that 0 < d < k < S and t,n > 0, then reads the matrix, ignoring any character
 * that is not '0' or '1'.
 *
 * @param filename Path of the file to read.
 * @return Newly allocated CFF on success, or NULL on error.
 */
CFF *read_cff(const char *filename)
{
    FILE *f = fopen(filename, "r");
    if (!f) { perror(filename); return NULL; }

    CFF *cff = malloc(sizeof(CFF));
    if (!cff) { fclose(f); return NULL; }
    cff->chain_len = 0;

    /* Read the entire header line (supports long chain= strings) */
    char header[HEADER_MAX];
    if (!fgets(header, HEADER_MAX, f)) {
        fprintf(stderr, "read_cff: empty file '%s'\n", filename);
        free(cff); fclose(f); return NULL;
    }

    /* Parse mandatory fields */
    if (sscanf(header, "d=%d k=%d S=%d t=%lld n=%lld",
               &cff->d, &cff->k, &cff->S, &cff->t, &cff->n) != 5) {
        fprintf(stderr, "read_cff: malformed header in '%s'\n", filename);
        free(cff); fclose(f); return NULL;
    }
    if (cff->d <= 0 || cff->k <= cff->d || cff->S <= cff->k ||
        cff->t <= 0 || cff->n <= 0) {
        fprintf(stderr, "read_cff: invalid parameters in '%s'\n", filename);
        free(cff); fclose(f); return NULL;
    }

    /* Parse optional chain= field */
    char *cp = strstr(header, "chain=");
    if (cp) {
        cp += 6;   /* skip "chain=" */
        while (cff->chain_len < MAX_CHAIN) {
            int kv, sv;
            if (sscanf(cp, "%d:%d", &kv, &sv) != 2) break;
            cff->chain_k[cff->chain_len] = kv;
            cff->chain_S[cff->chain_len] = sv;
            cff->chain_len++;
            char *comma = strchr(cp, ',');
            if (!comma) break;
            cp = comma + 1;
        }
    }

    cff->words_per_row = (cff->n + 63) / 64;

    cff->rows = calloc((size_t)cff->t, sizeof(uint64_t *));
    if (!cff->rows) { free(cff); fclose(f); return NULL; }

    for (long long i = 0; i < cff->t; i++) {
        cff->rows[i] = calloc((size_t)cff->words_per_row, sizeof(uint64_t));
        if (!cff->rows[i]) {
            for (long long p = 0; p < i; p++) free(cff->rows[p]);
            free(cff->rows); free(cff); fclose(f);
            return NULL;
        }
    }

    /* Read matrix: ignore everything that is not '0'/'1' */
    char   buf[IO_BUF_SIZE];
    long long row = 0, col = 0;

    while (row < cff->t) {
        size_t got = fread(buf, 1, IO_BUF_SIZE, f);
        if (got == 0) break;
        for (size_t p = 0; p < got && row < cff->t; p++) {
            char c = buf[p];
            if (c != '0' && c != '1') continue;
            if (c == '1')
                cff->rows[row][col / 64] |= 1ULL << (col % 64);
            if (++col == cff->n) { col = 0; row++; }
        }
    }

    if (row != cff->t) {
        fprintf(stderr,
                "read_cff: expected %lld rows but parsed %lld in '%s'\n",
                cff->t, row, filename);
        free_cff(cff); fclose(f); return NULL;
    }

    fclose(f);
    return cff;
}
