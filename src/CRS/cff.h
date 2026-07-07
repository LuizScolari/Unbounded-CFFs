#ifndef CFF_H
#define CFF_H

#include <stdint.h>
#include <stddef.h>

typedef struct {
    int      k;
    int     *q;
    int      n;
    int      d;
    int      t;
    size_t   row_words;
    uint64_t *data;
} CFF;

CFF  *cff_alloc(int k, const int *q, int n, int d);
void  cff_generate(CFF *cff, int nthreads);
void  cff_free(CFF *cff);

int   cff_save(const CFF *cff, const char *dir);
int   cff_load_header(const char *filename,
                      int *k, int **q, int *n, int *d);

CFF  *cff_embed(const char *old_file,
                int add_count, const int *add_q,
                int new_n, int new_d,
                int nthreads);

#endif /* CFF_H */
