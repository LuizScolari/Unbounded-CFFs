#ifndef GENERATE_CFF_H
#define GENERATE_CFF_H

#include <stdint.h>

#define MAX_CHAIN 64

typedef struct {
    int       d, k, S;
    long long t;
    long long n;
    long long words_per_row;
    uint64_t **rows;

    int chain_len;
    int chain_k[MAX_CHAIN];
    int chain_S[MAX_CHAIN];
} CFF;

long long binom(int n, int k);

void      colex_unrank(long long r, int k, int S, int *a);
long long colex_rank(const int *a, int k);

CFF *generate_cff(int d, int k, int S,
                  int chain_len, const int *chain_k, const int *chain_S);

void free_cff(CFF *cff);

#endif /* GENERATE_CFF_H */
