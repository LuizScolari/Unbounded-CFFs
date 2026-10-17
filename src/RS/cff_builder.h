#ifndef CFF_BUILDER_H
#define CFF_BUILDER_H

#include <stdint.h>
#include "flint/fq_nmod.h"
#include "flint/fq_nmod_poly.h"

#define BITS_PER_WORD 64
#define WORD_OFFSET(b) ((b) / BITS_PER_WORD)
#define BIT_OFFSET(b)  ((b) % BITS_PER_WORD)
#define SET_BIT(row_ptr, col) ((row_ptr)[WORD_OFFSET(col)] |= (1ULL << BIT_OFFSET(col)))
#define GET_BIT(row_ptr, col) (((row_ptr)[WORD_OFFSET(col)] >> BIT_OFFSET(col)) & 1ULL)
#define WORDS_FOR_BITS(bits) (((bits) + BITS_PER_WORD - 1) / BITS_PER_WORD)

typedef struct {
    long q;
    fq_nmod_t* all_elements;
    long count_all;
    long capacity_all;
    fq_nmod_t* only_elements;
    long count_only;
    long capacity_only;
} subfield_partition;

typedef struct {
    fq_nmod_t x;
    fq_nmod_t y;
} element_pair;

typedef struct {
    element_pair* combos_old;
    long count_old;
    long capacity_old;
    element_pair* combos_new;
    long count_new;
    long capacity_new;
} combination_partitions;

typedef struct {
    uint64_t** cff_old_new;
    long rows_old_new;
    long cols_old_new;

    uint64_t** cff_new_old;
    long rows_new_old;
    long cols_new_old;

    uint64_t** cff_new;
    long rows_new;
    long cols_new;
} generated_cffs;

typedef struct {
    fq_nmod_poly_t* old_polys;
    long num_old_polys;
    fq_nmod_poly_t* new_polys;
    long num_new_polys;
    fq_nmod_poly_t* all_polys;
    long num_all_polys;
} polynomial_partition;

void embed_cff(char construction, char block_size, const char *cff_file, int d, long Fq, long k);
void rs_generate_cff(char construction, char block_size, int d, long fq, long k);

/* Wall-clock breakdown of one generate_new_cff_blocks call, in seconds.
   setup: sub-field partition, polynomial partition and (x,y) combinations;
   index: the inverted evaluation indexes (old + new polynomials);
   fill : the three generate_single_cff calls that set the bits. */
typedef struct {
    double setup;
    double index;
    double fill;
} rs_phase_times;

/* Exposed for the benchmark harness (bench/): generates only the new CFF blocks
   (old_new, new_old, new_new) of one embedding step, and frees a bit-packed
   matrix. These are the same functions used internally by embed_cff. The
   _timed variant additionally reports the per-phase wall-clock split (pass NULL
   to skip); generate_new_cff_blocks is exactly that call with times == NULL. */
generated_cffs generate_new_cff_blocks(char construction, char block_size, int d,
                                       long *Fq_steps, long *k_steps, int num_steps);
generated_cffs generate_new_cff_blocks_timed(char construction, char block_size, int d,
                                             long *Fq_steps, long *k_steps, int num_steps,
                                             rs_phase_times *times);
void free_matrix(uint64_t **matrix, long rows);

#endif /* CFF_BUILDER_H */
