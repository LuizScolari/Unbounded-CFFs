#ifndef CRS_GENERATE_H
#define CRS_GENERATE_H

#include <stdint.h>
#include <stddef.h>

typedef struct {
    uint64_t *data;
    size_t    row_words;
    int       row_offset;
    int       q;
    int       col_start;
    int       col_end;
} BlockArg;

/* Fills one prime-power block of the matrix (thread body). Exposed so the
   benchmark harness (bench/) can time filling only the new columns of a block,
   reusing the exact same inner loop used by cff_generate. */
void *generate_block(void *arg);

#endif
