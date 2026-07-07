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


#endif
