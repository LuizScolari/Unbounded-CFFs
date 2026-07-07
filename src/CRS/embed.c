#include "cff.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Extends an existing CFF with new prime powers and/or a larger n.
 *
 * Because the construction is deterministic, the extended matrix is fully
 * defined by the merged parameter set: only the old file's header is read, its
 * prime powers are merged with the new ones, and the full matrix is generated.
 * The top-left t×n submatrix equals the old matrix by construction, so no data
 * is copied.
 *
 * @param old_file Path of the existing CFF file (header is read).
 * @param add_count Number of prime powers to add.
 * @param add_q Array of prime powers to add (must be disjoint from the old ones).
 * @param new_n New number of columns (must be >= old n).
 * @param new_d New strength (must be >= old d).
 * @param nthreads Thread count for parallel generation.
 * @return Newly generated extended CFF, or NULL on error.
 */
CFF *cff_embed(const char *old_file,
               int add_count, const int *add_q,
               int new_n, int new_d,
               int nthreads)
{
    /* 1. Read old parameters */
    int  old_k, old_n, old_d;
    int *old_q = NULL;

    if (cff_load_header(old_file, &old_k, &old_q, &old_n, &old_d) != 0) {
        fprintf(stderr, "Error: cannot read header from '%s'\n", old_file);
        return NULL;
    }

    /* 2. Validate */
    if (new_n < old_n) {
        fprintf(stderr, "Error: new n (%d) < old n (%d)\n", new_n, old_n);
        free(old_q);
        return NULL;
    }
    if (new_d < old_d) {
        fprintf(stderr, "Error: new d (%d) < old d (%d)\n", new_d, old_d);
        free(old_q);
        return NULL;
    }

    /* Check for duplicate prime powers */
    for (int i = 0; i < add_count; i++) {
        for (int j = 0; j < old_k; j++) {
            if (add_q[i] == old_q[j]) {
                fprintf(stderr, "Error: prime power %d already exists\n",
                        add_q[i]);
                free(old_q);
                return NULL;
            }
        }
    }

    /* 3. Merge prime powers */
    int  new_k = old_k + add_count;
    int *merged = malloc((size_t)new_k * sizeof(int));
    if (!merged) { free(old_q); return NULL; }

    memcpy(merged, old_q, (size_t)old_k * sizeof(int));
    memcpy(merged + old_k, add_q, (size_t)add_count * sizeof(int));
    free(old_q);

    /* 4. Allocate and generate */
    CFF *cff = cff_alloc(new_k, merged, new_n, new_d);
    free(merged);
    if (!cff) return NULL;

    cff_generate(cff, nthreads);

    return cff;
}
