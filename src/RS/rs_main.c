#include <stdio.h>
#include <stdlib.h>
#include "cff_builder.h"
#include <sys/stat.h>

/**
 * Entry point of the RS construction (dispatched from src/main.c).
 *
 * Processes command-line arguments to generate an initial CFF or to embed an
 * existing one. argv[0] is the construction name ("RS").
 *
 * Usage:
 *   - Polynomial first: ./cff RS p f <m|f> <d> <q> <k>
 *   - Embedding CFFs:   ./cff RS p g <m|f> <cff_file> <d> <q> <k>
 *   - Monotone CFFs:    ./cff RS m g <cff_file> <d> <q> <k>
 *
 * @param argc Number of arguments.
 * @param argv Array of arguments.
 * @return 0 on success, 1 on error.
 */
int rs_main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Error: Insufficient arguments.\n");
        return 1;
    }

    char construction = argv[1][0];
    char action = argv[2][0];
    char block_size = '\0';

    if (construction == 'p') {
        if (argc < 4) {
            fprintf(stderr, "Error: Missing block size for polynomial construction.\n");
            return 1;
        }
        block_size = argv[3][0];
    }

    mkdir("CFFs", 0777);

    if (action == 'g') {
        if (construction == 'p') {
            if (argc != 8) {
                fprintf(stderr, "Error (p g): Incorrect number of arguments.\n");
                fprintf(stderr, "Usage: ./cff RS p g <m|f> <cff_file> <d> <q> <k>\n");
                return 1;
            }
            if (block_size != 'm' && block_size != 'f') {
                fprintf(stderr, "Error: Block size must be 'm' (minimum) or 'f' (full).\n");
                return 1;
            }
        }

        if (construction == 'm') {
            if (argc != 7) {
                fprintf(stderr, "Error (m g): Incorrect number of arguments.\n");
                fprintf(stderr, "Usage: ./cff RS m g <cff_file> <d> <q> <k>\n");
                return 1;
            }
        }

        int d = 0;
        long Fq = 0, k = 0;
        char *cff_file = NULL;
        if (construction == 'p') {
            cff_file = argv[4];
            d = atoi(argv[5]);
            Fq = atol(argv[6]);
            k = atol(argv[7]);
        } else if (construction == 'm') {
            cff_file = argv[3];
            d = atoi(argv[4]);
            Fq = atol(argv[5]);
            k = atol(argv[6]);
        }

        embed_cff(construction, block_size, cff_file, d, Fq, k);

    } else if (action == 'f') {
        if (construction != 'p') {
            fprintf(stderr, "Error: Action 'f' is only available for polynomial construction ('p').\n");
            fprintf(stderr, "For monotone construction ('m'), use only action 'g' (embedding).\n");
            return 1;
        }

        if (argc != 7) {
            fprintf(stderr, "Error (p f): Incorrect number of arguments.\n");
            fprintf(stderr, "Usage: ./cff RS p f <m|f> <d> <q> <k>\n");
            return 1;
        }

        if (block_size != 'm' && block_size != 'f') {
            fprintf(stderr, "Error: Block size must be 'm' (minimum) or 'f' (full).\n");
            return 1;
        }

        int d = atoi(argv[4]);
        long Fq = atol(argv[5]);
        long k = atol(argv[6]);

        rs_generate_cff(construction, block_size, d, Fq, k);

    } else {
        fprintf(stderr, "Error: Unknown action '%c'. Use 'g' or 'f'.\n", action);
        return 1;
    }

    return 0;
}
