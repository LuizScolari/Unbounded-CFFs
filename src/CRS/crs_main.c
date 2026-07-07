#include "cff.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>

/**
 * Returns the number of online CPUs, or 4 if it cannot be determined.
 *
 * @return Thread count to use for parallel generation.
 */
static int get_nthreads(void)
{
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return (n > 0) ? (int)n : 4;
}

/**
 * Prints the usage message for the CRS construction.
 *
 * @param prog Construction name to embed in the usage text (argv[0]).
 */
static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage:\n"
        "  ./cff %s generate <d> <n> <q1> <q2> ... <qk>\n"
        "      Generate a d-CFF(t,n) from prime powers q1..qk.\n\n"
        "  ./cff %s embed <old_file> <d'> <n'> <q_new1> <q_new2> ...\n"
        "      Extend existing CFF with new prime powers / larger n.\n\n",
        prog, prog);
}

/**
 * Handles the "generate" command: builds a d-CFF(t,n) from prime powers.
 *
 * Parses d, n and the prime powers, generates the matrix in parallel, saves it
 * to CFFs/ and prints "d-CFF(t,n) saved!".
 *
 * @param argc Number of sub-arguments (after "generate").
 * @param argv Sub-arguments: <d> <n> <q1> ... <qk>.
 * @return 0 on success, 1 on error.
 */
static int cmd_generate(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr, "generate: need at least d, n, and one prime power\n");
        return 1;
    }

    int d = atoi(argv[0]);
    int n = atoi(argv[1]);
    int k = argc - 2;

    int *q = malloc((size_t)k * sizeof(int));
    for (int i = 0; i < k; i++)
        q[i] = atoi(argv[i + 2]);

    CFF *cff = cff_alloc(k, q, n, d);
    free(q);
    if (!cff) return 1;

    cff_generate(cff, get_nthreads());

    mkdir("CFFs", 0755);
    char filename[256];
    snprintf(filename, sizeof(filename), "CFFs/%d-CFF(%d,%d).txt",
             cff->d, cff->t, cff->n);

    if (cff_save(cff, filename) != 0) {
        cff_free(cff);
        return 1;
    }

    printf("%d-CFF(%d,%d) saved!\n", cff->d, cff->t, cff->n);

    cff_free(cff);
    return 0;
}

/**
 * Handles the "embed" command: extends an existing CFF.
 *
 * Reads the old file's header, adds new prime powers and/or a larger n, and
 * regenerates the extended matrix (the old one stays in the upper-left block).
 * Saves it to CFFs/ and prints "d-CFF(t,n) -> d-CFF(t,n) saved!".
 *
 * @param argc Number of sub-arguments (after "embed").
 * @param argv Sub-arguments: <old_file> <d'> <n'> <q_new1> ...
 * @return 0 on success, 1 on error.
 */
static int cmd_embed(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr, "embed: need <old_file> <d'> <n'> <q_new> ...\n");
        return 1;
    }

    const char *old_file = argv[0];
    int new_d = atoi(argv[1]);
    int new_n = atoi(argv[2]);
    int add_k = argc - 3;

    int *add_q = NULL;
    if (add_k > 0) {
        add_q = malloc((size_t)add_k * sizeof(int));
        for (int i = 0; i < add_k; i++)
            add_q[i] = atoi(argv[i + 3]);
    }

    CFF *cff = cff_embed(old_file, add_k, add_q, new_n, new_d, get_nthreads());
    free(add_q);
    if (!cff) return 1;

    mkdir("CFFs", 0755);
    char filename[256];
    snprintf(filename, sizeof(filename), "CFFs/%d-CFF(%d,%d).txt",
             cff->d, cff->t, cff->n);

    if (cff_save(cff, filename) != 0) {
        cff_free(cff);
        return 1;
    }

    /* Old dimensions (t = sum of prime powers) for the status line. */
    int old_k, old_n, old_d, *old_q = NULL, old_t = 0;
    if (cff_load_header(old_file, &old_k, &old_q, &old_n, &old_d) == 0) {
        for (int j = 0; j < old_k; j++) old_t += old_q[j];
        printf("%d-CFF(%d,%d) -> %d-CFF(%d,%d) saved!\n",
               old_d, old_t, old_n, cff->d, cff->t, cff->n);
        free(old_q);
    }

    cff_free(cff);
    return 0;
}

/**
 * Entry point of the CRS construction (dispatched from src/main.c).
 *
 * Routes to the generate or embed command. argv[0] is the construction
 * name ("CRS").
 *
 * @param argc Number of arguments.
 * @param argv Array of arguments; argv[1] selects the command.
 * @return 0 on success, 1 on error or unknown command.
 */
int crs_main(int argc, char **argv)
{
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    const char *cmd = argv[1];

    if (strcmp(cmd, "generate") == 0)
        return cmd_generate(argc - 2, argv + 2);
    else if (strcmp(cmd, "embed") == 0)
        return cmd_embed(argc - 2, argv + 2);
    else {
        fprintf(stderr, "Unknown command: %s\n", cmd);
        usage(argv[0]);
        return 1;
    }
}