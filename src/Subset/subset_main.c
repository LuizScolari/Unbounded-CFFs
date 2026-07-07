#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>

#include "generate_cff.h"
#include "cff_io.h"

#define OUT_DIR  "CFFs"

/**
 * Prints the usage message for the Subset construction.
 *
 * @param prog Construction name to embed in the usage text (argv[0]).
 */
static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage:\n"
        "  ./cff %s generate <d> <k> <S>\n"
        "      Build M(d,k,S): a d-CFF(C(S,d), C(S,k)).\n\n"
        "  ./cff %s embed <old_file> <k'> <S'>\n"
        "      Embed an existing CFF into a larger M(d,k',S') so every\n"
        "      predecessor occupies its upper-left block.\n\n",
        prog, prog);
}

/**
 * Builds the output path CFFs/d-CFF(t,n).txt for a CFF.
 *
 * @param cff CFF whose dimensions name the file.
 * @param buf Output buffer for the path.
 * @param bufsz Size of buf.
 */
static void make_outpath(const CFF *cff, char *buf, size_t bufsz)
{
    snprintf(buf, bufsz, "%s/%d-CFF(%lld,%lld).txt",
             OUT_DIR, cff->d, cff->t, cff->n);
}

/**
 * Creates the output directory OUT_DIR if it does not exist yet.
 */
static void ensure_out_dir(void)
{
    if (mkdir(OUT_DIR, 0755) != 0 && errno != EEXIST)
        fprintf(stderr, "Warning: could not create %s/\n", OUT_DIR);
}

/**
 * Generates M(d,k,S) (fresh or embedded), saves it to CFFs/ and prints a
 * one-line status: "d-CFF(t,n) saved!" for a fresh matrix, or
 * "d-CFF(t,n) -> d-CFF(t,n) saved!" when it grows from a predecessor.
 *
 * @param d Strength.
 * @param k Column subset size.
 * @param S Ground-set size.
 * @param chain_len Number of predecessors in the embedding chain.
 * @param chain_k Predecessor k values (oldest first), or NULL.
 * @param chain_S Predecessor S values (oldest first), or NULL.
 * @return 0 on success, 1 on error.
 */
static int build_and_save(int d, int k, int S,
                          int chain_len, const int *chain_k, const int *chain_S)
{
    CFF *cff = generate_cff(d, k, S, chain_len, chain_k, chain_S);
    if (!cff) { fputs("Error: allocation failed.\n", stderr); return 1; }

    ensure_out_dir();
    char outpath[512];
    make_outpath(cff, outpath, sizeof(outpath));

    if (write_cff(cff, outpath) != 0) {
        fputs("Error: write failed.\n", stderr);
        free_cff(cff);
        return 1;
    }

    if (chain_len > 0) {
        int Sp = chain_S[chain_len - 1];
        int kp = chain_k[chain_len - 1];
        printf("%d-CFF(%lld,%lld) -> %d-CFF(%lld,%lld) saved!\n",
               d, binom(Sp, d), binom(Sp, kp), d, cff->t, cff->n);
    } else {
        printf("%d-CFF(%lld,%lld) saved!\n", d, cff->t, cff->n);
    }

    free_cff(cff);
    return 0;
}

/**
 * Handles the "generate" command: builds a fresh M(d,k,S).
 *
 * Validates that the parameters satisfy 0 < d < k < S before building.
 *
 * @param argc Number of sub-arguments (after "generate").
 * @param argv Sub-arguments: <d> <k> <S>.
 * @return 0 on success, 1 on error.
 */
static int cmd_generate(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr, "generate: need <d> <k> <S>\n");
        return 1;
    }

    int d = atoi(argv[0]);
    int k = atoi(argv[1]);
    int S = atoi(argv[2]);

    if (d <= 0 || k <= d || S <= k) {
        fprintf(stderr,
                "Error: parameters must satisfy 0 < d < k < S. "
                "Got d=%d, k=%d, S=%d.\n", d, k, S);
        return 1;
    }

    return build_and_save(d, k, S, 0, NULL, NULL);
}

/**
 * Handles the "embed" command: grows an existing CFF into M(d,k',S').
 *
 * Reads M(d,k,S) from file (the full chain is stored in the file), requires
 * k' >= k and S' > S, appends (k,S) to the chain, and generates M(d,k',S') so
 * that every predecessor occupies its upper-left block.
 *
 * @param argc Number of sub-arguments (after "embed").
 * @param argv Sub-arguments: <old_file> <k'> <S'>.
 * @return 0 on success, 1 on error.
 */
static int cmd_embed(int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr, "embed: need <old_file> <k'> <S'>\n");
        return 1;
    }

    const char *infile = argv[0];
    int k_new = atoi(argv[1]);
    int S_new = atoi(argv[2]);

    CFF *base = read_cff(infile);
    if (!base) return 1;

    if (k_new < base->k) {
        fprintf(stderr, "Error: k' (%d) must be >= k (%d)\n", k_new, base->k);
        free_cff(base); return 1;
    }
    if (S_new <= base->S) {
        fprintf(stderr, "Error: S' (%d) must be > S (%d)\n", S_new, base->S);
        free_cff(base); return 1;
    }
    if (k_new >= S_new) {
        fprintf(stderr, "Error: need k' (%d) < S' (%d)\n", k_new, S_new);
        free_cff(base); return 1;
    }

    /* Build new chain: base's chain + base's (k,S) */
    int d = base->d;
    int new_chain_len = base->chain_len + 1;
    if (new_chain_len > MAX_CHAIN) {
        fprintf(stderr, "Error: chain depth exceeds MAX_CHAIN (%d)\n", MAX_CHAIN);
        free_cff(base); return 1;
    }
    int new_chain_k[MAX_CHAIN], new_chain_S[MAX_CHAIN];
    memcpy(new_chain_k, base->chain_k, (size_t)base->chain_len * sizeof(int));
    memcpy(new_chain_S, base->chain_S, (size_t)base->chain_len * sizeof(int));
    new_chain_k[base->chain_len] = base->k;
    new_chain_S[base->chain_len] = base->S;
    free_cff(base);

    return build_and_save(d, k_new, S_new, new_chain_len, new_chain_k, new_chain_S);
}

/**
 * Entry point of the Subset construction (dispatched from src/main.c).
 *
 * Routes to the generate or embed command. argv[0] is the construction name
 * ("Subset").
 *
 * @param argc Number of arguments.
 * @param argv Array of arguments; argv[1] selects the command.
 * @return 0 on success, 1 on error or unknown command.
 */
int subset_main(int argc, char **argv)
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
