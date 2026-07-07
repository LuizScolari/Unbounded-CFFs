#include <stdio.h>
#include <strings.h>

/* Entry points implemented in src/RS, src/Subset and src/CRS. */
int rs_main(int argc, char *argv[]);
int subset_main(int argc, char *argv[]);
int crs_main(int argc, char *argv[]);

/**
 * Prints the top-level usage message listing the available constructions.
 *
 * @param prog Program name to embed in the usage text (argv[0]).
 */
static void usage(const char *prog)
{
    fprintf(stderr,
        "Unbounded Cover-Free Families\n\n"
        "Usage: %s <construction> <args...>\n\n"
        "Constructions:\n"
        "  RS      Polynomial (Reed-Solomon) construction over finite fields.\n"
        "          %s RS <p|m> <f|g> ...\n"
        "  Subset  Subset-inclusion construction.\n"
        "          %s Subset <generate|embed> ...\n"
        "  CRS     Chinese Remainder Sieve construction.\n"
        "          %s CRS <generate|embed> ...\n\n"
        "Run a construction without further arguments to see its own help.\n",
        prog, prog, prog, prog);
}

/**
 * Selects a construction from the first argument and forwards the rest to it.
 *
 * The chosen entry point receives argv shifted by one, so its argv[0] is the
 * construction name (e.g. "RS") and argv[1] is its first sub-argument.
 *
 * @param argc Number of command-line arguments.
 * @param argv Command-line arguments; argv[1] selects the construction.
 * @return 0 on success, 1 on error or unknown construction.
 */
int main(int argc, char *argv[])
{
    if (argc < 2) {
        usage(argv[0]);
        return 1;
    }

    const char *construction = argv[1];

    /* Forward argv shifted by one: the entry point sees the construction
       name as its argv[0] and its own sub-arguments from argv[1] on. */
    if (strcasecmp(construction, "RS") == 0)
        return rs_main(argc - 1, argv + 1);
    else if (strcasecmp(construction, "Subset") == 0)
        return subset_main(argc - 1, argv + 1);
    else if (strcasecmp(construction, "CRS") == 0)
        return crs_main(argc - 1, argv + 1);

    fprintf(stderr, "Unknown construction: '%s'\n\n", construction);
    usage(argv[0]);
    return 1;
}
