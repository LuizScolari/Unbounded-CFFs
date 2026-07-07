#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glib.h>
#include <omp.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fq_nmod.h"
#include "flint/fq_nmod_poly.h"
#include "flint/nmod_poly.h"
#include "cff_builder.h"
#include "cff_file_generator.h"

static fq_nmod_ctx_t global_ctx;
static int global_ctx_initialized = 0;

void rs_generate_cff(char construction, char block_size, int d, long fq, long k);
void embed_cff(char construction, char block_size, const char *cff_file, int d, long Fq, long k);
generated_cffs generate_new_cff_blocks(char construction, char block_size, int d, long* Fq_steps, long* k_steps, int num_steps);

uint64_t** generate_single_cff(long* num_rows, const element_pair* combos, long num_combos, GHashTable* inverted_evals, long num_polys);

guint fq_nmod_hash_func(gconstpointer key);
gboolean fq_nmod_equal_func(gconstpointer a, gconstpointer b);
GHashTable* create_inverted_evaluation_index(long num_points, long num_polys, const fq_nmod_t* points, const fq_nmod_poly_t* polys, const fq_nmod_ctx_t ctx);
void g_array_destroy_wrapper(gpointer data);
void g_hash_table_destroy_wrapper(gpointer data);

combination_partitions generate_combinations(char construction, int dk_size, const subfield_partition* partitions, int num_partitions, const fq_nmod_ctx_t ctx);
void add_pair_to_list(element_pair** list, long* count, long* capacity, const fq_nmod_t x, const fq_nmod_t y, const fq_nmod_ctx_t ctx);

polynomial_partition partition_polynomials(const subfield_partition* partitions, const long* k_steps, int num_steps, const fq_nmod_ctx_t ctx);
fq_nmod_poly_t* generate_polynomials_from_coeffs(long* poly_count, long max_degree, const fq_nmod_t* coeffs, long num_coeffs, const fq_nmod_ctx_t ctx);
void generate_recursive_sorted(fq_nmod_poly_t* poly_list, long* current_index, fq_nmod_poly_t current_poly, long degree, const fq_nmod_t* elements, long num_elements, const fq_nmod_ctx_t ctx);
int fq_nmod_poly_is_in_list(const fq_nmod_poly_t poly, const fq_nmod_poly_t* list, long list_count, const fq_nmod_ctx_t ctx);
void add_poly_to_list(fq_nmod_poly_t** list, long* count, long* capacity, const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx);

void get_element_by_arithmetic(fq_nmod_t result, ulong i, const fq_nmod_ctx_t ctx);
subfield_partition* partition_by_subfields(const long* Fq_steps, int num_steps, const fq_nmod_ctx_t ctx);
long find_element_index(const fq_nmod_t element, const fq_nmod_t* list, long list_count, const fq_nmod_ctx_t ctx);
void add_element_to_list(fq_nmod_t** list, long* count, long* capacity, const fq_nmod_t element, const fq_nmod_ctx_t ctx);

static int is_prime(long n);
static int decompose_prime_power(long q, long* p_out, long* n_out);

void free_matrix(uint64_t** matrix, long rows);
void free_subfield_partitions(subfield_partition* partitions, int num_steps, const fq_nmod_ctx_t ctx);
void free_combination_partitions(combination_partitions* combos, const fq_nmod_ctx_t ctx);
static void free_generated_cffs(generated_cffs* cffs);
static void free_polynomial_partition(polynomial_partition* poly_part, const fq_nmod_ctx_t ctx);

/**
 * Generates an initial CFF from scratch and saves it to a file.
 *
 * @param construction Construction type ('p' for polynomial, 'm' for monotone).
 * @param block_size Defines the size of the CFF rows.
 * @param d CFF parameter d.
 * @param fq Finite field size.
 * @param k Maximum polynomial degree.
 */
void rs_generate_cff(char construction, char block_size, int d, long fq, long k) {
    long t0, n0;

    if(block_size == 'f'){
        t0 = fq*fq;
    } else {
        t0 = (d*k+1) * fq;
    }
    n0 = (long)pow(fq, k + 1);

    char filename0[100];
    snprintf(filename0, sizeof(filename0), "CFFs/%d-CFF(%ld,%ld).txt", d, t0, n0);

    long fq_array[1] = { fq };
    long k_array[1] = { k };
    int num_steps = 1;

    generated_cffs new_blocks = generate_new_cff_blocks(construction, block_size, d, fq_array, k_array, num_steps);

    uint64_t** final_cff = new_blocks.cff_new;
    long final_rows = new_blocks.rows_new;
    long final_cols = new_blocks.cols_new;

    if (final_cff == NULL) {
        printf("Error: Failed to generate initial CFF matrix.\n");
        return;
    }
    write_cff_to_file(filename0, construction, d, fq_array, 1, k_array, 1, final_cff, final_rows, final_cols);

    printf("%d-CFF(%ld,%ld) saved!\n", d, t0, n0);

    new_blocks.cff_old_new = NULL;
    new_blocks.cff_new_old = NULL;

    free_generated_cffs(&new_blocks);
}

/**
 * Embeds an existing CFF into a larger finite field.
 *
 * Reads a CFF from file, generates the new blocks required to grow the field,
 * assembles them together with the original matrix (kept in the upper-left
 * block), and saves the expanded CFF.
 *
 * @param construction Construction type ('p' for polynomial, 'm' for monotone).
 * @param block_size Defines the size of the CFF rows.
 * @param cff_file Path of the CFF file to read.
 * @param d CFF parameter d.
 * @param Fq New finite field size.
 * @param k New maximum polynomial degree.
 */
void embed_cff(char construction, char block_size, const char *cff_file, int d, long Fq, long k){

    long old_rows = 0, old_cols = 0;
    struct cff_parameters* params = read_parameters(cff_file);
    if (params == NULL) {
        printf("Error reading parameters from file %s\n", cff_file);
        return;
    }

    uint64_t** cff_old_old = read_cff_from_file(cff_file, &old_rows, &old_cols);

    int new_fqs_count = params->fqs_count + 1;
    int new_ks_count = params->ks_count + 1;
    long* new_Fq_steps = (long*) malloc(new_fqs_count * sizeof(long));
    long* new_k_steps  = (long*) malloc(new_ks_count * sizeof(long));

    if (new_Fq_steps == NULL || new_k_steps == NULL) {
        printf("Error: Failed to allocate memory for step conversion.\n");
        free(params->Fqs);
        free(params->ks);
        free(params);
        free(new_Fq_steps);
        free(new_k_steps);
        free_matrix(cff_old_old, old_rows);
        return;
    }

    for (int i = 0; i < new_fqs_count; i++) {
        if (i < params->fqs_count) {
            new_Fq_steps[i] = (long)params->Fqs[i];
        } else {
            new_Fq_steps[i] = Fq;
        }
    }

    for (int i = 0; i < new_ks_count; i++) {
        if (i < params->ks_count) {
            new_k_steps[i] = (long)params->ks[i];
        } else {
            new_k_steps[i] = k;
        }
    }

    generated_cffs new_blocks = generate_new_cff_blocks(construction, block_size, d, new_Fq_steps, new_k_steps, new_fqs_count);

    long new_total_rows = old_rows + new_blocks.rows_new_old;

    long width_top = old_cols + new_blocks.cols_old_new;
    long width_bottom = new_blocks.cols_new_old + new_blocks.cols_new;
    long new_total_cols = (width_top > width_bottom) ? width_top : width_bottom;

    long words_per_row = WORDS_FOR_BITS(new_total_cols);
    uint64_t** final_cff = (uint64_t**) malloc(new_total_rows * sizeof(uint64_t*));
    for (long i = 0; i < new_total_rows; i++) {
        final_cff[i] = (uint64_t*) calloc(words_per_row, sizeof(uint64_t));
    }

    for(long i=0; i < old_rows; i++) {
        for (long j = 0; j < old_cols; j++) {
            if (GET_BIT(cff_old_old[i], j)) {
                SET_BIT(final_cff[i], j);
            }
        }
    }

    for(long i=0; i < new_blocks.rows_old_new; i++) {
        for (long j = 0; j < new_blocks.cols_old_new; j++) {
            if (GET_BIT(new_blocks.cff_old_new[i], j)) {
                SET_BIT(final_cff[i], old_cols + j);
            }
        }
    }

    for(long i=0; i < new_blocks.rows_new_old; i++) {
        for (long j = 0; j < new_blocks.cols_new_old; j++) {
            if (GET_BIT(new_blocks.cff_new_old[i], j)) {
                SET_BIT(final_cff[old_rows + i], j);
            }
        }
    }

    for(long i=0; i < new_blocks.rows_new; i++) {
        for (long j = 0; j < new_blocks.cols_new; j++) {
            if (GET_BIT(new_blocks.cff_new[i], j)) {
                SET_BIT(final_cff[old_rows + i], new_blocks.cols_new_old + j);
            }
        }
    }

    long t1 = 0, n1 = 0;

    if (construction == 'm') {
        t1 = (d * k + 1) * Fq;
        n1 = (long)pow(Fq, k + 1);
    } else {
        if(d < ((Fq-1)/k)){
            t1 = (d*k+1)*Fq;
        } else {
            if(block_size == 'f'){
                t1 = (Fq * Fq);
            } else {
                t1 = (d*k+1)*Fq;
            }
        }
        n1 = (long)pow(Fq, k + 1);
    }

    char filename1[100];
    snprintf(filename1, sizeof(filename1), "CFFs/%d-CFF(%ld,%ld).txt", d, t1, n1);

    write_cff_to_file(filename1, construction, d, new_Fq_steps, new_fqs_count, new_k_steps, new_ks_count, final_cff, new_total_rows, new_total_cols);

    printf("%d-CFF(%ld,%ld) -> %d-CFF(%ld,%ld) saved!\n", d, old_rows, old_cols, d, t1, n1);

    free_matrix(cff_old_old, old_rows);
    free_generated_cffs(&new_blocks);
    free_matrix(final_cff, new_total_rows);
    free(new_Fq_steps);
    free(new_k_steps);
    free(params->Fqs);
    free(params->ks);
    free(params);
}

/**
 * Generates the new CFF blocks needed for one embedding step.
 *
 * Orchestrates the generation of the three blocks that grow a CFF: old_new,
 * new_old and new_new, by partitioning field elements and polynomials into
 * "old" (inherited from previous steps) and "new" sets.
 *
 * @param construction Construction type ('p' or 'm').
 * @param block_size Defines the size of the CFF rows.
 * @param d CFF parameter d (used for the monotone construction).
 * @param Fq_steps Array with the finite field sizes of each step.
 * @param k_steps Array with the maximum polynomial degrees of each step.
 * @param num_steps Number of steps.
 * @return Structure containing the three generated blocks.
 */
generated_cffs generate_new_cff_blocks(char construction, char block_size, int d, long* Fq_steps, long* k_steps, int num_steps) {
    generated_cffs result = {0};

    long q_final = Fq_steps[num_steps - 1];
    long p, n;

    if (!decompose_prime_power(q_final, &p, &n)) {
        fprintf(stderr, "Error: %ld is not a prime power!\n", q_final);
        return result;
    }

    fmpz_t pz;
    fmpz_init(pz);
    fmpz_set_ui(pz, (ulong)p);

    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init_ui(ctx, (ulong)p, (slong)n, "a");

    subfield_partition* partitions = partition_by_subfields(Fq_steps, num_steps, ctx);

    polynomial_partition poly_part = partition_polynomials(partitions, k_steps, num_steps, ctx);

    long num_new_rows = 0;
    int dk_size = 0;
    if(construction == 'p'){
        if (num_steps == 1) {
            if(d < ((Fq_steps[0]-1)/k_steps[0])){
                num_new_rows = (d*k_steps[0]+1)*Fq_steps[0];
                dk_size = (d*k_steps[0])+1;
            } else {
                if(block_size == 'f'){
                    num_new_rows = (Fq_steps[0] * Fq_steps[0]);
                    dk_size = Fq_steps[0];
                } else {
                    num_new_rows = (d*k_steps[0]+1)*Fq_steps[0];
                    dk_size = (d*k_steps[0])+1;
                }
            }
        } else {
            if(Fq_steps[num_steps-1] != Fq_steps[num_steps-2]){
                if(d < ((Fq_steps[num_steps-1]-1)/k_steps[num_steps-1])){
                    num_new_rows = ((d*k_steps[num_steps-1]+1)*Fq_steps[num_steps-1]) - (Fq_steps[num_steps-2] * Fq_steps[num_steps-2]);
                    dk_size = d*k_steps[num_steps-1]+1;
                } else {
                    if(block_size == 'f'){
                        num_new_rows = (Fq_steps[num_steps-1] * Fq_steps[num_steps-1]) - (Fq_steps[num_steps-2] * Fq_steps[num_steps-2]);
                        dk_size = Fq_steps[num_steps-1];
                    } else {
                        num_new_rows = ((d*k_steps[num_steps-1]+1)*Fq_steps[num_steps-1]) - (Fq_steps[num_steps-2] * Fq_steps[num_steps-2]);
                        dk_size = d*k_steps[num_steps-1]+1;
                    }
                }
            }
        }
    } else if (construction == 'm') {
        num_new_rows = (d*k_steps[0]+1)*Fq_steps[num_steps-1] - (d*k_steps[num_steps-2]+1)*Fq_steps[num_steps-2];
        dk_size = d*k_steps[0]+1;
    };

    combination_partitions combos = generate_combinations(construction, dk_size, partitions, num_steps, ctx);

    subfield_partition all_partition = partitions[num_steps - 1];
    fq_nmod_t* points_for_eval = all_partition.all_elements;
    long num_points = all_partition.count_all;
    GHashTable* inverted_index_old = create_inverted_evaluation_index(num_points, poly_part.num_old_polys, points_for_eval, poly_part.old_polys, ctx);
    GHashTable* inverted_index_new = create_inverted_evaluation_index(num_points, poly_part.num_new_polys, points_for_eval, poly_part.new_polys, ctx);

    if (!global_ctx_initialized) {
        fq_nmod_ctx_init_modulus(global_ctx, fq_nmod_ctx_modulus(ctx), "a");
        global_ctx_initialized = 1;
    }

    result.cff_old_new = generate_single_cff(&result.rows_old_new, combos.combos_old, combos.count_old, inverted_index_new, poly_part.num_new_polys);
    result.cols_old_new = poly_part.num_new_polys;

    result.cff_new_old = generate_single_cff(&result.rows_new_old, combos.combos_new, num_new_rows, inverted_index_old, poly_part.num_old_polys);
    result.cols_new_old = poly_part.num_old_polys;

    result.cff_new = generate_single_cff(&result.rows_new, combos.combos_new, num_new_rows, inverted_index_new, poly_part.num_new_polys);
    result.cols_new = poly_part.num_new_polys;

    g_hash_table_destroy(inverted_index_old);
    g_hash_table_destroy(inverted_index_new);
    free_combination_partitions(&combos, ctx);
    free_polynomial_partition(&poly_part, ctx);
    free_subfield_partitions(partitions, num_steps, ctx);
    fq_nmod_ctx_clear(ctx);
    fq_nmod_ctx_clear(global_ctx);
    global_ctx_initialized = 0;

    return result;
}

/**
 * Builds a single CFF block as a bit-packed matrix.
 *
 * Each row corresponds to a pair (x, y) and each column to a polynomial.
 * Bit (i, j) is set when polynomial j evaluates to y at point x of pair i,
 * looked up through the inverted evaluation index.
 *
 * @param num_rows Output: number of rows in the block.
 * @param combos Array of element pairs (one per row).
 * @param num_combos Number of pairs.
 * @param inverted_evals Inverted evaluation index.
 * @param num_polys Number of polynomials (columns).
 * @return CFF matrix in bit-packed format, or NULL when there are no rows.
 */
uint64_t** generate_single_cff(long* num_rows, const element_pair* combos, long num_combos, GHashTable* inverted_evals, long num_polys) {
    *num_rows = num_combos;
    if (num_combos == 0) return NULL;

    long words_per_row = WORDS_FOR_BITS(num_polys);

    uint64_t** cff_matrix = (uint64_t**) malloc(num_combos * sizeof(uint64_t*));
    if (cff_matrix == NULL) exit(EXIT_FAILURE);

    #pragma omp parallel for schedule(dynamic)
    for (long i = 0; i < num_combos; i++) {

        const fq_nmod_t* x = &combos[i].x;
        const fq_nmod_t* y = &combos[i].y;

        cff_matrix[i] = (uint64_t*) calloc(words_per_row, sizeof(uint64_t));

        GHashTable* inner_hash = g_hash_table_lookup(inverted_evals, x);

        if (inner_hash != NULL) {
            GArray* indices = g_hash_table_lookup(inner_hash, y);

            if (indices != NULL) {
                for (guint k = 0; k < indices->len; k++) {
                    long j = g_array_index(indices, long, k);
                    SET_BIT(cff_matrix[i], j);
                }
            }
        }
    }

    return cff_matrix;
}

/**
 * Builds an inverted index of polynomial evaluations.
 *
 * For each point x it stores a mapping y -> list of polynomial indices j such
 * that p_j(x) = y, so the matrix filling can look up in average O(1) which
 * polynomials pass through a given pair (x, y).
 *
 * @param num_points Number of evaluation points.
 * @param num_polys Number of polynomials.
 * @param points Array of evaluation points.
 * @param polys Array of polynomials.
 * @param ctx Finite field context.
 * @return Hash table holding the inverted index.
 */
GHashTable* create_inverted_evaluation_index(long num_points, long num_polys, const fq_nmod_t* points, const fq_nmod_poly_t* polys, const fq_nmod_ctx_t ctx) {
    GHashTable* inverted_index = g_hash_table_new_full(fq_nmod_hash_func, fq_nmod_equal_func, g_free, g_hash_table_destroy_wrapper);

    #pragma omp parallel
    {
        fq_nmod_t y_eval_local;
        fq_nmod_init(y_eval_local, ctx);

        fq_nmod_t* x_key_copy;
        fq_nmod_t* y_key_copy;

        #pragma omp for schedule(dynamic)
        for (long i = 0; i < num_points; i++) {
            const fq_nmod_t* x = &points[i];

            GHashTable* inner_hash_local = g_hash_table_new_full(fq_nmod_hash_func, fq_nmod_equal_func, g_free, g_array_destroy_wrapper);

            for (long j = 0; j < num_polys; j++) {
                fq_nmod_poly_evaluate_fq_nmod(y_eval_local, polys[j], *x, ctx);

                GArray* indices = g_hash_table_lookup(inner_hash_local, y_eval_local);
                if (indices == NULL) {
                    indices = g_array_new(FALSE, FALSE, sizeof(long));
                    y_key_copy = (fq_nmod_t*) malloc(sizeof(fq_nmod_t));
                    fq_nmod_init(*y_key_copy, ctx);
                    fq_nmod_set(*y_key_copy, y_eval_local, ctx);
                    g_hash_table_insert(inner_hash_local, y_key_copy, indices);
                }
                g_array_append_val(indices, j);
            }

            x_key_copy = (fq_nmod_t*) malloc(sizeof(fq_nmod_t));
            fq_nmod_init(*x_key_copy, ctx);
            fq_nmod_set(*x_key_copy, *x, ctx);

            #pragma omp critical
            g_hash_table_insert(inverted_index, x_key_copy, inner_hash_local);
        }
        fq_nmod_clear(y_eval_local, ctx);
    }

    return inverted_index;
}

/**
 * Hash function for fq_nmod_t elements, based on the degree-zero coefficient.
 *
 * @param key Pointer to the fq_nmod_t element.
 * @return Computed hash value.
 */
guint fq_nmod_hash_func(gconstpointer key) {
    const fq_nmod_t* element = (const fq_nmod_t*)key;
    return nmod_poly_get_coeff_ui(*element, 0);
}

/**
 * Equality comparison for two fq_nmod_t elements, using the global context.
 *
 * @param a Pointer to the first element.
 * @param b Pointer to the second element.
 * @return TRUE if the elements are equal, FALSE otherwise.
 */
gboolean fq_nmod_equal_func(gconstpointer a, gconstpointer b) {
    const fq_nmod_t* elem1 = (const fq_nmod_t*)a;
    const fq_nmod_t* elem2 = (const fq_nmod_t*)b;
    return fq_nmod_equal(*elem1, *elem2, global_ctx);
}

/**
 * Destroys a GArray stored as a value inside a hash table.
 *
 * @param data Pointer to the GArray to free.
 */
void g_array_destroy_wrapper(gpointer data) {
    g_array_free((GArray*)data, TRUE);
}

/**
 * Destroys an inner GHashTable stored as a value inside another hash table.
 *
 * @param data Pointer to the GHashTable to free.
 */
void g_hash_table_destroy_wrapper(gpointer data) {
    g_hash_table_destroy((GHashTable*)data);
}

/**
 * Builds the "old" and "new" lists of element pairs (x, y).
 *
 * Uses the subfield partitions to enumerate the pairs, classifying them as
 * belonging to previous steps ("old") or to the current step ("new").
 *
 * @param construction Construction type ('p' for polynomial, 'm' for monotone).
 * @param dk_size Size of the dk block (used only for the monotone construction).
 * @param partitions Array of subfield partitions.
 * @param num_partitions Number of partitions.
 * @param ctx Finite field context.
 * @return Structure containing the partitioned combinations.
 */
combination_partitions generate_combinations(char construction, int dk_size, const subfield_partition* partitions, int num_partitions, const fq_nmod_ctx_t ctx) {
    combination_partitions result = {0};

    fq_nmod_t* all_accumulated_elements = NULL;
    long accumulated_count = 0;
    long accumulated_capacity = 0;
    fq_nmod_t* dk_block_elements = NULL;

    if(construction == 'm'){
        dk_block_elements = malloc(dk_size * sizeof(fq_nmod_t));
        int start_index = 0;

        fq_nmod_t* current_only_elements = partitions[0].only_elements;
        for(int i=0; i<dk_size; i++){
            fq_nmod_init(dk_block_elements[i], ctx);
            fq_nmod_set(dk_block_elements[i], current_only_elements[start_index + i], ctx);
        }
    }

    for (int i = 0; i < num_partitions; i++) {
        fq_nmod_t* current_only_elements = partitions[i].only_elements;
        long current_only_count = partitions[i].count_only;

        element_pair** target_list = (i == num_partitions - 1) ? &result.combos_new : &result.combos_old;
        long* target_count = (i == num_partitions - 1) ? &result.count_new : &result.count_old;
        long* target_capacity = (i == num_partitions - 1) ? &result.capacity_new : &result.capacity_old;

        if (i == 0) {
            for (long ix = 0; ix < current_only_count; ix++) {
                for (long iy = 0; iy < current_only_count; iy++) {
                    add_pair_to_list(target_list, target_count, target_capacity, current_only_elements[ix], current_only_elements[iy], ctx);
                }
            }
        } else {
            if(construction == 'p'){
                for (long ix = 0; ix < accumulated_count; ix++) {
                    for (long iy = 0; iy < current_only_count; iy++) {
                        add_pair_to_list(target_list, target_count, target_capacity, all_accumulated_elements[ix], current_only_elements[iy], ctx);
                    }
                }
                for (long j = 0; j < current_only_count; j++) {
                    add_element_to_list(&all_accumulated_elements, &accumulated_count, &accumulated_capacity, current_only_elements[j], ctx);
                }
                for (long ix = 0; ix < current_only_count; ix++) {
                    for (long iy = 0; iy < accumulated_count; iy++) {
                        add_pair_to_list(target_list, target_count, target_capacity, current_only_elements[ix], all_accumulated_elements[iy], ctx);
                    }
                }
            } else if (construction == 'm'){
                for (long ix = 0; ix < dk_size; ix++) {
                    for (long iy = 0; iy < current_only_count; iy++) {
                        add_pair_to_list(target_list, target_count, target_capacity, dk_block_elements[ix], current_only_elements[iy], ctx);
                    }
                }
            }
        }

        if (i == 0) {
             for (long j = 0; j < current_only_count; j++) {
                add_element_to_list(&all_accumulated_elements, &accumulated_count, &accumulated_capacity, current_only_elements[j], ctx);
            }
        }
    }

    for (long i = 0; i < accumulated_count; i++) {
        fq_nmod_clear(all_accumulated_elements[i], ctx);
    }
    free(all_accumulated_elements);

    return result;
}

/**
 * Appends an element pair (x, y) to a dynamic array of pairs.
 *
 * Doubles the array capacity when needed, initialising the newly reserved
 * slots in the finite field context.
 *
 * @param list Pointer to the pair array.
 * @param count Pointer to the pair counter.
 * @param capacity Pointer to the current array capacity.
 * @param x First element of the pair.
 * @param y Second element of the pair.
 * @param ctx Finite field context.
 */
void add_pair_to_list(element_pair** list, long* count, long* capacity, const fq_nmod_t x, const fq_nmod_t y, const fq_nmod_ctx_t ctx) {
    if (*count >= *capacity) {
        *capacity = (*capacity == 0) ? 8 : (*capacity) * 2;
        *list = realloc(*list, (*capacity) * sizeof(element_pair));
        for (long i = *count; i < *capacity; i++) {
            fq_nmod_init((*list)[i].x, ctx);
            fq_nmod_init((*list)[i].y, ctx);
        }
    }
    fq_nmod_set((*list)[*count].x, x, ctx);
    fq_nmod_set((*list)[*count].y, y, ctx);
    (*count)++;
}

/**
 * Splits the polynomials into "old" and "new" sets across subfield steps.
 *
 * Generates the polynomials of each step and keeps only those not already seen
 * in previous steps as "old"; the polynomials exclusive to the last step form
 * the "new" set.
 *
 * @param partitions Array of subfield partitions.
 * @param k_steps Array with the maximum degree of each step.
 * @param num_steps Number of steps.
 * @param ctx Finite field context.
 * @return Structure containing the partitioned polynomials.
 */
polynomial_partition partition_polynomials(const subfield_partition* partitions, const long* k_steps, int num_steps, const fq_nmod_ctx_t ctx) {
    polynomial_partition result = {0};

    fq_nmod_poly_t* old_polys = NULL;
    long num_old_polys_total = 0;

    for (int i = 0; i < num_steps - 1; i++) {
        long num_polys_in_step = 0;
        subfield_partition old_partition = partitions[i];

        fq_nmod_poly_t* old_polys_bySteps = generate_polynomials_from_coeffs(
            &num_polys_in_step, k_steps[i], old_partition.all_elements,
            old_partition.count_all, ctx);

        if (num_polys_in_step > 0) {
            fq_nmod_poly_t* unique_polys_in_step = (fq_nmod_poly_t*) malloc(num_polys_in_step * sizeof(fq_nmod_poly_t));
            if (unique_polys_in_step == NULL) { exit(EXIT_FAILURE); }
            long num_unique_in_step = 0;

            for (long k = 0; k < num_polys_in_step; k++) {
                fq_nmod_poly_init(unique_polys_in_step[k], ctx);
            }

            for (long j = 0; j < num_polys_in_step; j++) {
                if (!fq_nmod_poly_is_in_list(old_polys_bySteps[j], old_polys, num_old_polys_total, ctx)) {
                    fq_nmod_poly_set(unique_polys_in_step[num_unique_in_step], old_polys_bySteps[j], ctx);
                    num_unique_in_step++;
                }
            }

            if (num_unique_in_step > 0) {
                fq_nmod_poly_t* temp = realloc(old_polys, (num_old_polys_total + num_unique_in_step) * sizeof(fq_nmod_poly_t));
                if (temp == NULL) {
                    fprintf(stderr, "Error reallocating memory!\n");
                    for (long k = 0; k < num_polys_in_step; k++) fq_nmod_poly_clear(unique_polys_in_step[k], ctx);
                    free(unique_polys_in_step);
                    exit(EXIT_FAILURE);
                }
                old_polys = temp;

                for (long j = 0; j < num_unique_in_step; j++) {
                    long dest_index = num_old_polys_total + j;
                    fq_nmod_poly_init(old_polys[dest_index], ctx);
                    fq_nmod_poly_set(old_polys[dest_index], unique_polys_in_step[j], ctx);
                }

                num_old_polys_total += num_unique_in_step;
            }

            for (long k = 0; k < num_polys_in_step; k++) {
                fq_nmod_poly_clear(unique_polys_in_step[k], ctx);
            }
            free(unique_polys_in_step);
        }

        for (long j = 0; j < num_polys_in_step; j++) {
            fq_nmod_poly_clear(old_polys_bySteps[j], ctx);
        }
        free(old_polys_bySteps);
    }

    subfield_partition all_partition = partitions[num_steps - 1];
    long num_all_polys = 0;
    fq_nmod_poly_t* all_polys = generate_polynomials_from_coeffs(&num_all_polys, k_steps[num_steps-1], all_partition.all_elements, all_partition.count_all, ctx);

    fq_nmod_poly_t* new_polys = (fq_nmod_poly_t*) malloc(num_all_polys * sizeof(fq_nmod_poly_t));
    long num_new_polys = 0;
    for (long i = 0; i < num_all_polys; i++) {
        if (!fq_nmod_poly_is_in_list(all_polys[i], old_polys, num_old_polys_total, ctx)) {
            fq_nmod_poly_init(new_polys[num_new_polys], ctx);
            fq_nmod_poly_set(new_polys[num_new_polys], all_polys[i], ctx);
            num_new_polys++;
        }
    }

    result.old_polys = old_polys;
    result.num_old_polys = num_old_polys_total;
    result.new_polys = new_polys;
    result.num_new_polys = num_new_polys;
    result.all_polys = all_polys;
    result.num_all_polys = num_all_polys;

    return result;
}

/**
 * Generates every polynomial up to a maximum degree over a set of coefficients.
 *
 * @param poly_count Output: number of generated polynomials.
 * @param max_degree Maximum polynomial degree.
 * @param coeffs Array of admissible coefficients.
 * @param num_coeffs Number of admissible coefficients.
 * @param ctx Finite field context.
 * @return Array of generated polynomials.
 */
fq_nmod_poly_t* generate_polynomials_from_coeffs(long* poly_count, long max_degree, const fq_nmod_t* coeffs, long num_coeffs, const fq_nmod_ctx_t ctx) {
    fmpz_t num_coeffs_z, total_polys_z;
    fmpz_init(num_coeffs_z);
    fmpz_init(total_polys_z);
    fmpz_set_si(num_coeffs_z, num_coeffs);

    fmpz_pow_ui(total_polys_z, num_coeffs_z, max_degree + 1);
    *poly_count = fmpz_get_si(total_polys_z);

    fq_nmod_poly_t* poly_list = (fq_nmod_poly_t*) malloc((*poly_count) * sizeof(fq_nmod_poly_t));
    for(long i=0; i < (*poly_count); i++) {
        fq_nmod_poly_init(poly_list[i], ctx);
    }

    fq_nmod_poly_t temp_poly;
    fq_nmod_poly_init(temp_poly, ctx);
    long start_index = 0;

    generate_recursive_sorted(poly_list, &start_index, temp_poly, max_degree, coeffs, num_coeffs, ctx);

    fq_nmod_poly_clear(temp_poly, ctx);
    fmpz_clear(num_coeffs_z);
    fmpz_clear(total_polys_z);

    return poly_list;
}

/**
 * Recursively builds polynomials, fixing coefficients from the highest degree
 * down to the constant term so the resulting list is in a stable order.
 *
 * @param poly_list Array where the generated polynomials are stored.
 * @param current_index Pointer to the current write index in poly_list.
 * @param current_poly Polynomial being built.
 * @param degree Degree currently being assigned.
 * @param elements Array of admissible coefficients.
 * @param num_elements Number of admissible coefficients.
 * @param ctx Finite field context.
 */
void generate_recursive_sorted(
    fq_nmod_poly_t* poly_list,
    long* current_index,
    fq_nmod_poly_t current_poly,
    long degree,
    const fq_nmod_t* elements,
    long num_elements,
    const fq_nmod_ctx_t ctx)
{
    if (degree < 0) {
        fq_nmod_poly_set(poly_list[*current_index], current_poly, ctx);
        (*current_index)++;
        return;
    }

    for (long i = 0; i < num_elements; i++) {
        fq_nmod_poly_set_coeff(current_poly, degree, elements[i], ctx);
        generate_recursive_sorted(poly_list, current_index, current_poly, degree - 1, elements, num_elements, ctx);
    }
}

/**
 * Checks whether a polynomial is present in an array of polynomials.
 *
 * @param poly Polynomial to search for.
 * @param list Array of polynomials.
 * @param list_count Number of polynomials in the array.
 * @param ctx Finite field context.
 * @return 1 if found, 0 otherwise.
 */
int fq_nmod_poly_is_in_list(const fq_nmod_poly_t poly, const fq_nmod_poly_t* list, long list_count, const fq_nmod_ctx_t ctx) {
    for (long i = 0; i < list_count; i++) {
        if (fq_nmod_poly_equal(poly, list[i], ctx)) {
            return 1;
        }
    }
    return 0;
}

/**
 * Appends a polynomial to a dynamic array, doubling capacity when needed.
 *
 * @param list Pointer to the polynomial array.
 * @param count Pointer to the polynomial counter.
 * @param capacity Pointer to the current array capacity.
 * @param poly Polynomial to append.
 * @param ctx Finite field context.
 */
void add_poly_to_list(fq_nmod_poly_t** list, long* count, long* capacity, const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx) {
    if (*count >= *capacity) {
        *capacity = (*capacity == 0) ? 8 : (*capacity) * 2;
        *list = realloc(*list, (*capacity) * sizeof(fq_nmod_poly_t));
        for (long i = *count; i < *capacity; i++) {
            fq_nmod_poly_init((*list)[i], ctx);
        }
    }
    fq_nmod_poly_set((*list)[*count], poly, ctx);
    (*count)++;
}

/**
 * Builds a finite field element from its arithmetic index.
 *
 * Interprets the index in base p (the field characteristic): each base-p digit
 * is the coefficient of a power of the field generator.
 *
 * @param result Resulting element.
 * @param i Element index (0 to q-1).
 * @param ctx Finite field context.
 */
void get_element_by_arithmetic(fq_nmod_t result, ulong i, const fq_nmod_ctx_t ctx) {
    fq_nmod_zero(result, ctx);
    if (i == 0) return;

    fq_nmod_t gen;
    fq_nmod_init(gen, ctx);
    fq_nmod_gen(gen, ctx);

    fq_nmod_t power_of_a;
    fq_nmod_init(power_of_a, ctx);
    fq_nmod_one(power_of_a, ctx);

    ulong p = fq_nmod_ctx_prime(ctx);
    ulong temp_i = i;

    while (temp_i > 0) {
        ulong remainder = temp_i % p;
        if (remainder != 0) {
            fq_nmod_t term;
            fq_nmod_init(term, ctx);
            fq_nmod_mul_ui(term, power_of_a, remainder, ctx);
            fq_nmod_add(result, result, term, ctx);
            fq_nmod_clear(term, ctx);
        }
        fq_nmod_mul(power_of_a, power_of_a, gen, ctx);
        temp_i /= p;
    }
    fq_nmod_clear(power_of_a, ctx);
    fq_nmod_clear(gen, ctx);
}

/**
 * Partitions the field elements by subfield membership.
 *
 * For each requested subfield F_q it collects the elements satisfying x^q = x,
 * separating those seen for the first time ("only") from all matches ("all").
 *
 * @param Fq_steps Array with the subfield sizes.
 * @param num_steps Number of subfields.
 * @param ctx Finite field context.
 * @return Array of subfield partitions.
 */
subfield_partition* partition_by_subfields(const long* Fq_steps, int num_steps, const fq_nmod_ctx_t ctx) {
    fmpz_t order_z;
    fmpz_init(order_z);
    fq_nmod_ctx_order(order_z, ctx);
    long field_size = fmpz_get_si(order_z);
    fmpz_clear(order_z);

    subfield_partition* partitions = malloc(num_steps * sizeof(subfield_partition));
    char* is_seen = calloc(field_size, sizeof(char));

    fq_nmod_t elem, elem_pow_q;
    fmpz_t q_fmpz;
    fq_nmod_init(elem, ctx);
    fq_nmod_init(elem_pow_q, ctx);
    fmpz_init(q_fmpz);

    for (int i = 0; i < num_steps; i++) {
        long q = Fq_steps[i];

        partitions[i] = (subfield_partition){.q = q, .count_all = 0, .capacity_all = 0, .count_only = 0, .capacity_only = 0};
        fmpz_set_ui(q_fmpz, q);

        for (long j = 0; j < field_size; j++) {
            get_element_by_arithmetic(elem, j, ctx);

            fq_nmod_pow(elem_pow_q, elem, q_fmpz, ctx);

            if (fq_nmod_equal(elem_pow_q, elem, ctx)) {
                add_element_to_list(&partitions[i].all_elements, &partitions[i].count_all, &partitions[i].capacity_all, elem, ctx);

                if (is_seen[j] == 0) {
                    add_element_to_list(&partitions[i].only_elements, &partitions[i].count_only, &partitions[i].capacity_only, elem, ctx);
                    is_seen[j] = 1;
                }
            }
        }
    }

    free(is_seen);
    fq_nmod_clear(elem, ctx);
    fq_nmod_clear(elem_pow_q, ctx);
    fmpz_clear(q_fmpz);

    return partitions;
}

/**
 * Finds the index of an element in an array of elements.
 *
 * @param element Element to search for.
 * @param list Array of elements.
 * @param list_count Number of elements in the array.
 * @param ctx Finite field context.
 * @return Index of the element if found, -1 otherwise.
 */
long find_element_index(const fq_nmod_t element, const fq_nmod_t* list, long list_count, const fq_nmod_ctx_t ctx) {
    for (long i = 0; i < list_count; i++) {
        if (fq_nmod_equal(element, list[i], ctx)) {
            return i;
        }
    }
    return -1;
}

/**
 * Appends a field element to a dynamic array, doubling capacity when needed.
 *
 * @param list Pointer to the element array.
 * @param count Pointer to the element counter.
 * @param capacity Pointer to the current array capacity.
 * @param element Element to append.
 * @param ctx Finite field context.
 */
void add_element_to_list(fq_nmod_t** list, long* count, long* capacity, const fq_nmod_t element, const fq_nmod_ctx_t ctx) {
    if (*count >= *capacity) {
        *capacity = (*capacity == 0) ? 8 : (*capacity) * 2;
        *list = realloc(*list, (*capacity) * sizeof(fq_nmod_t));
        for (long i = *count; i < *capacity; i++) {
            fq_nmod_init((*list)[i], ctx);
        }
    }
    fq_nmod_set((*list)[*count], element, ctx);
    (*count)++;
}

/**
 * Checks whether a number is prime.
 *
 * @param n Number to test.
 * @return 1 if n is prime, 0 otherwise.
 */
static int is_prime(long n) {
    if (n <= 1) return 0;
    if (n <= 3) return 1;
    if (n % 2 == 0 || n % 3 == 0) return 0;
    for (long i = 5; i * i <= n; i += 6) {
        if (n % i == 0 || n % (i + 2) == 0) return 0;
    }
    return 1;
}

/**
 * Decomposes q into p^n with p prime.
 *
 * @param q Finite field size (must be a prime power).
 * @param p_out Output: the prime characteristic p.
 * @param n_out Output: the exponent n.
 * @return 1 on success, 0 if q is not a prime power.
 */
static int decompose_prime_power(long q, long* p_out, long* n_out) {
    if (q <= 1) return 0;

    if (is_prime(q)) {
        *p_out = q;
        *n_out = 1;
        return 1;
    }

    for (long p = 2; p * p <= q; p++) {
        if (!is_prime(p)) continue;

        long temp = q;
        long n = 0;

        while (temp % p == 0) {
            temp /= p;
            n++;
        }

        if (temp == 1 && n > 0) {
            *p_out = p;
            *n_out = n;
            return 1;
        }
    }

    return 0;
}

/**
 * Frees a bit-packed uint64_t matrix.
 *
 * @param matrix Matrix to free.
 * @param rows Number of rows in the matrix.
 */
void free_matrix(uint64_t** matrix, long rows) {
    if (!matrix) return;
    for (long i = 0; i < rows; i++) free(matrix[i]);
    free(matrix);
}

/**
 * Frees an array of subfield partitions.
 *
 * @param partitions Array of partitions to free.
 * @param num_steps Number of partitions in the array.
 * @param ctx Finite field context.
 */
void free_subfield_partitions(subfield_partition* partitions, int num_steps, const fq_nmod_ctx_t ctx) {
    if (!partitions) return;
    for (int i = 0; i < num_steps; i++) {
        for (long j = 0; j < partitions[i].capacity_all; j++) fq_nmod_clear(partitions[i].all_elements[j], ctx);
        for (long j = 0; j < partitions[i].capacity_only; j++) fq_nmod_clear(partitions[i].only_elements[j], ctx);
        free(partitions[i].all_elements);
        free(partitions[i].only_elements);
    }
    free(partitions);
}

/**
 * Frees a combination-partition structure.
 *
 * @param combos Pointer to the structure to free.
 * @param ctx Finite field context.
 */
void free_combination_partitions(combination_partitions* combos, const fq_nmod_ctx_t ctx) {
    if (!combos) return;
    for (long i = 0; i < combos->capacity_old; i++) { fq_nmod_clear(combos->combos_old[i].x, ctx); fq_nmod_clear(combos->combos_old[i].y, ctx); }
    free(combos->combos_old);
    for (long i = 0; i < combos->capacity_new; i++) { fq_nmod_clear(combos->combos_new[i].x, ctx); fq_nmod_clear(combos->combos_new[i].y, ctx); }
    free(combos->combos_new);
}

/**
 * Frees a generated-CFF-blocks structure.
 *
 * @param cffs Pointer to the structure to free.
 */
static void free_generated_cffs(generated_cffs* cffs) {
    if (!cffs) return;
    if (cffs->cff_old_new) { for (long i = 0; i < cffs->rows_old_new; i++) free(cffs->cff_old_new[i]); free(cffs->cff_old_new); }
    if (cffs->cff_new_old) { for (long i = 0; i < cffs->rows_new_old; i++) free(cffs->cff_new_old[i]); free(cffs->cff_new_old); }
    if (cffs->cff_new) { for (long i = 0; i < cffs->rows_new; i++) free(cffs->cff_new[i]); free(cffs->cff_new); }
}

/**
 * Frees a polynomial-partition structure.
 *
 * @param poly_part Pointer to the structure to free.
 * @param ctx Finite field context.
 */
static void free_polynomial_partition(polynomial_partition* poly_part, const fq_nmod_ctx_t ctx) {
    if (!poly_part) return;
    if (poly_part->old_polys != NULL) {
        for (long i = 0; i < poly_part->num_old_polys; i++) fq_nmod_poly_clear(poly_part->old_polys[i], ctx);
        free(poly_part->old_polys);
    }
    if (poly_part->new_polys != NULL) {
        for (long i = 0; i < poly_part->num_new_polys; i++) fq_nmod_poly_clear(poly_part->new_polys[i], ctx);
        free(poly_part->new_polys);
    }
    if (poly_part->all_polys != NULL) {
        for (long i = 0; i < poly_part->num_all_polys; i++) fq_nmod_poly_clear(poly_part->all_polys[i], ctx);
        free(poly_part->all_polys);
    }
}
