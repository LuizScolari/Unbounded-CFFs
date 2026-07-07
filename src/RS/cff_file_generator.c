#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "cff_file_generator.h"
#include "cff_builder.h"

int* parse_int_list(char* str, int* count);

/**
 * Reads the CFF parameters from a file's first line.
 *
 * Extracts the construction type, the d parameter (for the monotone
 * construction) and the arrays of field sizes (Fqs) and degrees (ks).
 *
 * @param filename Path of the file to read.
 * @return Pointer to the parameter structure, or NULL on error.
 */
struct cff_parameters* read_parameters(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("Input file '%s' not found.\n", filename);
        return NULL;
    }

    struct cff_parameters* params = (struct cff_parameters*) malloc(sizeof(struct cff_parameters));
    if (params == NULL) {
        printf("Error: Failed to allocate memory for parameters.\n");
        fclose(file);
        return NULL;
    }

    params->Fqs = NULL;
    params->ks = NULL;
    params->fqs_count = 0;
    params->ks_count = 0;
    params->d = 0;

    char* line = NULL;
    size_t len = 0;

    if (getline(&line, &len, file) == -1) {
        printf("Error: File is empty or could not read the line.\n");
        free(params);
        fclose(file);
        if (line) free(line);
        return NULL;
    }

    char fqs_str[256];
    char ks_str[256];
    int items_scanned;

    if (sscanf(line, "%c", &params->construction) != 1) {
        printf("Error: Could not read construction type.\n");
        free(params);
        fclose(file);
        free(line);
        return NULL;
    }

    if (params->construction == 'm') {
        items_scanned = sscanf(line, "%c %d [%[^]]] [%[^]]]",
                               &params->construction, &params->d, fqs_str, ks_str);
        if (items_scanned != 4) {
            printf("Error: Invalid first line format for 'm' construction.\n");
            free(params);
            fclose(file);
            free(line);
            return NULL;
        }
    } else {
        items_scanned = sscanf(line, "%c [%[^]]] [%[^]]]",
                               &params->construction, fqs_str, ks_str);
        if (items_scanned != 3) {
            printf("Error: Invalid first line format.\n");
            free(params);
            fclose(file);
            free(line);
            return NULL;
        }
    }

    params->Fqs = parse_int_list(fqs_str, &params->fqs_count);
    params->ks  = parse_int_list(ks_str,  &params->ks_count);

    fclose(file);
    free(line);

    if ((params->fqs_count > 0 && params->Fqs == NULL) ||
        (params->ks_count > 0 && params->ks == NULL)) {
        printf("Error: Failed to allocate memory for Fqs/ks lists.\n");
        free(params->Fqs);
        free(params->ks);
        free(params);
        return NULL;
    }

    return params;
}

/**
 * Reads a CFF matrix from a file into bit-packed format.
 *
 * Skips the parameter line, then reads the binary matrix, storing each row as
 * 64-bit words for efficiency.
 *
 * @param filename Path of the file to read.
 * @param rows Output: number of rows.
 * @param cols Output: number of columns.
 * @return CFF matrix in bit-packed format, or NULL if the file cannot be read.
 */
uint64_t** read_cff_from_file(const char* filename, long* rows, long* cols) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("Input file '%s' not found. Starting from scratch.\n", filename);
        *rows = 0; *cols = 0;
        return NULL;
    }

    *rows = 0; *cols = 0;
    char* line = NULL;
    size_t len = 0;

    if (getline(&line, &len, file) == -1) {
        fclose(file);
        if(line) free(line);
        printf("File '%s' is empty or has no data after the first line.\n", filename);
        return NULL;
    }

    if (getline(&line, &len, file) != -1) {
        (*rows)++;
        char* token = strtok(line, " \t\n");
        while (token != NULL) {
            if (strlen(token) > 0) (*cols)++;
            token = strtok(NULL, " \t\n");
        }
    }
    while (getline(&line, &len, file) != -1) {
        (*rows)++;
    }

    rewind(file);
    if (getline(&line, &len, file) == -1) {
        fclose(file);
        if (line) free(line);
        return NULL;
    }

    long words_per_row = WORDS_FOR_BITS(*cols);
    uint64_t** matrix = (uint64_t**) malloc(*rows * sizeof(uint64_t*));

    for (long i = 0; i < *rows; i++) {
        matrix[i] = (uint64_t*) calloc(words_per_row, sizeof(uint64_t));
        for (long j = 0; j < *cols; j++) {
            int bit;
            if (fscanf(file, "%d", &bit) == 1 && bit == 1) {
                SET_BIT(matrix[i], j);
            }
        }
    }

    fclose(file);
    if(line) free(line);
    return matrix;
}

/**
 * Writes a CFF matrix to a file.
 *
 * The first line holds the parameters and each subsequent line is one binary
 * matrix row.
 *
 * @param filename Output file path.
 * @param construction Construction type ('p' or 'm').
 * @param d CFF parameter d (used for the monotone construction).
 * @param Fq_steps Array with the finite field sizes.
 * @param fqs_count Number of elements in Fq_steps.
 * @param K_steps Array with the maximum polynomial degrees.
 * @param ks_count Number of elements in K_steps.
 * @param matrix CFF matrix in bit-packed format.
 * @param rows Number of rows in the matrix.
 * @param cols Number of columns in the matrix.
 */
void write_cff_to_file(const char* filename, char construction, int d, long* Fq_steps, int fqs_count, long* K_steps, int ks_count, uint64_t** matrix, long rows, long cols){
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error opening file '%s' for writing.\n", filename);
        return;
    }

    fprintf(file, "%c ", construction);

    if (construction == 'm') {
        fprintf(file, "%d ", d);
    }

    fprintf(file, "[");
    for (int i = 0; i < fqs_count; i++) {
        fprintf(file, "%ld", Fq_steps[i]);
        if (i < fqs_count - 1) {
            fprintf(file, ",");
        }
    }
    fprintf(file, "] ");

    fprintf(file, "[");
    for (int i = 0; i < ks_count; i++) {
        fprintf(file, "%ld", K_steps[i]);
        if (i < ks_count - 1) {
            fprintf(file, ",");
        }
    }
    fprintf(file, "]\n");

    for (long i = 0; i < rows; i++) {
        for (long j = 0; j < cols; j++) {
            fprintf(file, "%d", GET_BIT(matrix[i], j) ? 1 : 0);
            if (j < cols - 1) {
                fprintf(file, " ");
            }
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

/**
 * Parses a comma-separated list of integers into a newly allocated array.
 *
 * @param str String with comma-separated integers (e.g. "2,4,16").
 * @param count Output: number of integers parsed.
 * @return Dynamically allocated integer array, or NULL if the string is empty.
 */
int* parse_int_list(char* str, int* count) {
    *count = 0;

    if (str == NULL || strlen(str) == 0) {
        return NULL;
    }

    int commas = 0;
    for (int i = 0; str[i] != '\0'; i++) {
        if (str[i] == ',') {
            commas++;
        }
    }
    *count = commas + 1;

    int* list = (int*) malloc(*count * sizeof(int));
    if (list == NULL) {
        *count = 0;
        return NULL;
    }

    int i = 0;
    char* token = strtok(str, ",");

    while (token != NULL) {
        list[i] = atoi(token);
        i++;
        token = strtok(NULL, ",");
    }

    return list;
}
