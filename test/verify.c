/**
 * @file verify.c
 * @brief CFF (Cover-Free Family) verification program.
 * 
 * This program reads a binary matrix from a file and verifies
 * whether it satisfies the d-CFF property.
 * 
 * run: gcc -O3 -o prog verify.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

//Update accordingly to the CFF parameters. 
#define MAX_ROWS 1000
#define MAX_COLS 1000
#define MAX_LINE_LENGTH 2000

typedef struct {
    int *indices;
    int size;
    int capacity;
} Block;

typedef struct {
    int rows;
    int cols;
    int **data;
} Matrix;

Matrix* read_matrix_from_file(const char *file_path) {
    FILE *file = fopen(file_path, "r");
    if (!file) {
        perror("Error opening file");
        return NULL;
    }
    
    Matrix *matrix = (Matrix*)malloc(sizeof(Matrix));
    matrix->rows = 0;
    matrix->cols = 0;
    matrix->data = (int**)malloc(MAX_ROWS * sizeof(int*));
    
    char line[MAX_LINE_LENGTH];
    
    if (fgets(line, sizeof(line), file)) {
        char *token = strtok(line, " \t\n");
        while (token) {
            matrix->cols++;
            token = strtok(NULL, " \t\n");
        }
        rewind(file);
    }
    
    while (fgets(line, sizeof(line), file) && matrix->rows < MAX_ROWS) {
        matrix->data[matrix->rows] = (int*)calloc(matrix->cols, sizeof(int));
        
        char *token = strtok(line, " \t\n");
        int col = 0;
        while (token && col < matrix->cols) {
            matrix->data[matrix->rows][col] = atoi(token);
            token = strtok(NULL, " \t\n");
            col++;
        }
        matrix->rows++;
    }
    
    fclose(file);
    return matrix;
}

Block* process_columns(Matrix *matrix, int *num_blocks) {
    *num_blocks = matrix->cols;
    Block *blocks = (Block*)malloc(matrix->cols * sizeof(Block));
    
    for (int col = 0; col < matrix->cols; col++) {
        blocks[col].capacity = matrix->rows;
        blocks[col].indices = (int*)malloc(matrix->rows * sizeof(int));
        blocks[col].size = 0;
        
        for (int row = 0; row < matrix->rows; row++) {
            if (matrix->data[row][col] == 1) {
                blocks[col].indices[blocks[col].size++] = row + 1;
            }
        }
    }
    
    return blocks;
}

void compute_union_fast(Block *blocks, int *selected, int selected_count, bool *seen) {
    for (int i = 0; i < selected_count; i++) {
        Block *block = &blocks[selected[i]];
        for (int j = 0; j < block->size; j++) {
            seen[block->indices[j]] = true;
        }
    }
}

bool is_subset_fast(Block *block1, bool *seen) {
    for (int i = 0; i < block1->size; i++) {
        int idx = block1->indices[i];
        if (!seen[idx]) {
            return false;
        }
    }
    return true;
}

bool check_combinations(Block *blocks, int n, int d, int i, int start, int *selected, int selected_count, bool *seen_array) {
    if (selected_count == d) {
        memset(seen_array, 0, (MAX_ROWS + 1) * sizeof(bool));
        compute_union_fast(blocks, selected, selected_count, seen_array);
        if (is_subset_fast(&blocks[i], seen_array)) {
            return true;
        }
        return false;
    }
    
    for (int j = start; j < n; j++) {
        if (j == i) continue;
        
        selected[selected_count] = j;
        if (check_combinations(blocks, n, d, i, j + 1, selected, selected_count + 1, seen_array)) {
            return true;
        }
    }
    
    return false;
}

bool is_cff(Block *blocks, int n, int d) {
    int *selected = (int*)malloc(n * sizeof(int));
    bool *seen_array = (bool*)calloc(MAX_ROWS + 1, sizeof(bool));
    
    for (int i = 0; i < n; i++) {
        if (check_combinations(blocks, n, d, i, 0, selected, 0, seen_array)) {
            printf("No\n");
            free(selected);
            free(seen_array);
            return false;
        }
    }
    
    printf("Yes\n");
    free(selected);
    free(seen_array);
    return true;
}

void free_matrix(Matrix *matrix) {
    for (int i = 0; i < matrix->rows; i++) {
        free(matrix->data[i]);
    }
    free(matrix->data);
    free(matrix);
}

void free_blocks(Block *blocks, int num_blocks) {
    for (int i = 0; i < num_blocks; i++) {
        free(blocks[i].indices);
    }
    free(blocks);
}

/**
 * @brief Main function of the CFF verification program.
 * 
 * Reads a matrix from 'input.txt' and verifies if it is a valid d-CFF.
 * 
 * @return 0 on success, 1 on error.
 */
int main() {
    const char *file_path = "input.txt";
    
    Matrix *matrix = read_matrix_from_file(file_path);
    if (!matrix) {
        return 1;
    }
    
    int num_blocks;
    Block *blocks = process_columns(matrix, &num_blocks);

    int d = 2;
    is_cff(blocks, num_blocks, d);
    
    free_blocks(blocks, num_blocks);
    free_matrix(matrix);
    
    return 0;
}
