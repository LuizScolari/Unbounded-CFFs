"""
verify.py - Matrix sum verification script.

This script reads a binary matrix from a file and verifies if the sum of 1s
in each row and column matches the expected values.
"""

def verify_uniform_sum_matrix(file_path, expected_row_sum, expected_col_sum):
    matrix = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                matrix.append([int(x) for x in line.strip().split()])
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return
    except ValueError:
        print("Error: File contains non-integer characters.")
        return

    if not matrix:
        print("Matrix is empty.")
        return

    num_rows = len(matrix)
    num_cols = len(matrix[0])

    print("--- Checking Rows ---")
    check = True
    for i, row in enumerate(matrix):
        current_row_sum = sum(row)
        if current_row_sum != expected_row_sum:
            print(f"Row {i}: Expected: {expected_row_sum}, Found: {current_row_sum}")
            check = False
    if check == True:
        print("Valid!")
            
    check = True
    print("\n--- Checking Columns ---")
    for j in range(num_cols):
        current_col_sum = sum(matrix[i][j] for i in range(num_rows))
        if current_col_sum != expected_col_sum:
            print(f"Column {j}: Expected: {expected_col_sum}, Found: {current_col_sum}")
            check = False
    if check == True:
        print("Valid!")


def main():
    """
    Main function that runs the matrix verification.
    
    Reads a matrix from 'input.txt' and checks if all rows sum to expected_row_sum
    and all columns sum to expected_col_sum.
    """
    file_path = 'input.txt'
    expected_row_sum = 3
    expected_col_sum = 3

    print(f"Verifying matrix with expected sum of '{expected_row_sum}' for rows and '{expected_col_sum}' for columns.\n")

    verify_uniform_sum_matrix(
        file_path,
        expected_row_sum,
        expected_col_sum
    )


if __name__ == "__main__":
    main()
