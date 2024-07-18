import numpy as np
from itertools import combinations,permutations,product

def count_columns_with_odd_zeros(table):
    return np.sum(np.sum(table == 0, axis=0) % 2 == 1)

def is_table_valid(table_input):

    max_element = np.max(table_input)
    num_occurrences = np.zeros(max_element, dtype=int)
    
    for i in range(1, max_element + 1):
        num_occurrences[i-1] = np.sum(table_input == i)
    print(num_occurrences)
    odd_columns = count_columns_with_odd_zeros(table_input)
    num_triples = np.sum(num_occurrences == 3) + np.sum(num_occurrences == 1)
    num_extra_columns = num_triples - odd_columns
    
    if num_extra_columns % 2 == 1 or num_extra_columns < 0:
        print('EXCEPTION: extra column')
        return False
    else:
        num_extra_columns //= 2
        return True

# first generate all the possible n marked tables
def generate_all(n):
    rows = 6
    all_matrices = []
    for cols in range(0,n+1):
     
        # Generate all possible positions to place n nonzero elements in a matrix with 6 rows and n columns
        positions = list(combinations(range(rows * cols), n))
        
        print(positions)
        # Generate all possible values combinations for the nonzero elements
        value_combinations = list(product(range(1, n + 1), repeat=n))
        print(value_combinations)
        for pos in positions:
            for values in value_combinations:
                # Create a flat matrix with zeros
                flat_matrix = [0] * (rows * cols)
                for i, val in zip(pos, values):
                    flat_matrix[i] = val
                
                # Reshape the flat matrix into the required 6-row matrix
                matrix = np.array(flat_matrix).reshape(rows, cols)
                
                # Check if there is any column with only zeros
                if not any(np.all(matrix[:, col] == 0) for col in range(cols)):
                    all_matrices.append(matrix)
    
    return all_matrices


# Main function
table_input = np.transpose(np.array([
    [1,1,1,0,0,0],[0,0,0,1,0,0]
]))

result = is_table_valid(table_input)
print(table_input)
print(result)

n = 2
matrices = generate_all(n)
print(len(matrices))
# for i, matrix in enumerate(matrices):
#     print(f"Matrix {i+1}:\n{matrix}\n")


