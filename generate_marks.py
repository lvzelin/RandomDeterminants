import numpy as np
from itertools import combinations,permutations,product
import math

def count_columns_with_odd_zeros(table):
    return np.sum(np.sum(table == 0, axis=0) % 2 == 1)

def is_table_valid(table_input):

    max_element = np.max(table_input)
    num_occurrences = np.zeros(max_element, dtype=int)
    
    for i in range(1, max_element + 1):
        num_occurrences[i-1] = np.sum(table_input == i)
    
    odd_columns = count_columns_with_odd_zeros(table_input)
    num_triples = np.sum(num_occurrences == 3) + np.sum(num_occurrences == 1)
    num_extra_columns = num_triples - odd_columns
    
    if num_extra_columns % 2 == 1 or num_extra_columns < 0:
        # print('EXCEPTION: extra column')
        return False
    else:
        num_extra_columns //= 2
        return True

# a matrix is sequential if num of element i >= num of element j if i<j 
def is_sequential(matrix):
    n=matrix.max()
    for i in range(1,n):
        if np.sum(matrix == i) < np.sum(matrix == i+1):
            return False

    return True

# two matrices are isomorphic if we can obtain one from the other by permuting rows, permuting columns and permuting elements
def are_isomorphic(matrix1, matrix2):
    if matrix1.shape != matrix2.shape:
        return False
    
    rows, cols = matrix1.shape
    unique_elements = np.unique(matrix1)

    for row_perm in permutations(range(rows)):
        permuted_rows_matrix1 = matrix1[row_perm, :]
        
        for col_perm in permutations(range(cols)):
            permuted_matrix1 = permuted_rows_matrix1[:, col_perm]
            for elem_perm in permutations(unique_elements):
                mapping = {old: new for old, new in zip(unique_elements, elem_perm)}
                transformed_matrix = np.vectorize(mapping.get)(permuted_matrix1)
                
                if np.array_equal(transformed_matrix, matrix2):
                    return True
    return False

# generate all the isomorphic matrices of the input matrix
def generate_isomorphic_matrices(matrix):
    rows, cols = matrix.shape
    unique_elements, counts = np.unique(matrix, return_counts=True)
    
    # Group elements by their counts
    element_groups = {count: [] for count in counts}
    for element, count in zip(unique_elements, counts):
        if element != 0:
            element_groups[count].append(element)
    
    isomorphic_matrices = set()
    
    for row_perm in permutations(range(rows)):
        permuted_rows_matrix = matrix[row_perm, :]
        
        for col_perm in permutations(range(cols)):
            permuted_matrix = permuted_rows_matrix[:, col_perm]
            
            # Permute non-zero elements within groups of the same count
            for count, elements in element_groups.items():
                for elem_perm in permutations(elements):
                    mapping = {old: new for old, new in zip(elements, elem_perm)}
                    transformed_matrix = np.vectorize(lambda x: mapping.get(x, x))(permuted_matrix)
                    
                    isomorphic_matrices.add(tuple(map(tuple, transformed_matrix)))
    
    return [np.array(matrix) for matrix in isomorphic_matrices]


def is_isomorphic_to_any(matrix, matrix_list):
    for ref_matrix in matrix_list:
        if are_isomorphic(matrix, ref_matrix):
            return True
    return False

def matrix_to_type(matrix, matrix_list):
    i = 0
    for ref_matrix in matrix_list:
        if are_isomorphic(matrix, ref_matrix):
            return i
        i=i+1
    return -1

# first generate all the possible n marked tables
def generate_all(n):
    rows = 6
    all_matrices = []
    for cols in range(1,n+1):
    
        # Generate all possible positions to place n nonzero elements in a matrix with 6 rows and n columns
        positions = list(combinations(range(rows * cols), n))
        
        # print(positions)
        # Generate all possible values combinations for the nonzero elements
        value_combinations = list(product(range(1, n + 1), repeat=n))
        for pos in positions:
            for values in value_combinations:
                # Create a flat matrix with zeros
                flat_matrix = [0] * (rows * cols)
                for i, val in zip(pos, values):
                    flat_matrix[i] = val
                
                # Reshape the flat matrix into the required 6-row matrix
                matrix = np.array(flat_matrix).reshape(rows, cols)
                
                # Check if there is any column with only zeros
                if any(np.all(matrix[:, col] == 0) for col in range(cols)):
                    continue
            
                # Check if any row contains more than one nonzero element
                if any(np.count_nonzero(matrix[row, :]) > 1 for row in range(rows)):
                    continue

                # num of element i >= num of element j if i<j
                if not is_sequential(matrix):
                    continue

                if not is_table_valid(matrix):
                    continue


                all_matrices.append(matrix)
    
    return all_matrices

def get_extra_coef(M):
    # Flatten the array to a 1D array
    flat_M = M.flatten()
    flat_M = flat_M[flat_M != 0]

    # Count the occurrences of each element
    unique_elements, counts = np.unique(flat_M, return_counts=True)
    
    element_perm=1
    for i in range(6):
        elements_appearing_i_times = unique_elements[counts == i]
        num_elements_appearing_i_times = len(elements_appearing_i_times)
        element_perm=element_perm*math.factorial(num_elements_appearing_i_times)

    num_cols=M.shape[1]
    return element_perm*math.factorial(num_cols)


# Main function
table_input = np.transpose(np.array([
    [1,0,0,1,0,0],[1,0,2,0,0,0]
]))

# Define the matrices
one_mark_list = [
    np.array([[1, 0, 0, 0, 0, 0]]).T,
]

# [15, 30, 60]
two_mark_list = [
    np.array([[1, 1, 0, 0, 0, 0]]).T,
    np.array([[1, 2, 0, 0, 0, 0]]).T,
    np.array([[1, 0, 0, 0, 0, 0], [0, 2, 0, 0, 0, 0]]).T,
]

# [20, 60, 120, 120, 120, 240, 720, 720]
three_mark_list = [
    np.array([[1, 1, 1, 0, 0, 0]]).T,
    np.array([[1, 1, 2, 0, 0, 0]]).T,
    np.array([[1, 2, 3, 0, 0, 0]]).T,
    np.array([[1, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0]]).T,
    np.array([[1, 1, 0, 0, 0, 0], [0, 0, 2, 0, 0, 0]]).T,
    np.array([[1, 2, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0]]).T,
    np.array([[1, 2, 0, 0, 0, 0], [0, 0, 3, 0, 0, 0]]).T,
    np.array([[1, 0, 0, 0, 0, 0], [0, 2, 0, 0, 0, 0], [0, 0, 3, 0, 0, 0]]).T
] 

four_mark_list = [
    np.array([[1, 1, 1, 1, 0, 0]]).T,
    np.array([[1, 1, 1, 2, 0, 0]]).T,
    np.array([[1, 1, 2, 3, 0, 0]]).T,
    np.array([[1, 2, 3, 4, 0, 0]]).T,
    np.array([[1, 1, 0, 0, 0, 0], [0, 0, 1, 1, 0, 0]]).T,
    np.array([[1, 1, 0, 0, 0, 0], [0, 0, 2, 2, 0, 0]]).T,
    np.array([[1, 2, 0, 0, 0, 0], [0, 0, 1, 2, 0, 0]]).T,
    np.array([[1, 1, 0, 0, 0, 0], [0, 0, 2, 3, 0, 0]]).T,
    np.array([[1, 2, 0, 0, 0, 0], [0, 0, 1, 3, 0, 0]]).T,
    np.array([[1, 1, 2, 0, 0, 0], [0, 0, 0, 3, 0, 0]]).T,
    np.array([[1, 2, 3, 0, 0, 0], [0, 0, 0, 1, 0, 0]]).T,
    np.array([[1, 2, 3, 0, 0, 0], [0, 0, 0, 4, 0, 0]]).T,
    np.array([[1, 2, 0, 0, 0, 0], [0, 0, 3, 4, 0, 0]]).T,
    np.array([[1, 1, 0, 0, 0, 0], [0, 0, 2, 0, 0, 0], [0, 0, 0, 3, 0, 0]]).T,
    np.array([[1, 2, 0, 0, 0, 0], [0, 0, 3, 0, 0, 0], [0, 0, 0, 4, 0, 0]]).T,
    np.array([[1, 2, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 3, 0, 0]]).T,
    np.array([[1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 2, 3, 0, 0]]).T,
    np.array([[1, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 2, 0, 0]]).T,
    np.array([[1, 2, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0], [0, 0, 0, 1, 0, 0]]).T,
    np.array([[1, 0, 0, 0, 0, 0], [0, 2, 0, 0, 0, 0], [0, 0, 3, 0, 0, 0], [0, 0, 0, 4, 0, 0]]).T,
    np.array([[1, 1, 1, 0, 0, 0], [0, 0, 0, 2, 0, 0]]).T,
    np.array([[1, 1, 2, 0, 0, 0], [0, 0, 0, 1, 0, 0]]).T,
    np.array([[1, 1, 2, 2, 0, 0]]).T,
    np.array([[1, 1, 0, 0, 0, 0], [0, 0, 1, 2, 0, 0]]).T
]

five_mark_list = [
    # at least have two columns?
    np.array([[1, 1, 1, 1, 0, 0], [0, 0, 0, 0, 2, 0]]).T,
    np.array([[1, 1, 1, 2, 0, 0], [0, 0, 0, 0, 3, 0]]).T,
]

mark_lists=[[],one_mark_list,two_mark_list,three_mark_list,four_mark_list,five_mark_list]

n = 5
all_matrices = generate_all(n)
coefs=[]
type_list=[0]*len(mark_lists[n])
coefs_list=[0]*len(mark_lists[n])
total_type_automorphism=[]
i=0
for type_matrix in mark_lists[n]:
    print(i)
    num_automorphism=len(generate_isomorphic_matrices(type_matrix))
    type_list[i]=num_automorphism
    coefs_list[i]=num_automorphism/get_extra_coef(type_matrix)
    i=i+1
    
# for matrix in matrices:
#     i=i+1
#     if i%100==0:
#         print(i)
#     type_matrix=matrix_to_type(matrix,mark_lists[n])
#     if type_matrix <0:
#         print(matrix)
#         mark_lists[n].append(matrix)
#         type_list.append[1]
#     else:
#         type_list[type_matrix]=type_list[type_matrix]+1

print('number of all possible matrices: ',len(all_matrices))
print('total number from automorphisms: ',sum(type_list))
print('raw numbers:')
print(type_list)
print('coefficients:')
print(coefs_list)




