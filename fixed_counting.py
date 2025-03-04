import numpy as np
import itertools
import time
from collections import Counter 
import random
from bisect import bisect_left


def generate_matrices(vectors):
    # Get all permutations of each vector
    permutations = [list(itertools.permutations(vec)) for vec in vectors]
    # Get all combinations of permutations where each row is one permutation of a vector
    all_combinations = itertools.product(*permutations)

    # Convert each combination of permutations into a numpy matrix
    matrices = [np.array(combination) for combination in all_combinations]
    
    return matrices

# vectors = [
#     np.array([-1,2,3]),
#     np.array([-1,0,0])
# ]

def count_appearance_frequencies(matrix):
    # Number of rows in the matrix
    num_rows = matrix.shape[0]

    # Initialize a 1-indexed array of size num_rows (1 to num_rows inclusive)
    result = [0] * (num_rows + 1)
    
    # Iterate over each column
    for col in range(matrix.shape[1]):
        # Get the column values
        col_values = matrix[:, col]

        # don't marked elements
        col_values = col_values[col_values >= 0]
        # Count the frequency of each element in the column
        frequency_count = Counter(col_values)

        # Count how many elements appear exactly i times
        for value, count in frequency_count.items():
            result[count] += 1
    
    return result

def perm_parity(row):
    # Count the number of inversions in the row
    inversions = 0
    for i in range(len(row)):
        for j in range(i + 1, len(row)):
            if row[i] > row[j]:
                inversions += 1
    return inversions

def matrix_parity(M):
    overall_sign = 0
    for row in M:
        absrow=np.abs(row)
        overall_sign+=perm_parity(absrow)
    
    return (-1)**overall_sign

def all_possible_rows(vector):
    k = vector.size
    full_set = set(range(1, k + 1))
    vector_set = set(np.absolute(vector))
    difference = full_set - vector_set
    zero_indices = np.where(vector == 0)[0]
    B_permutations = itertools.permutations(difference)
    
    # Create a list to store all the possible vectors
    all_possible_vectors = []
    
    for perm in B_permutations:     
        A_copy = vector.copy()  # Make a copy of A
        A_copy[zero_indices] = perm  # Replace zeros with the current permutation
        all_possible_vectors.append(A_copy)  # Store the result

    return all_possible_vectors

def generate_fixed_matrices(vectors):
    # Get all permutations of each vector
    permutations = [all_possible_rows(vec) for vec in vectors]
    # Get all combinations of permutations where each row is one permutation of a vector
    all_combinations = itertools.product(*permutations)

    # Convert each combination of permutations into a numpy matrix
    matrices = [np.array(combination) for combination in all_combinations]
    
    return matrices

# works for k=4 or k=6
def num_core_mu_k(matrix,k):
    # Step 1: Find the intersection of all rows (elements that appear in all rows) to make sure it's a core element?
    common_elements = set(matrix[0])
    for row in matrix[1:]:
        common_elements.intersection_update(row)
    
    # Step 2: Check for elements that appear 4 times in any column
    common_elements_with_four_in_column = 0
    
    for element in common_elements:
        # For each common element, check each column
        for col in range(matrix.shape[1]):
            if np.sum(matrix[:, col] == element) == k:
                common_elements_with_four_in_column += 1
                break  
    
    return common_elements_with_four_in_column

def has_unique_element(matrix: np.ndarray) -> bool:
    """
    Returns True if there exists at least one column in the matrix
    where some positive element appears exactly once.
    Returns False if in every column all positive elements appear at least twice.
    """
    # Iterate over each column of the matrix
    for col in range(matrix.shape[1]):
        # Extract the current column
        current_col = matrix[:, col]
        # Filter to only consider positive elements (greater than 0)
        positive_elements = current_col[current_col > 0]
        
        # If there are positive elements in this column, count occurrences
        if positive_elements.size > 0:
            unique_vals, counts = np.unique(positive_elements, return_counts=True)
            # If any positive element appears exactly once, return True immediately
            if np.any(counts == 1):
                return True

    # If no column has a positive element appearing exactly once, return False
    return False

vectors = [
    np.array([-1,5,2,3,4,6]),
    np.array([-2,6,0,0,0,0]),
    np.array([3,-1,0,0,0,0]),
    np.array([4,-2,0,0,0,0]),
    np.array([3,5,0,0,0,0]),
    np.array([4,6,0,0,0,0])
]

# vectors = [
#     np.array([1,3,0,0]),
#     np.array([1,4,0,0]),
#     np.array([2,4,0,0]),
#     np.array([2,3,0,0]),
#     np.array([2,3,0,0]),
#     np.array([2,3,0,0])
# ]

matrices = generate_fixed_matrices(vectors)

pos_info={}
neg_info={}
total=0
print(vectors)
for i, matrix in enumerate(matrices):
    if has_unique_element(matrix):
        continue
    # coefs=count_appearance_frequencies(matrix)
    pm=matrix_parity(matrix)
    coefs=count_appearance_frequencies(matrix[:,2:])
    # print(matrix)
    
    # print('------')
    
    if pm>=1:
        total=total+3**coefs[4]
        if str(coefs[2:]) in pos_info.keys():
            pos_info[str(coefs[2:])]+=1

        else:
            pos_info[str(coefs[2:])]=1
    else:
        total=total-3**coefs[4]
        if str(coefs[2:]) in neg_info.keys():
            neg_info[str(coefs[2:])]+=1
        else:
            neg_info[str(coefs[2:])]=1

    # if coefs[1]==0 and coefs[3]==0 and coefs[5]==0 and coefs[6]==0 and num_core_mu_k(matrix[:,2:],4)==0:
        # further restricted
        # if matrix[2,0]==matrix[4,0] and matrix[3,0]==matrix[5,0] and matrix[0,1]==matrix[4,1] and matrix[1,1]==matrix[5,1]:
            # pm=matrix_parity(matrix)
            # # coefs=count_appearance_frequencies(matrix[2:,:])
            # print(matrix)
            # coefs=count_appearance_frequencies(matrix[:,2:])
            # print(coefs[2:])
            # print('------')
            
            # if pm>=1:
            #     if str(coefs[2:]) in pos_info.keys():
            #         pos_info[str(coefs[2:])]+=1
            #     else:
            #         pos_info[str(coefs[2:])]=1
            # else:
            #     if str(coefs[2:]) in neg_info.keys():
            #         neg_info[str(coefs[2:])]+=1
            #     else:
            #         neg_info[str(coefs[2:])]=1


print('positive',pos_info)
print('negative',neg_info)
print('total',total)
# Define each row as a separate NumPy array
# row1 = np.array([-1, 5, 2, 3, 4, 6])
# row2 = np.array([-2, 5, 3, 6, 1, 4])
# row3 = np.array([ 3, -1, 5, 6, 2, 4])
# row4 = np.array([ 4, -2, 3, 5, 1, 6])
# row5 = np.array([ 3, 6, 2, 5, 4, 1])
# row6 = np.array([ 4, 6, 5, 3, 2, 1])

# # Stack the rows into a single NumPy matrix
# matrix = np.array([row1, row2, row3, row4, row5, row6])


# print(matrix[:,2:])
# coefs=count_appearance_frequencies(matrix[:,2:])
# print(coefs[2:])
# print(num_core_mu_k(matrix[:,3:],4))