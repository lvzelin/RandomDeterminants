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

t = time.time()

vectors = [
    np.array([-1,2,3,4]),
    np.array([-1,2,3,4]),
    np.array([1,-2,3,4]),
    np.array([1,-2,3,4]),
    np.array([1,2,3,4]),
    np.array([1,2,3,4])
]



# print('parity: ',matrix_parity(matrix))

# Generate all matrices
matrices = generate_matrices(vectors)
elapsed = time.time() - t
print('time generation: ',elapsed)
# random.shuffle(matrices)

pos_info={}
neg_info={}


total=0
subcase_total=0
for i, matrix in enumerate(matrices):
    coefs=count_appearance_frequencies(matrix)
    if coefs[1]==0 and coefs[3]+coefs[5]==0 and coefs[6]==0:        
        if np.array_equal(matrix[:, 0], [-1,-1,-2,-2,1,1]):
            
            pm=matrix_parity(matrix)
            if pm>=1:
                if str(coefs[2:]) in pos_info.keys():
                    pos_info[str(coefs[2:])]+=1
                else:
                    pos_info[str(coefs[2:])]=1
            else:
                if str(coefs[2:]) in neg_info.keys():
                    neg_info[str(coefs[2:])]+=1
                else:
                    neg_info[str(coefs[2:])]=1
            if str(coefs[2:])=='[4, 0, 3, 0, 0]' and matrix[0,1]==1:
                subcase_total+=1
                print('------------------------------')
                print(subcase_total)
                print(i,'- matrix')
                print(matrix)
                print('parity: ',matrix_parity(matrix))

        
            

elapsed = time.time() - t
# print('constant weight: ',total)
# print('positive mu34: ',pos_num_mu34)
# print('negative mu34: ',neg_num_mu34)
print('positive: ',pos_info)
print('negative: ',neg_info)
print('time total: ',elapsed)
