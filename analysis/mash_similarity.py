import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from fractions import Fraction

filename = sys.argv[1]

with open(filename) as f:
    mash_results = f.readlines()

# print(len(mash_results))

num_hashes = 1000.0

# Infer number of elements in file
first_genome = mash_results[0].split()[0]

second_size = 0

for line in mash_results:
    if line.split()[0] != first_genome:
        break
    else:
        second_size += 1

first_size =  int(len(mash_results) / second_size)
# print(first_size, second_size)

mash_matrix = [[0.0 for x in range(second_size)] for y in range(first_size)]
# print(len(mash_matrix), len(mash_matrix[0]))

# Fill out this matrix
for i in range(first_size):
    for j in range(second_size):
        # Value for [i][j] is found at line i * second_size + j
        pos = (i * second_size) + j
        matches = mash_results[pos].split()[-1]
        mash_matrix[i][j] = float(Fraction(matches)) * num_hashes
        # print(matches_val)

mash_matrix = np.asarray(mash_matrix)

# TODO - create annotation mask for interesting genomes
heatmap_annotations = [["" for x in range(second_size)] for y in range(first_size)]
for i in range(first_size):
    for j in range(second_size):
        if mash_matrix[i][j] > (0.001*num_hashes) and mash_matrix[i][j] < (num_hashes-5):
            heatmap_annotations[i][j] = int(mash_matrix[i][j])
heatmap_annotations = np.asarray(heatmap_annotations)
# TODO - organism names? Instead of just numbers


# Heatmap
# all vs all
if (first_size == second_size):
    ax = sns.heatmap(mash_matrix, square = True, cmap="Blues", fmt = "s", annot = heatmap_annotations)
    # ax = sns.heatmap(mash_matrix, square = True, cmap="Blues", annot=True)
    plt.title("Mash Similarity (" + str(int(num_hashes)) + " hashes)")
    plt.show()
else:
    # ax = sns.heatmap(mash_matrix, square = True, cmap="Blues")
    ax = sns.heatmap(mash_matrix, square = True, cmap="Blues", fmt = "s", annot = heatmap_annotations)
    # ax = sns.heatmap(mash_matrix, square = True, cmap="Blues", annot=True)
    plt.title("Mash Similarity (" + str(int(num_hashes)) + " hashes)")
    plt.xlabel('Community #1')
    plt.ylabel('Community #2')
    plt.show()

# Log of Heatmap
# mash_matrix_log = np.log(mash_matrix)
# ax = sns.heatmap(mash_matrix_log, square = True, cmap="Blues")
# # ax = sns.heatmap(mash_matrix, square = True, cmap="Blues", annot=True)
# plt.title("Log(Mash Similarity)")
# plt.xlabel('Community #1')
# plt.ylabel('Community #2')
# plt.show()
