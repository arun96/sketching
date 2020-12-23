import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import pickle

def generate_kmers(sequence, k):
    kmer_counts = {}
    num_kmers = len(sequence) - k + 1
    for i in range(num_kmers):
        kmer = sequence[i:i+k]
        if kmer not in kmer_counts:
            kmer_counts[kmer] = 0
        kmer_counts[kmer] += 1
    return kmer_counts

def generate_kmers_list(sequence, k):
    kmers = []
    num_kmers = len(sequence) - k + 1
    for i in range(num_kmers):
        kmer = sequence[i:i+k]
        kmers.append(kmer)
    return kmers

def resemblance(a, b):
    top = len(list(set(a) & set(b)))
    bottom = len(list(set(a).union(b)))
    return float(top)/float(bottom)


def containment(a, b):
    top = len(list(set(a) & set(b)))
    bottom = len(a)
    return float(top)/float(bottom)


filename = sys.argv[1]
experiment = sys.argv[2]

with open(filename, 'rb') as fp:
    data_matrix = pickle.load(fp)

c_matrix = [[0 for x in range(len(data_matrix))] for y in range(len(data_matrix))]
r_matrix = [[0 for x in range(len(data_matrix))] for y in range(len(data_matrix))]

for x in range(0, len(data_matrix)):
    for y in range(0, len(data_matrix)):
        c_matrix[x][y] = containment(data_matrix[x], data_matrix[y])
        r_matrix[x][y] = resemblance(data_matrix[x], data_matrix[y])

print(len(c_matrix), len(c_matrix[0]))
print(len(r_matrix), len(r_matrix[0]))

with open('./' + experiment + '_containment.pkl', 'wb') as fp:
    pickle.dump(c_matrix, fp)

with open('./' + experiment + '_resemblance.pkl', 'wb') as fp2:
    pickle.dump(r_matrix, fp2)
