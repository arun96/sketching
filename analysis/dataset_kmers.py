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

# ANALYSIS

path = sys.argv[1]
ext = sys.argv[2]
experiment = sys.argv[3]

files = glob.glob(path + "*." + ext)
filenames = [x[len(path):(len(x) - len(ext) - 1)] for x in files]

genomes = []

# Get Genomes
for file in files:
    with open(file) as f:
        all_lines = f.readlines()
    all_lines = [x.rstrip() for x in all_lines]
    genome_lines = [x for x in all_lines if x[0] != '>']
    genome = "".join(genome_lines)
    genomes.append(genome)

print(len(genomes[0]))

# genomes = [[hash(k) for k in list(generate_kmers(g, 21))] for g in genomes]
genomes = [[hash(k) for k in generate_kmers_list(g, 21)] for g in genomes]
print(len(genomes))

print(len(genomes[0]))

with open('./' + experiment + '_genome_kmers.pkl', 'wb') as fp:
    pickle.dump(genomes, fp)
