import os
import sys
import random
from heapq import nsmallest
import itertools
import numpy as np
import mmh3
import re
import pickle


pickle_flag = False

kmer_size = 21
# Default - if we want to dynamically adjust this, we enable the flag below
sketch_size = 10000

# varied sketch sizes for the genomes
calculate_sketch = True
# target number of matches - used when calculating sketch size
target_matches = 30

# TODO
weighted = False

def generate_kmers(sequence, k):
    kmer_counts = {}
    num_kmers = len(sequence) - k + 1
    for i in range(num_kmers):
        kmer = sequence[i:i+k]
        if kmer not in kmer_counts:
            kmer_counts[kmer] = 0
        kmer_counts[kmer] += 1
    return kmer_counts

def get_reads(read_file):
    reads = []
    with open(read_file) as f:
        for line in itertools.islice(f, 1, None, 2):
            r = line.rstrip()
            # removing all non-ACGTN characters
            r = r.upper()
            r = re.sub('[^ACTGN]','', r)
            reads.append(r)
    return reads

def get_sketch_size(genome, target_matches, read_len, read_err, kmer_size):
    multiplier = (1 - read_err)**21
    denom = read_len/(target_matches/multiplier)
    num_sketches = len(genome)/denom
    # return num_sketches
    num_sketches = (target_matches * len(genome))/(read_len*multiplier)
    return num_sketches

# Num sketches = Genome size/ (read Len /(Num Matches/multiplier)

def max_index(input_list):
    return np.argmax(input_list)

# TODO - Why is there a Y??? Or a ;????
reverse_complement = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N', 'Y':'Y', 'F':'F', ';':';','S':'S', 'R':'R'}[B] for B in x][::-1])

def canonical(a, b):
    if a < b:
        return a
    else:
        return b


microbes = ['Bacillus_subtilis',
    'Cryptococcus_neoformans',
    'Enterococcus_faecalis',
    'Escherichia_coli',
    'Lactobacillus_fermentum',
    'Listeria_monocytogenes',
    'Pseudomonas_aeruginosa',
    'Saccharomyces_cerevisiae',
    'Salmonella_enterica',
    'Staphylococcus_aureus'
]

genomes = ['./Genomes/Bacillus_subtilis_complete_genome.fasta',
    './Genomes/Cryptococcus_neoformans_draft_genome.fasta',
    './Genomes/Enterococcus_faecalis_complete_genome.fasta',
    './Genomes/Escherichia_coli_complete_genome.fasta',
    './Genomes/Lactobacillus_fermentum_complete_genome.fasta',
    './Genomes/Listeria_monocytogenes_complete_genome.fasta',
    './Genomes/Pseudomonas_aeruginosa_complete_genome.fasta',
    './Genomes/Saccharomyces_cerevisiae_draft_genome.fasta',
    './Genomes/Salmonella_enterica_complete_genome.fasta',
    './Genomes/Staphylococcus_aureus_complete_genome.fasta'
]


# one microbe
# microbes = ['Bacillus_subtilis']
# genomes = ['./Genomes/Bacillus_subtilis_complete_genome.fasta']

# read directory
experiment = './reads/reads_10000_looped_1/'

# TODO - change this for each run
read_len = 10000
read_err = 0.01

# hashes = []
# stores sketches for each genome
sketch = []

#TODO: add option for weighting
for g in genomes:
    genome_sketch = []
    genome = ""
    if g == './Genomes/Escherichia_coli_complete_genome.fasta' or  g == './Genomes/Salmonella_enterica_complete_genome.fasta':
        # print(g)
        with open(g) as f:
            g_temp = f.readlines()
            genome1 = g_temp[1].rstrip()
            genome2 = g_temp[3].rstrip()
            genome = genome1 + genome2
    elif g == './Genomes/Cryptococcus_neoformans_draft_genome.fasta' or g == './Genomes/Saccharomyces_cerevisiae_draft_genome.fasta' or './Genomes/Staphylococcus_aureus_complete_genome.fasta':
        with open(g) as f:
            g_temp = f.readlines()
            for t in g_temp:
                if t[0] != '>':
                    genome += t.rstrip()
    else:
        with open(g) as f:
            genome = f.readlines()[1].rstrip()

    genome_kmer_breakdown = generate_kmers(genome, kmer_size)

    # hash_values = [hash(kmer) for kmer in list(genome_kmer_breakdown)]
    hash_values = [hash(canonical(kmer, reverse_complement(kmer))) for kmer in list(genome_kmer_breakdown)]

    # default version
    # genome_sketch = nsmallest(sketch_size, hash_values)

    if calculate_sketch:
        calculated_sketch_size = int(get_sketch_size(genome, target_matches, read_len, read_err, kmer_size))
        print(g, len(genome), calculated_sketch_size)
        genome_sketch = nsmallest(calculated_sketch_size, hash_values)
    else:
        print(g, len(genome), sketch_size)
        genome_sketch = nsmallest(sketch_size, hash_values)


    # hashes.append(hashes)
    sketch.append(genome_sketch)

# get read hashes
all_reads = []
# all_reads_hashes = []

# get all the reads
for microbe in microbes:
    x = experiment + microbe + '.fasta'
    reads = get_reads(x)
    all_reads.append(reads)

read_matches_all = []
correct_counts = []
mis_counts = []
insuf_counts = []
# get all hashes of all reads
# TODO - compare against the hashes of genome
for r in range(0, len(all_reads)):
    print(microbes[r])
    correct = r
    correct_count = 0
    mis_count = 0
    insuf_count = 0
    read_matches = []

    readset = all_reads[r]
    # readset_hashes = []
    for read in readset:
        read_kmers = generate_kmers(read, 21)
        # read_hashes = [hash(kmer) for kmer in list(read_kmers)]
        read_hashes = [hash(canonical(kmer, reverse_complement(kmer))) for kmer in list(read_kmers)]

        hash_matches = [list(set(read_hashes) & set(g_s)) for g_s in sketch]
        count_matches = [len(match) for match in hash_matches]
        # print(count_matches)
        pred = max_index(count_matches)
        if correct == pred and count_matches[pred] != 0:
            correct_count += 1
        elif count_matches[pred] == 0:
            insuf_count += 1

            # debugging
            # print(count_matches, pred, len(read), read)
        else:
            mis_count += 1

            #debugging
            # print(count_matches, pred, len(read), read)

        read_matches.append(count_matches[pred])

    read_matches_all.append(read_matches)
    correct_counts.append(correct_count)
    mis_counts.append(mis_count)
    insuf_counts.append(insuf_count)

    print(correct_count, len(readset), mis_count, insuf_count)

    # readset_hashes.append(read_hashes)
    # all_reads_hashes.append(readset_hashes)

if pickle_flag:
    with open("./dynamic_results/30_matches_looped_10/read_matches_30.pkl", "wb") as fp1:
        pickle.dump(read_matches_all, fp1)

    with open("./dynamic_results/30_matches_looped_10/correct_counts_30.pkl", "wb") as fp2:
        pickle.dump(correct_counts, fp2)

    with open("./dynamic_results/30_matches_looped_10/mis_counts_30.pkl", "wb") as fp3:
        pickle.dump(mis_counts, fp3)

    with open("./dynamic_results/30_matches_looped_10/insuf_counts_30.pkl", "wb") as fp4:
        pickle.dump(insuf_counts, fp4)
