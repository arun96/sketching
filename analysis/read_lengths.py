# Plotting
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# I/O
import glob
import ast
import os
import sys

file_path = sys.argv[1]

ext = '.fasta'

all_files = glob.glob(file_path + '*' + ext)

read_lengths = []
for read_file in all_files:
    with open(read_file) as f:
        f = open(read_file)
        i = 0
        for line in f.readlines():
            if i % 2 == 0 :
                s = line.split('_')
                read_len = int(s[-2])
                read_lengths.append(read_len)
            i += 1

# Plotting
target_length = 10000
hist_range = 500
steps = 10
plt.hist(read_lengths, list(range(target_length-hist_range, target_length+hist_range, steps)))
plt.xlabel('Read Length')
plt.ylabel('Frequency')
plt.title("Read Length Distribution")
plt.show()
