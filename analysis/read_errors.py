# Plotting
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# I/O
import glob
import ast
import os
import sys


reads_file = sys.argv[1]

read_len = float(sys.argv[2])
read_len_KB = int(read_len/1000)

read_error = float(sys.argv[3])

expected_errors = read_len * (read_error/100)
print(expected_errors)

with open(reads_file) as f:
    reads_data = f.readlines()
# print(len(reads_data))

cigars = reads_data[2::3]
# print(len(cigars))

errors = []
SNP = []
indel = []

for c in cigars:
    s = c.count("X")
    id = c.count("D") + c.count("I")
    errors.append(s + id)
    SNP.append(s)
    indel.append(id)

bins = np.arange(2*expected_errors) - 0.5
plt.hist(errors, bins)
plt.xlabel("Errors per Read")
plt.ylabel("Frequency")
plt.title("# of Errors per Read - " + str(read_len_KB) + "KB, " + str(read_error) + "% Error")
plt.axvline(expected_errors, color='k', linestyle='dashed', linewidth=1)
min_ylim, max_ylim = plt.ylim()
plt.text(expected_errors*1.1, max_ylim*0.9, 'Expected Number of Errors: {:.2f}'.format(expected_errors))
plt.show()
