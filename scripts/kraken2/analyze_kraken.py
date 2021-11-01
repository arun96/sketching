import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
from collections import Counter

def get_read_prediction(assignments):
    read_results = {}
    for a in assignments:
        a_list = a.split(":")
        tax_id = str(a_list[0])
        count = int(a_list[1])
        if tax_id not in read_results:
            read_results[tax_id] = count
        else:
            read_results[tax_id] = read_results[tax_id] + count

    # Discount unclassified
    read_results.pop('0', None)
    read_results.pop('A', None)

    if len(read_results) == 0:
        prediction = 0
    else:
        prediction = max(read_results, key=read_results.get)
    return read_results, prediction

def get_file_predictions(file):
    file_results = []

    for read in file:
        # get the k-mer predictions
        assignments = list(read.split())[4:]
        classification_dict, prediction = get_read_prediction(assignments)
        file_results.append(prediction)

    return file_results

# File with list of sequences per file
kraken_folder = sys.argv[1]
ext = ".log"

# Iterate through files
kraken_files = glob.glob(kraken_folder + '*' + ext)
# print(kraken_files)

predictions = []

for reads_file in kraken_files:
    with open(reads_file) as f:
        reads_lines = f.readlines()
    reads_lines = [x.strip() for x in reads_lines]

    read_predictions = get_file_predictions(reads_lines)
    predictions.append(read_predictions)
    # print(Counter(read_predictions))

np.save("kraken_classification.npy", predictions)
