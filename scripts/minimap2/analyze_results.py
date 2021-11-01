import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# File with list of sequences per file
contigs_file = sys.argv[1]
contigs_prefix = ("Sequences in")

# File with classification breakdown
classification_file = sys.argv[2]
classification_prefix = (">")

with open(contigs_file) as f:
    contigs_lines = f.readlines()
contigs_lines = [x.strip() for x in contigs_lines]

# Build dictionary of correct sequences per file
genome_contigs = {}
curr_genome = contigs_lines[0].split()[len(contigs_prefix.split())]

for c in contigs_lines:
    if c.startswith(contigs_prefix):
        curr_genome = c.split()[len(contigs_prefix.split())]
        genome_contigs[curr_genome] = []
    else:
        genome_contigs[curr_genome].append(c.split()[0][1:])

# print(genome_contigs)

# Read classification results
with open(classification_file) as f:
    classification_lines = f.readlines()

# Strip
classification_lines = [x.strip() for x in classification_lines]

curr_file = ""
correct_contigs = []

# Track results
files = []
correct = []
curr_correct = 0
incorrect = []
curr_incorrect = 0


# For each line
for r in classification_lines:

    # If we're on a new organism
    if r.startswith(classification_prefix):

        # Save last correct/incorrect values
        if curr_correct > 0 or curr_incorrect > 0:
            correct.append(curr_correct)
            incorrect.append(curr_incorrect)

        curr_correct = 0
        curr_incorrect = 0

        # Get the name of this organism
        curr_file = r[1:]
        # Get the correct contigs for this genome
        correct_contigs = genome_contigs[curr_file]
        #
        files.append(curr_file)

    else:

        counts = int(r.split()[0])
        source = r.split()[1]
        if source in correct_contigs:
            curr_correct += counts
        else:
            curr_incorrect += counts
            # print(curr_file, source, counts)

# Add last counts
correct.append(curr_correct)
incorrect.append(curr_incorrect)

# Output
print(correct)
# print(incorrect)

np.save("correct_counts.npy", correct)
np.save("incorrect_counts.npy", incorrect)

# numerator = 0
# denominator = 0
# scores = []
# for i in range(len(correct)-22):
#     numerator += correct[i]
#     denominator += correct[i]
#     denominator += incorrect[i]
#     scores.append(correct[i]/(correct[i] + incorrect[i]))
#     print(correct[i], incorrect[i], correct[i]/(correct[i] + incorrect[i]))
# print(numerator, denominator, numerator/denominator)
# print(scores)
# # Counting how many have accuracy < threshold
# print(len([1 for i in scores if i < 0.7]))
