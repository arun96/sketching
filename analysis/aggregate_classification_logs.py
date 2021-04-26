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

prefix_count = 0
prefix = str(prefix_count) + "_"

predictions = []
matches = []

screen_size = 0

# File Formats
# Regular - just one line, array of values
# Cluster - four lines: number of members in screen, path, prediction, score

while (len(glob.glob(file_path + prefix  + '*.log')) > 0):

    selectedFiles = glob.glob(file_path + prefix  + '*.log')

    readset_predictions = []
    readset_matches = []

    for filename in selectedFiles:
        with open(filename) as f:
            results_tmp = f.readlines()

        # Regular
        if len(results_tmp) < 4:
            results_string = results_tmp[0].strip()
            # results = results_string.strip('][').split(',')
            results = ast.literal_eval(results_string)
            # print(results)
            results = [int(x) for x in results]

            screen_size = len(results)

            pred = np.argmax(results)
            score = results[pred]
            readset_predictions.append(pred)
            readset_matches.append(score)

        # Clustered
        else:
            # print(filename)
            screen_size = int(results_tmp[0].strip())
            pred = int(results_tmp[2].strip())
            score = int(results_tmp[3].strip())
            readset_predictions.append(pred)
            readset_matches.append(score)

    predictions.append(readset_predictions)
    matches.append(readset_matches)
    prefix_count += 1
    prefix = str(prefix_count) + "_"


### PLOTTING

# 1. Classification per read set
for i in range(0, prefix_count):
    plt.hist(predictions[i])
    # plt.hist(predictions[i], bins = list(range(0, screen_size)))
    plt.xlabel('Predicted Genome')
    plt.ylabel('Frequency')
    plt.xlim([0, prefix_count])
    plt.title("Readset " + str(i))
    plt.show()

# 2. Classification as a heatmap? Readsets vs input genomes (just to see classification as a matrix)

classification_matrix = [[0.0 for i in range(screen_size)] for j in range(prefix_count)]

for i in range(prefix_count):
    for j in range(screen_size):
        classification_matrix[i][j] = predictions[i].count(j)/len(predictions[i])

# print(classification_matrix)

# Super simple one:
# plt.imshow(classification_matrix, cmap='hot', interpolation='nearest')
# plt.xlabel('Predicted Genome')
# plt.ylabel('Readset')
# plt.show()

# Seaborn
ax = sns.heatmap(classification_matrix, square = True, cmap="Blues", annot=True)
plt.title("Classification Matrix")
plt.xlabel('Predicted Genome')
plt.ylabel('Readset')
plt.show()

# 3. Distribution of matches
