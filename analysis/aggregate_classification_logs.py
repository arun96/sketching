import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import ast

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
            pred = int(results_tmp[2].strip())
            score = int(results_tmp[3].strip())
            readset_predictions.append(pred)
            readset_matches.append(score)

    predictions.append(readset_predictions)
    matches.append(readset_matches)
    prefix_count += 1
    prefix = str(prefix_count) + "_"


### PLOTTING

# Plot of what the reads are classified as
for i in range(0, prefix_count):
    plt.hist(predictions[i])
    # plt.hist(predictions[i], bins = list(range(0, screen_size)))
    plt.xlabel('Predicted Genome')
    plt.ylabel('Frequency')
    plt.xlim([0, prefix_count])
    plt.title("Readset " + str(i))
    plt.show()
