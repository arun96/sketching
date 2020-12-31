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

screen_size = 0

while (len(glob.glob(file_path + prefix  + '*.log')) > 0):

    selectedFiles = glob.glob(file_path + prefix  + '*.log')

    readset_predictions = []

    for filename in selectedFiles:
        with open(filename) as f:
            results_tmp = f.readlines()

        results_string = results_tmp[0].strip()
        # results = results_string.strip('][').split(',')
        results = ast.literal_eval(results_string)
        # print(results)
        results = [int(x) for x in results]

        screen_size = len(results)

        pred = np.argmax(results)
        score = results[pred]
        readset_predictions.append(pred)

    predictions.append(readset_predictions)
    prefix_count += 1
    prefix = str(prefix_count) + "_"

# Plot of what the reads are classified as
for i in range(0, prefix_count):
    plt.hist(predictions[i])
    # plt.hist(predictions[i], bins = list(range(0, screen_size)))
    plt.xlabel('Predicted Genome')
    plt.ylabel('Frequency')
    plt.title("Readset " + str(i))
    plt.show()
