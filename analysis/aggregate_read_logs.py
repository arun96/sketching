import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import ast

def tied_preds(predictions, maxpos):
    max = predictions[maxpos]
    if predictions.count(maxpos) > 1:
        return True
    else:
        return False

file_path = sys.argv[1]

prefix_count = 0
prefix = str(prefix_count) + "_"

predictions = []

screen_size = 0

correct_counts = []
misclassified_counts = []
tied_counts = []
insufficient_counts = []
total_counts = []

### Get the number of members
all_files = glob.glob(file_path + '*.log')
test_file = all_files[0]
with open(test_file) as f:
    rt = f.readlines()
rs = rt[0].strip()
r = ast.literal_eval(rs)
members = len(r)
classification_matrix = [[0 for x in range(members)] for y in range(members)]

while (len(glob.glob(file_path + prefix  + '*.log')) > 0):

    selectedFiles = glob.glob(file_path + prefix  + '*.log')

    readset_predictions = []

    correct = 0
    misclassified = 0
    tied = 0
    insufficient = 0
    total = 0

    for filename in selectedFiles:
        total += 1
        with open(filename) as f:
            results_tmp = f.readlines()

        results_string = results_tmp[0].strip()
        # results = results_string.strip('][').split(',')
        results = ast.literal_eval(results_string)
        # print(results)
        results = [int(x) for x in results]

        source = int(results_tmp[1].strip())

        screen_size = len(results)

        pred = np.argmax(results)
        score = results[pred]
        readset_predictions.append(pred)

        if results[pred] == 0:
            insufficient += 1
        else:
            if pred == source:
                classification_matrix[source][pred] += 1
                correct += 1
            else:
                misclassified += 1

            if tied_preds(results, pred):
                tied += 1


    predictions.append(readset_predictions)
    total_counts.append(total)
    correct_counts.append(correct)
    misclassified_counts.append(misclassified)
    tied_counts.append(tied)
    insufficient_counts.append(insufficient)

    prefix_count += 1
    prefix = str(prefix_count) + "_"

# Output summary
print("Total Reads, # Correctly Classified, # Misclassified, # Unclassified, # of Ties")
for i in range(0, prefix_count):
    print("Readset " + str(i) + ":")
    print(total_counts[i], correct_counts[i], misclassified_counts[i], insufficient_counts[i], tied_counts[i])

print("Classification Matrix:")
print(classification_matrix)
