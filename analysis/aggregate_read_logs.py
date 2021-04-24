import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import ast

def tied_preds(predictions, maxpos):
    max = predictions[maxpos]
    if predictions.count(max) > 1:
        return True
    else:
        return False

file_path = sys.argv[1]

# Iterating through read sets
prefix_count = 0
prefix = str(prefix_count) + "_"

predictions = []
matches = []

correct_counts = []
misclassified_counts = []
tied_counts = []
insufficient_counts = []
total_counts = []

## TODO - add option to detect cluster vs regular
# Regular files have two lines - score list + Source
# Cluster files have five lines - screen members, clusters, prediction, number of matches, source

clustered = False

### Get the number of members
all_files = glob.glob(file_path + '*.log')

test_file = all_files[0]
with open(test_file) as f:
    rt = f.readlines()

if len(rt) == 2:
    rs = rt[0].strip()
    r = ast.literal_eval(rs)
    members = len(r)
    classification_matrix = [[0 for x in range(members)] for y in range(members)]

else:
    clustered = True
    members = int(rt[0].strip())
    classification_matrix = [[0 for x in range(members)] for y in range(members)]

# While there are reads of this read set left to be classified
while (len(glob.glob(file_path + prefix  + '*.log')) > 0):

    # Get the files
    selectedFiles = glob.glob(file_path + prefix  + '*.log')

    readset_predictions = []
    readset_matches = []

    # Counters
    correct = 0
    misclassified = 0
    tied = 0
    insufficient = 0
    total = 0

    for filename in selectedFiles:
        total += 1

        with open(filename) as f:
            results_tmp = f.readlines()

        # TODO - split into cluster and Regular

        if clustered:
            pred = int(results_tmp[2].strip())
            score = int(results_tmp[3].strip())
            source = int(results_tmp[4].strip())

            readset_predictions.append(pred)
            readset_matches.append(score)

            if score == 0:
                insufficient += 1
            else:
                if pred == source:
                    classification_matrix[source][pred] += 1
                    correct += 1
                else:
                    classification_matrix[source][pred] += 1
                    misclassified += 1

        else:
            results_string = results_tmp[0].strip()
            results = ast.literal_eval(results_string)
            results = [int(x) for x in results]

            source = int(results_tmp[1].strip())

            pred = np.argmax(results)
            score = results[pred]
            readset_predictions.append(pred)
            readset_matches.append(results[source])

            if results[pred] == 0:
                insufficient += 1
            else:
                if pred == source:
                    classification_matrix[source][pred] += 1
                    correct += 1
                else:
                    classification_matrix[source][pred] += 1
                    misclassified += 1

                if tied_preds(results, pred):
                    tied += 1


    predictions.append(readset_predictions)
    matches.append(readset_matches)
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

# TODO - PLOTTING
