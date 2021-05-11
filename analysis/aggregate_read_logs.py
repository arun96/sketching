# Plotting
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# I/O
import glob
import ast
import os
import sys

# MISC
import copy

def tied_preds(predictions, maxpos):
    max = predictions[maxpos]
    if predictions.count(max) > 1:
        return True
    else:
        return False

file_path = sys.argv[1]

target = sys.argv[2]
error_rate = sys.argv[3]

save_flag = False
if len(sys.argv) == 5:
    save_flag = True
    save_dir = sys.argv[4]


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
accuracy = []

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
    accuracy.append(correct/total)

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

# 1. Plots across the organisms (combined plot)
# 1a - Number of organisms vs accuracy

# plt.hist(accuracy)
# plt.xlabel('Accuracy')
# plt.ylabel('# of Organisms')
# plt.title("Classification Accuracy across all genomes")
# plt.show()

# 2. Number of incorrect attached reads (look at classification matrix) - could be done as a heatmap

# Seaborn Heatmap

classification_fraction = copy.deepcopy(classification_matrix)

for i in range(len(classification_fraction)):
    nr = sum(classification_fraction[i])
    for j in range(len(classification_fraction[i])):
        classification_fraction[i][j] = classification_fraction[i][j]/nr

# Read Counts
# ax = sns.heatmap(classification_matrix, square = True, cmap="Blues", annot=True)
# plt.title("Classification Matrix (number of reads)")
# plt.xlabel('Predicted Genome')
# plt.ylabel('Readset')
# plt.show()

if save_flag:
    classification_matrix_np = np.asarray(classification_matrix)
    np.save(save_dir + 'classification_matrix.npy', classification_matrix_np)

# Read Fractions
ax = sns.heatmap(classification_fraction, square = True, cmap="Blues")
# ax = sns.heatmap(classification_fraction, square = True, cmap="Blues", annot=True)
plt.title("Classification Matrix (Percentages)")
plt.xlabel('Predicted Genome')
plt.ylabel('Readset')

if save_flag:
    plt.savefig(save_dir + 'classification_matrix.png', dpi = 1200, bbox_inches='tight')
    plt.clf()
else:
    plt.show()

# 3. Matches per read (histograms)

all_read_matches = []
bins = np.arange(100) - 0.5

# a. Readset by Readset
for i in range(0, prefix_count):
    all_read_matches = all_read_matches + matches[i]
    # plt.hist(matches[i], bins)
    # plt.xlabel('Number of Matches with Source Genome')
    # plt.ylabel('Frequency')
    # plt.title("Readset " + str(i))
    # plt.show()


# b. All reads
all_read_matches_np = np.asarray(all_read_matches)
plt.hist(all_read_matches_np, bins)
plt.xlabel('Number of Matches with Source Genome')
plt.ylabel('Frequency')
plt.title("Read Matches (Target = " + str(target) + ", Error Rate = " + error_rate + ")")

# Mean Line
plt.axvline(all_read_matches_np.mean(), color='k', linestyle='dashed', linewidth=1)
min_ylim, max_ylim = plt.ylim()
plt.text(all_read_matches_np.mean()*1.1, max_ylim*0.9, 'Mean: {:.2f}'.format(all_read_matches_np.mean()))


if save_flag:
    plt.savefig(save_dir + 'read_matches.png', dpi = 1200, bbox_inches='tight')
    plt.clf()
else:
    plt.show()
