import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob

file_path = sys.argv[1]

members = int(sys.argv[2])

classification_matrix = [[0 for x in range(members)] for y in range(members)]

for i in range(0, members):

    prefix = str(i) + '_'
    selectedFiles = glob.glob(file_path + prefix  + '*.log')
    # selectedFiles.sort()
    # print(selectedFiles)

    for filename in selectedFiles:
        with open(filename) as f:
            results_tmp = f.readlines()
        # print(results_tmp)
        results = [x.strip().split() for x in results_tmp][0]
        # print(results)
        source = int(results[0])
        pred = int(results[1])
        source_matches = int(results[2])
        pred_matches = int(results[3])
        if source_matches > pred_matches:
            classification_matrix[source][source]+=1
        else:
            classification_matrix[source][pred]+=1

print(classification_matrix)
