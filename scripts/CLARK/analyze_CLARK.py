import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
from collections import Counter
import pandas as pd

CLARK_folder = sys.argv[1]
ext = ".csv"

CLARK_files = glob.glob(CLARK_folder + '*' + ext)
CLARK_files.sort()

results = []
for filename in CLARK_files:
    results_csv = pd.read_csv(filename, skipinitialspace = True)
    # print(list(results_csv))
    results.append(results_csv['Assignment'].to_list())

#TODO - handle N/A (nan) reads

# print(len(results))
for i in range(len(results)):
    # print(r)
    r = results[i]
    # print(CLARK_files[i])
    # print(r)
    # print([[x,r.count(x)] for x in set(r)])

# Save list of files, and assignments in each file
np.save("CLARK_files.npy", CLARK_files)
np.save("CLARK_results.npy", results)