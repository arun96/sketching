import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
from collections import Counter
import pandas as pd

centrifuge_folder = sys.argv[1]
ext = ".log"

centrifuge_files = glob.glob(centrifuge_folder + '*' + ext)
centrifuge_files.sort()
results = []
for filename in centrifuge_files:
    results_csv = pd.read_csv(filename, skipinitialspace = True, sep='\t')
    # print(list(results_csv))
    results.append(results_csv['taxID'].to_list())

#TODO - handle N/A (nan) reads

# print(len(results))
for i in range(len(results)):
    # print(r)
    r = results[i]
    print(centrifuge_files[i])
    print([[x,r.count(x)] for x in set(r)])
