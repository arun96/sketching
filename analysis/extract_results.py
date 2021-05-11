import os
import sys
import numpy as np
import matplotlib.pyplot as plt

file_name = sys.argv[1]

with open(file_name) as f:
    results = f.readlines()

members = int((len(results)-2)/3)

results = [x.strip() for x in results]

sketches = results[1:members+1]
total_size = 0

for s in sketches:
    s_list = s.split()
    total_size += int(s_list[1])

sketch_results = results[2+members:]
result_lines = sketch_results[1::2]

total_reads = []
correct_reads = []
insuf_reads = []
mis_reads = []
tie_reads = []

accuracy = []

for r in result_lines:
    r_list = r.split()
    total_reads.append(int(r_list[0]))
    correct_reads.append(int(r_list[1]))
    mis_reads.append(int(r_list[2]))
    insuf_reads.append(int(r_list[3]))
    tie_reads.append(int(r_list[4]))
    a = int(r_list[1])/int(r_list[0])
    accuracy.append(a)

print("Total size of the community, Total Number of reads, Number of Correctly Classified Reads, Number of Misclassified Reads, Number of Reads with too few matches to classify, Number of reads with tied numbers of matches:")
print(total_size, sum(total_reads), sum(correct_reads), sum(mis_reads), sum(insuf_reads), sum(tie_reads))

# Comment out - this is specifically for 34 microbes + 22 chromosomes
# print(sum(correct_reads[:34]), sum(total_reads[:34]), sum(correct_reads[:34])/sum(total_reads[:34]))
# print(sum(correct_reads[:34]) - correct_reads[5], sum(total_reads[:34]) - total_reads[5], (sum(correct_reads[:34]) - correct_reads[5])/(sum(total_reads[:34]) - total_reads[5]))
# print(sum(correct_reads[34:]), sum(total_reads[34:]), sum(correct_reads[34:])/sum(total_reads[34:]))

print("Accuracy for each genome:")
print(accuracy)

plt.hist(accuracy, bins = 10)
plt.xlabel('Accuracy')
plt.ylabel('Frequency')
plt.title("Classification Accuracy across all genomes")
plt.show()
