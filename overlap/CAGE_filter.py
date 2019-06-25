#!/usr/bin/env python3
# coding=utf-8

import pandas as pd
import datetime

start = datetime.datetime.now()

path = "/home/laverre/Documents/Regulatory_Landscape/data/human/CAGE/try_pandas/"

# Exclude all invalid sample ID
invalid_ID = []
with open(path+"CAGE_human_sample_info.csv", 'r') as f:
    for i in f.readlines()[1:]:
        i = i.strip('\n')
        i = i.split('\t')
        sample_ID = i[1]
        anomaly = i[3]
        if anomaly != "":
            invalid_ID.append(sample_ID)

sorted_ID = []
with open(path+"ordered_list.txt", 'r') as f:
    for i in f.readlines()[1:]:
        i = i.strip('\n')
        i = i.split('\t')
        if i[0] not in invalid_ID:
            sorted_ID.append(i[0])

print("Loading data...")
data = pd.read_csv(path+'human.enhancers.expression.tpm.matrix', sep="\t", index_col=0)

data = data.drop(invalid_ID, axis=1)
print("Drop invalid sample done ! ")

data = data.reindex(sorted_ID, axis=1)
print("Reorder samples done ! ")

# Select only enhancers with TPM >= 1 in at least 1 sample
def has_one_value_above_max(list_values, max_values):
    try:
        next(i for i in list_values if i > max_values)
        return True
    except StopIteration:
        return False


count = 0
for row in data.iterrows():
    index, value = row
    count += 1
    if count % 10000 == 0:
        print(count, 'sur 209912')

    if not has_one_value_above_max(value.tolist(), 1):
        data = data.drop(index, axis=0)

print("Drop invalid peaks done ! ")

print('Writting output...')
data.to_csv(path+'filtered_enhancer_matrix', sep='\t')

end = datetime.datetime.now()
print("Execution time:", str(end-start))
