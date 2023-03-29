import sys
import math
import os
import random
import numpy as np
import pandas as pd
import mykmeanssp as km

def distance(arr1, arr2):
    d = 0
    for index in range(0,len(arr1)):
        d += math.pow((float(arr1[index]) - float(arr2[index])), 2)
    d = math.sqrt(d)
    return d

def selectNextCentroid(data, centroids, indexes):
    min_distances = []
    distances_from_centroids = []
    for dataPoint in data:
        for centroid in centroids:
            distances_from_centroids.append(distance(dataPoint, centroid))
        min_distances.append(min(distances_from_centroids))
        distances_from_centroids = []
    probabilities = []
    for i in indexes:
        probabilities.append(min_distances[indexes.index(i)] / sum(min_distances))
    chosen = np.random.choice(indexes, p=probabilities)
    return data[indexes.index(chosen)], chosen

def kmeansPlusPlus(data, k, indexes):
    result_indexes = ""
    centroids = []
    centroid = data[indexes.index(np.random.choice(indexes))]
    centroids.append(centroid)
    result_indexes += f"{indexes[data.index(centroid)]}, "
    for i in range(1, k):
        new_centroid, index = selectNextCentroid(data, centroids, indexes)
        centroids.append(new_centroid)
        result_indexes += f"{index}, "
    print(result_indexes[:-2])
    return centroids

if len(sys.argv) < 6:
    epsi = float(sys.argv[2])
    file_name1 = sys.argv[3]
    file_name2 = sys.argv[4]
    iter = 300
else:
    epsi = float(sys.argv[3])
    file_name1 = sys.argv[4]
    file_name2 = sys.argv[5]
    iter = int(sys.argv[2])

k = int(sys.argv[1])
centroids = []
min_index = 0
min_d = sys.float_info.max
current_d = 0
members_count = []
count = 0
temp_clusters = []
data = []

with open(f"{file_name1}", "r") as file:
    N = len(file.readlines())
    file.seek(0)

if k >= N or k < 1:
    print("Invalid number of clusters!")
    sys.exit()

if int(iter) < 1 or int(iter) >= 1000:
    print("Invalid maximum iteration!")
    sys.exit()
float_vals = []

np.random.seed(0)


with open(f"{file_name1}") as file:
    for line in file:
        float_vals = []
        current_values = line.split(",")
        current_values[-1] = current_values[-1].replace(os.linesep, "")
        current_values[-1] = current_values[-1].replace("\n", "")
        for value in current_values:
            float_vals.append(float(value))
        data.append(float_vals)

dataframe1 = pd.DataFrame(data, columns=range(0, len(data[0])))
data = []
with open(f"{file_name2}") as file:
    for line in file:
        float_vals = []
        current_values = line.split(",")
        current_values[-1] = current_values[-1].replace(os.linesep, "")
        current_values[-1] = current_values[-1].replace("\n", "")
        for value in current_values:
            float_vals.append(float(value))
        data.append(float_vals)

dataframe2 = pd.DataFrame(data, columns=range(0, len(data[0])))
dataframe1 = dataframe1.merge(dataframe2, how='inner', on=0)
dataframe1.rename(columns = {0: "0"}, inplace=True)
dataframe1["0"] = dataframe1["0"].astype('int')
dataframe1 = dataframe1.set_index("0")
dataframe1 = dataframe1.sort_index()
indexes = dataframe1.index.to_list()
data = dataframe1.values.tolist()
centroids = kmeansPlusPlus(data, k, indexes)
results = km.fit(centroids, data, epsi)

for j in range(0, len(results)):
    for i in range(0, len(results[0])):
        results[j][i] = round(results[j][i], 4)

for result in results:
    print(','.join(str(x) for x in result))
