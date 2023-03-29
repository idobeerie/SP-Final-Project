import sys
import math
import os

def distance(arr1, arr2):
    d = 0
    for index in range(0,len(arr1)):
        d += math.pow((float(arr1[index]) - float(arr2[index])), 2)
    d = math.sqrt(d)
    return d


if len(sys.argv) < 4:
    file_name = sys.argv[2]
    k = int(sys.argv[1])
    iter = 200
else:
    iter = int(sys.argv[2])
    file_name = sys.argv[3]
    k = int(sys.argv[1])

epsi = 0.001
centroids = []
min_index = 0
min_d = sys.float_info.max
current_d = 0
members_count = []
count = 0
temp_clusters = []
data = []

with open(file_name) as file:
    N = len(file.readlines())
    file.seek(0)

if k >= N or k < 1:
    print("Invalid number of clusters!")
    sys.exit()

if iter < 1 or iter >= 1000:
    print("Invalid maximum iteration!")
    sys.exit()
float_vals = []

with open(file_name) as file:
    for i in range(0, k):
        float_vals = []
        members_count.append(0)
        line = file.readline()
        current_values = line.split(",")
        current_values[-1] = current_values[-1].replace(os.linesep, "")
        current_values[-1] = current_values[-1].replace("\n", "")
        for value in current_values:
            float_vals.append(float(value))
        centroids.append(float_vals)
        temp_clusters.append(float_vals)
    file.seek(0)

with open(file_name) as file:
    for line in file:
        float_vals = []
        current_values = line.split(",")
        current_values[-1] = current_values[-1].replace(os.linesep, "")
        current_values[-1] = current_values[-1].replace("\n", "")
        for value in current_values:
            float_vals.append(float(value))
        data.append(float_vals)
    file.seek(0)

line_count = 0
small_change = True
choice =[]

while count < iter:
    for i in range(0, k):
        temp_clusters[i] = centroids[i].copy()
    for point in data:
        min_d = sys.float_info.max
        for index in range(0, k):
            current_d = distance(centroids[index], point)
            if current_d < min_d:
                min_index = index
                min_d = current_d
        choice.append(min_index)
    for cluster_idx in range(0, len(centroids)):
        for cent_idx in range(0, len(centroids[cluster_idx])):
            centroids[cluster_idx][cent_idx] = 0
        members_count[cluster_idx] = 0
    choice_count = 0
    for point in data:
        members_count[choice[choice_count]] += 1
        for i in range(0, len(centroids[choice[choice_count]])):
            centroids[choice[choice_count]][i] += float(point[i])
        choice_count+=1
    for cluster_idx in range(0, len(centroids)):
        for cent_idx in range(0, len(centroids[cluster_idx])):
            centroids[cluster_idx][cent_idx] = float(centroids[cluster_idx][cent_idx]) / members_count[cluster_idx]

    choice = []
    small_change = True
    for j in range(0, k):
        if distance(temp_clusters[j], centroids[j]) >= epsi:
            small_change = False
    if small_change:
        break
    count+=1
    if count >= iter:
        break


for j in range(0, len(centroids)):
    for i in range(0, len(centroids[0])):
        centroids[j][i] = round(centroids[j][i],4)

for centroid in centroids:
    print(centroid)
    