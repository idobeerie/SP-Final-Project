import pandas as pd
import numpy as np
import sys
import math
import os
import random
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
np.random.seed(0)
eps = 0.0
k = 0



if len(sys.argv) < 3 or len(sys.argv) > 4:
    print("Error: missing arguments")
    sys.exit()

if len(sys.argv) == 4:
    k = int(sys.argv[1])
    goal = sys.argv[2]
    filename = sys.argv[3]
elif len(sys.argv) == 3:
    goal = sys.argv[1]
    filename = sys.argv[2]

centroids = []
min_index = 0
min_d = sys.float_info.max
current_d = 0
members_count = []
count = 0
temp_clusters = []
data = []



def create_centroids(data, k):   #the datais U
    centroids = []
    indices = []

    centroid_no1_idx = int(np.random.choice(data.shape[0], 1))
    centroid_no1 = data[centroid_no1_idx]
    centroids.append(centroid_no1)
    indices.append(centroid_no1_idx)
    min_idx = np.sqrt(np.sum((data - centroid_no1) ** 2, axis=1))
    for i in range(k - 1):
        for cent in  centroids:
            current_idx = np.sqrt(np.sum((data - cent) ** 2, axis=1))
            min_idx = np.minimum(min_idx, current_idx)
        
        p = min_idx / np.sum(min_idx)
        next_idx = np.random.choice(data.shape[0], 1, p=p)
        new_centroid = data[next_idx]
        centroids.append(new_centroid)
        indices.append(next_idx)
    
    return pd.DataFrame(centroids), indices

if goal =='spk':
    df = pd.read_csv(filename, header=None, dtype=float)
    X = df.values.tolist()
    L = km.gl(df.shape[0], df.shape[1], X)
    T = np.array(L)
    jacobi_vals = km.jacobi(df.shape[0], df.shape[0], T)
    eigen_vals = jacobi_vals[0].copy()
    eigen_vals = np.array(eigen_vals)
    eigen_vals = np.sort(eigen_vals)
    if k == 0:
        eigen_vals_cut = eigen_vals[:len(eigen_vals) /2]
        max_dif = sys.float_info.min
        for i in range(0,len(eigen_vals_cut)-1):
            if float(eigen_vals_cut[i+1] - eigen_vals_cut[i]) > max_dif:
                max_dif = float(eigen_vals_cut[i+1] - eigen_vals_cut[i])
        k = math.floor(max_dif)
    i = np.argsort(jacobi_vals[0])
    U = jacobi_vals[:,i]
    U = U[1:,:k+1]
    U_lst = U.tolist()
    start_indices, start_centroids = create_centroids(U, k)
    start_centroids = start_centroids.tolist()
    kmeans_result = km.spk(U_lst, start_centroids, 0.0)   

    # print(",".join(f'{x:.0f}' for x in start_indices))

    # for cluster in kmeans_result:
    #     print(",".join(f'{x:.0f}' for x in cluster))
    
    if goal in["gl", "ddg", "wam"]:
        df = pd.read_csv(filename, header=None, dtype=float)
        X = df.values.tolist()
        if "goal" == "gl":
            mat = km.gl(df.shape[0], df.shape[1], X)
        elif goal == "ddg":
            mat = km.ddg(df.shape[0], df.shape[1], X)
        elif goal == "wam":
            mat = km.wam(df.shape[0], df.shape[1], X)
        # for row in mat:
        #     print(",".join(f'{x:.4f}' for x in row))

if goal == "jacobi":
    sym_df = pd.read_csv(filename, header=None, dtype=float)
    sym_mat = sym_df.values.tolist()
    jacobi_res = km.jacobi(sym_mat, 0)
    jacobi_vals = jacobi_res[0]
    jacobi_vectors = jacobi_res[1:]
    print(",".join(f'{x:.4f}' for x in jacobi_vals))
    for row in jacobi_vectors:
        print(",".join(f'{x:.4f}' for x in row))
