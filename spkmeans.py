import pandas as pd
import numpy as np
import sys
import math
import os
import random
import mykmeanssp as km

###
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
   #print(result_indexes[:-2])
    return centroids





# def create_centroids(data, k):   #the datais U
#         centroids = []
#         indices = []
#         centroid_no1_idx = int(np.random.choice(data.shape[0], 1))
#         centroid_no1 = data[centroid_no1_idx]
#         centroids.append(centroid_no1)
#         indices.append(centroid_no1_idx)
#         min_idx = np.sqrt(np.sum((data - centroid_no1) ** 2, axis=1))
#         for i in range(k - 1):
#             for cent in  centroids:
#                 current_idx = np.sqrt(np.sum((data - cent) ** 2, axis=1))
#                 min_idx = np.minimum(min_idx, current_idx)
            
#             p = min_idx / np.sum(min_idx)
#             next_idx = np.random.choice(data.shape[0], 1, p=p)
#             new_centroid = data[next_idx]
#             centroids.append(new_centroid)
#             indices.append(next_idx)
        
#         return pd.DataFrame(centroids), indices

# def init_centroids(data, num_clusters):
#     centroids = []
#     centroids_indices = []

#     first_centroid_index = int(np.random.choice(data.shape[0]))
#     first_centroid = data.iloc[first_centroid_index]
#     centroids.append(first_centroid)
#     centroids_indices.append(first_centroid_index)

#     min_dx = np.sqrt(np.sum((data - first_centroid) ** 2, axis=1))

#     for i in range(num_clusters - 1):  # init num_clusters-1 more centroids
#         for c in centroids:
#             current_dx = np.sqrt(np.sum((data - c) ** 2, axis=1))
#             min_dx = np.minimum(current_dx, min_dx)

#         prob = min_dx / np.sum(min_dx)

#         next_centroid_index = int(np.random.choice(data.shape[0], p=prob))
#         new_centroid = data.iloc[next_centroid_index]
#         centroids.append(new_centroid)
#         centroids_indices.append(next_centroid_index)

#     return pd.DataFrame(centroids), centroids_indices



def init_centroids(data, num_clusters):
    centroids = []
    centroids_indices = []

    first_centroid_index = int(np.random.choice(data.shape[0]))
    first_centroid = data.iloc[first_centroid_index]
    centroids.append(first_centroid)
    centroids_indices.append(first_centroid_index)

    min_dx = np.sqrt(np.sum((data - first_centroid) ** 2, axis=1))

    for i in range(num_clusters - 1):  # init num_clusters-1 more centroids
        for c in centroids:
            current_dx = np.sqrt(np.sum((data - c) ** 2, axis=1))
            min_dx = np.minimum(current_dx, min_dx)

        prob = min_dx / np.sum(min_dx)

        next_centroid_index = int(np.random.choice(data.shape[0], p=prob))
        new_centroid = data.iloc[next_centroid_index]
        centroids.append(new_centroid)
        centroids_indices.append(next_centroid_index)
    return pd.DataFrame(centroids), centroids_indices
    


if __name__ == '__main__':

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



    
    if goal =='spk':
        dff = pd.read_csv(filename, header=None, dtype=float)
        X = dff.values.tolist()
        n = dff.shape[0]
        d = dff.shape[1]
        gl_res = km.gl(n, d, X, 0)
        jacobi_res = km.jacobi(n, n, gl_res, 0)
        eigen_vals = np.array(jacobi_res)
        eigenvaulues_arr = eigen_vals[0]
        eigenvector_arr = eigen_vals[1:]

        i = np.argsort(eigenvaulues_arr)
        eigenvector_arr =  eigenvector_arr[:,i]
        eigenvaulues_arr = np.sort(eigenvaulues_arr)
        max_dif = sys.float_info.min
        max_dif_idx = 0
        if k ==0:
            for i in range(0, int(len(eigenvaulues_arr) /2)+1):
                if math.fabs(eigenvaulues_arr[i] - eigenvaulues_arr[i+1]) > max_dif:
                    max_dif = math.fabs(eigenvaulues_arr[i] - eigenvaulues_arr[i+1])
                    max_dif_idx = i+1
            k = math.floor(max_dif_idx) 
        eigen_vectors = eigenvector_arr[:, :k]
        ev_df = pd.DataFrame(eigen_vectors)
        starting_centroids_df, starting_centroids_indices = init_centroids(ev_df, k)
        starting_centroids = starting_centroids_df.values.tolist()
        print(",".join(f'{x:.0f}' for x in starting_centroids_indices))
        res_spk = km.spk(eigen_vectors.tolist(), starting_centroids, 300, 0.0)
       
        # for j in range(0, len(res_spk)):
        #     for i in range(0, len(res_spk[0])):
        #         res_spk[j][i] = round(res_spk[j][i], 4)

        for result in res_spk:
            result = [round(num,4) for num in result]
            print(','.join(map(str, result)))
        
    if goal in["gl", "ddg", "wam"]:
        df = pd.read_csv(filename, header=None, dtype=float)
        X = df.values.tolist()
        if goal == "gl":
            mat = km.gl(df.shape[0], df.shape[1], X, 0)
        elif goal == "ddg":
            mat = km.ddg(df.shape[0], df.shape[1], X, 0)
        elif goal == "wam":
            mat = km.wam(df.shape[0], df.shape[1], X, 0)
        for row in mat:
            print(",".join(f'{x:.4f}' for x in row))

    if goal == "jacobi":
        sym_df = pd.read_csv(filename, header=None, dtype=float)
        sym_mat = sym_df.values.tolist()
        jacobi_res = km.jacobi(sym_df.shape[0], sym_df.shape[1], sym_mat, 1)
        jacobi_vals = jacobi_res[0]
        jacobi_vectors = jacobi_res[1:]
        # print(",".join(f'{x:.4f}' for x in jacobi_vals))
        # for row in jacobi_vectors:
        #     print(",".join(f'{x:.4f}' for x in row))
