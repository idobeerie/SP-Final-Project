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
    print(result_indexes[:-2])
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
        # 1. get jacobi result
        #2. get k: if k = 0 then k = max_dif else k = k
        #3. sort by first row of jacobi (eigenvalues) and get the first k columns (each row here is a point for kmeans++)
        #3.1 create centroids for kmeans++ with the first k columns of jacobi
        #4 do kmeans++ with k and the first k columns of jacobi
        jacobi_res = km.jacobi(n, d, X, 0)
        eigen_vals = np.array(jacobi_res)
        #sorted_j = eigen_vals[:,eigen_vals[0, :].argsort()]
        eigenvaulues_arr = eigen_vals[0]
       
        eigenvector_arr = eigen_vals[1:]

        i = np.argsort(eigenvaulues_arr)
        eigenvector_arr =  eigenvector_arr[:,i]
        eigenvaulues_arr = np.sort(eigenvaulues_arr)
        max_dif = sys.float_info.min
        max_dif_idx = 0
        if k ==0:
            for i in range(0, int(len(eigenvaulues_arr) /2)+1):
                if eigenvaulues_arr[i+1] - eigenvaulues_arr[i] > max_dif:
                    max_dif = eigenvaulues_arr[i+1] - eigenvaulues_arr[i]
                    max_dif_idx = i
            k = math.floor(max_dif_idx) + 1
        eigen_vectors = eigenvector_arr[:, :k]
        ev_df = pd.DataFrame(eigen_vectors)
       
        starting_centroids_df, starting_centroids_indices = init_centroids(ev_df, k)
        starting_centroids = starting_centroids_df.values.tolist()
        # res_spk = km.spk(k, n, k, starting_centroids, eigen_vectors, 1)
        res_spk = km.spk(k,n,k,starting_centroids, eigen_vectors.tolist())
        # df = pd.DataFrame({'coordinates':[i for i in eigen_vectors]})
        # chosen = np.random.choice(df.index)
        # centroids = pd.DataFrame(df['coordinates'].filter(items = [chosen]))
        # scatters = df.drop(chosen,axis=0)
        # while (len(centroids) != k):
        #         probs = []
        #         for scatter in scatters["coordinates"]:
        #             def euc_dist(x): return np.linalg.norm(x - scatter)
        #             probs += [min(centroids["coordinates"].map(euc_dist))]
        #         probs = np.array(probs)
        #         probs = probs / sum(probs)
        #         chosen = (np.random.choice(scatters.index, p=probs))
        #         added_to_centroids = pd.DataFrame(
        #             scatters["coordinates"].filter(items=[chosen]))
        #         centroids = pd.concat([centroids, added_to_centroids])
        #         scatters = scatters.drop(chosen, axis=0)
        # dataPointList = [arr.tolist() for arr in scatters["coordinates"]]
        # centroidsList = [arr.tolist() for arr in centroids["coordinates"]]
        # dataPointList = dataPointList+centroidsList
        # Prints centroid's indices.
        # for i in centroids.index[0: -1]:
        #     print(i, end=",")
        # print(centroids.index[len(centroids.index) - 1])
       
        print(",".join(f'{x:.0f}' for x in starting_centroids_indices))

        for j in range(0, len(res_spk)):
            for i in range(0, len(res_spk[0])):
                res_spk[j][i] = round(res_spk[j][i], 4)

        for result in res_spk:
            print(','.join(str(x) for x in res_spk))
        
    if goal in["gl", "ddg", "wam"]:
        df = pd.read_csv(filename, header=None, dtype=float)
        X = df.values.tolist()

        if goal == "gl":
            mat = km.gl(df.shape[0], df.shape[1], X, 1)
        elif goal == "ddg":
            mat = km.ddg(df.shape[0], df.shape[1], X, 1)
        elif goal == "wam":
            mat = km.wam(df.shape[0], df.shape[1], X, 1)
        # for row in mat:
        #     print(",".join(f'{x:.4f}' for x in row))

    if goal == "jacobi":
        sym_df = pd.read_csv(filename, header=None, dtype=float)
        sym_mat = sym_df.values.tolist()
        jacobi_res = km.jacobi(sym_df.shape[0], sym_df.shape[0], sym_mat, 1)
        jacobi_vals = jacobi_res[0]
        jacobi_vectors = jacobi_res[1:]
        # print(",".join(f'{x:.4f}' for x in jacobi_vals))
        # for row in jacobi_vectors:
            # print(",".join(f'{x:.4f}' for x in row)
