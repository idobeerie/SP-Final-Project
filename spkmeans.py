import pandas as pd
import numpy as np
import sys
import kymeanssp as km

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


def create_centroids(data, k):
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
    df = pd.read_csv(filename, header=None)
    X = df.values.tolist()
    L = mk.gl(X)
    if k == 0:
        jacobi_vals = km.jacobi(L, 0)
        eigen_vals = jacobi_vals[0].copy()
        eigen_vals = np.array(eigen_vals)
        eigen_vals = np.sort(eigen_vals)
        eigen_vals_cut = eigen_vals[:len(eigen_vals) /2]
        max_dif = sys.float_info.min
        for i in range(0,len(eigen_vals_cut)-1):
            if float(eigen_vals_cut[i+1] - eigen_vals_cut[i]) > max_dif:
                max_dif = float(eigen_vals_cut[i+1] - eigen_vals_cut[i])
        k = math.floor(max_dif)
        i = np.argsort(jacobi_vals[0])
        U = jacobi_vals[:,i]
        U = U[1:,:k+1]
        start_centroids, start_indices = create_centroids(pd.DataFrame(U), k)
    starting_centroids = start_centroids.values.tolist()

    kmeans_result = km.spk(U, starting_centroids)

    print(",".join(f'{x:.0f}' for x in start_indices))

    for cluster in kmeans_result:
        print(",".join(f'{x:.0f}' for x in cluster))
    
    if goal in["gl", "ddg", "wam"]:
        df = pd.read_csv(filename, header=None)
        X = df.values.tolist()
        if "goal" == "gl":
            mat = km.gl(X)
        elif goal == "ddg":
            mat = km.ddg(X)
        elif goal == "wam":
            mat = km.wam(X)
        for row in mat:
            print(",".join(f'{x:.4f}' for x in row))

if goal == "jacobi":
    sym_df = pd.read_csv(filename, header=None)
    sym_mat = sym_df.values.tolist()
    jacobi_vals, jacobi_vectors, k = km.jacobi(sym_mat, 0)
    print(",".join(f'{x:.4f}' for x in jacobi_vals))
    for row in jacobi_vectors:
        print(",".join(f'{x:.4f}' for x in row))
