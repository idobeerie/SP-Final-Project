#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double** create_adj_matrix(double** centroids, int n, int d) {
    double** adj_matrix = (double*) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        adj_matrix[i] = (double*) malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) {
            if (i == j) {
                adj_matrix[i][j] = 0;
            } else {
                double dist_squared = 0;
                for (int k = 0; k < d; k++) {
                    double diff = centroids[i][k] - centroids[j][k];
                    dist_squared += diff * diff;
                }
                adj_matrix[i][j] = exp(-dist_squared / 2);
            }
        }
    }
    return adj_matrix;
}

double** create_degree_matrix(double** adj_matrix, int n) {
    double** degree_matrix = (double*) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        degree_matrix[i] = (double*) malloc(n * sizeof(double));
        double degree = 0;
        for (int j = 0; j < n; j++) {
            degree += adj_matrix[i][j];
        }
        degree_matrix[i][i] = degree;
    }
    return degree_matrix;
}


double** create_laplacian_matrix(double** degree_matrix, double** adj_matrix, int n) {
    double** laplacian_matrix = (double*) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        laplacian_matrix[i] = (double*) malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) {
            laplacian_matrix[i][j] = degree_matrix[i][j] - adj_matrix[i][j];
        }
    }
    return laplacian_matrix;
}