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

void jacobi_rotation(double** A, int n, double** V) {
    int i, j, k;
    double p, y, t, c, s;
    double* b = (double*) malloc(n * sizeof(double));
    double* z = (double*) malloc(n * sizeof(double));
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            V[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    while (1) {
        p = 0.0;
        for (i = 0; i < n-1; i++) {
            for (j = i+1; j < n; j++) {
                p += fabs(A[i][j]);
            }
        }
        if (p < TOLERANCE) {
            break;
        }
        
        for (i = 0; i < n-1; i++) {
            for (j = i+1; j < n; j++) {
                if (fabs(A[i][j]) > TOLERANCE) {
                    t = (A[j][j] - A[i][i]) / (2.0 * A[i][j]);
                    if (t >= 0.0) {
                        y = 1.0 / (t + sqrt(1.0 + t*t));
                    } else {
                        y = -1.0 / (-t + sqrt(1.0 + t*t));
                    }
                    c = 1.0 / sqrt(1.0 + y*y);
                    s = y * c;
                    for (k = 0; k < n; k++) {
                        b[k] = c * A[i][k] - s * A[j][k];
                        A[j][k] = s * A[i][k] + c * A[j][k];
                    }
                    for (k = 0; k < n; k++) {
                        A[i][k] = b[k];
                    }
                    for (k = 0; k < n; k++) {
                        b[k] = c * A[k][i] - s * A[k][j];
                        A[k][j] = s * A[k][i] + c * A[k][j];
                    }
                    for (k = 0; k < n; k++) {
                        A[k][i] = b[k];
                    }
                    for (k = 0; k < n; k++) {
                        b[k] = c * V[k][i] - s * V[k][j];
                        V[k][j] = s * V[k][i] + c * V[k][j];
                    }
                    for (k = 0; k < n; k++) {
                        V[k][i] = b[k];
                    }
                }
            }
        }
    }
    free(b);
    free(z);
}
    
double** create_eigenvectors(double** laplacian_matrix, int n) {
    double** eigenvectors = (double*) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        eigenvectors[i] = (double*) malloc(n * sizeof(double));
    }
    jacobi_rotation(laplacian_matrix, n, eigenvectors);
    return eigenvectors;
}


