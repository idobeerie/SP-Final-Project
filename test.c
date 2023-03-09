#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double** create_adj_matrix(double** centroids, int n, int d) {
    double** adj_matrix =  malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        adj_matrix[i] = malloc(n * sizeof(double));
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
    double** degree_matrix = malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        degree_matrix[i] = malloc(n * sizeof(double));
        double degree = 0;
        for (int j = 0; j < n; j++) {
            if adj_matrix[i][j] >0{
                degree += 1;
            }
        }
        degree_matrix[i][i] = degree;
    }
    return degree_matrix;
}


double** create_laplacian_matrix(double** degree_matrix, double** adj_matrix, int n) {
    double** laplacian_matrix = malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        laplacian_matrix[i] = malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) {
            laplacian_matrix[i][j] = degree_matrix[i][j] - adj_matrix[i][j];
        }
    }
    return laplacian_matrix;
}

void jacobi_rotation(double** A, int n) {
    int i, j, k;
    double p, y, t, c, s;
    double** V = malloc(n * sizeof(double*));
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
    return V;
}


double** matrix_multiplication(double** A, double** B, int n) {
    double** C = malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        C[i] = malloc(n * sizeof(double));
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = 0;
            for (int k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

//1.2.1 B
void apply_rotation(double** A, double** V, int n) {
    int i, j, k;
    double** B = malloc(n * sizeof(double*));
    for (i = 0; i < n; i++) {
        B[i] = malloc(n * sizeof(double));
    }
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            B[i][j] = 0.0;
            for (k = 0; k < n; k++) {
                B[i][j] += V[k][i] * A[k][j];
            }
        }
    }
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            A[i][j] = 0.0;
            for (k = 0; k < n; k++) {
                A[i][j] += B[i][k] * V[k][j];
            }
        }
    }
    
    for (i = 0; i < n; i++) {
        free(B[i]);
    }
    free(B);
}

int is_diagonal(double** A, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j && A[i][j] != 0.0) {
                return 0; // Matrix is not diagonal
            }
        }
    }
    return 1; // Matrix is diagonal
}

