#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double** wam(double** centroids, int n, int d) {   // we need to remmeber to free this and also create X to be the centroids as a matrix
    double** adj_matrix = (double*) malloc(n * sizeof(double));  
    double diff; 
    for (int i = 0; i < n; i++) {
        adj_matrix[i] = malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) {
            if (i == j) {
                adj_matrix[i][j] = 0;
            } else {
                double dist_squared = 0;
                for (int k = 0; k < d; k++) {
                    diff = centroids[i][k] - centroids[j][k];
                    dist_squared += diff * diff;
                }
                adj_matrix[i][j] = exp(-dist_squared / 2);
            }
        }
    }
    return adj_matrix;
}

void set_matrix_to_zero(double** matrix, int n) {    //helper function
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = 0;
        }
    }
}

double** ddg(double** adj_matrix, int n) {    // free this also!!!
    double** degree_matrix = (double*) malloc(n * sizeof(double));
    set_matrix_to_zero(degree_matrix, n);  
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


double** gl(double** D, double** W, int n) {     // free this also!!!or not i dont have a clue rami is a benzona
    double** laplacian_matrix = malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        laplacian_matrix[i] = malloc(n * sizeof(double));
        laplacian_matrix[i] = malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) {
            laplacian_matrix[i][j] = D[i][j] - W[i][j];
        }
    }
    return laplacian_matrix;
}


int* find_pivot(double** matrix, int n) {    // free this
    double max_abs = 0.0;
    int pivot[2] = {0, 0};   
    int max_i = 0, max_j = 0;
    for (int i = 0; i < n; i++) {   
        for (int j = i; j < n; j++) {
            if (fabs(matrix[i][j]) > max_abs) {
                max_abs = fabs(matrix[i][j]);
                max_i = i;
                max_j = j;
            }
        }
    }
    pivot[0] = max_i;
    pivot[1] = max_j;
    return pivot;
}

double* calculate_c_s(double ii, double jj, double ij) {   // free this
    double* c_s = malloc(2 * sizeof(double));
    double teta, c, s, t;
    teta = (jj - ii) / (2 * ij);
    if (teta >= 0) {
        t = 1 / (teta + sqrt(teta * teta + 1));
    } else {
        t = -1 / (-teta + sqrt(teta * teta + 1));
    }
    c = 1 / sqrt(t * t + 1);
    s = c * t;
    c_s[0] = c;
    c_s[1] = s;
    return c_s;
}

double** create_p(double** A, int n){   // create the rotation matrix, we need to see how we want to use it and free it, maybe we want to put this in the main function
    double teta, c, s, t;
    double** p = malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        p[i] = (double*) malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) {
            if (i == j) {
                p[i][j] = 1;
            } else {
                p[i][j] = 0;
            }
        }
    }
    int* pivot = find_pivot(A, n);
    double* c_s = calculate_c_s(A[pivot[0]][pivot[0]], A[pivot[1]][pivot[1]], A[pivot[0]][pivot[1]]);
    p[pivot[0]][pivot[0]] = c_s[0];
    p[pivot[0]][pivot[1]] = c_s[1];
    p[pivot[1]][pivot[0]] = -c_s[1];
    p[pivot[1]][pivot[1]] = c_s[0];
    return p;   
}


double off(double** A, int n){   // calculate the off diagonal sum of the matrix
    double sum = 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            if(i != j){
                sum += A[i][j] * A[i][j];
            }
        }
    }
    return sum;
}

double** transpose(double** A, int n){   // transpose a matrix, free this also
    double** A_t = (double*) malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        A_t[i] = (double*) malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) {
            A_t[i][j] = A[j][i];
        }
    }
    return A_t;
}

double** matrix_multiply(double** A, double** B, int n){   // multiply two matrices
    double** C = (double*) malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        C[i] = (double*) malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) {
            C[i][j] = 0;
            for (int k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

double cmpfunc (const void * a, const void * b) {
   return (*(double*)a - *(double*)b);
}

int find_number_k(int n + 1, double ** res_d){   // find the number of clusters
    double k;    
    k = 0.0;   //maybe nir is right and this is not allowed
    double* eigenVals = (double*) malloc(n * sizeof(double));
    for(int i = 0; i < n; i++){
        eigenVals[i] = res_d[0][i];
    }
    qsort(arr, n, sizeof(double), cmpfunc);
    for(int i = 1; i <= n / 2; i++){
        if(fabs(eigenVals[i] - eigenVals[i - 1]) > k){
            k = fabs(eigenVals[i] - eigenVals[i - 1]);
        }
    }
    // free eigenVals
    return (int)k;
}







void jacobi(double** L, int n){
    int iter = 0, i, j, rows, cols;
    double c, s;
    double current_off, prev_off = 0.0, epsi = 1.0*0.00001;  //global maybe?
    current_off = off(L, n);
    double** prev_L = L;
    double** rotation;
    double** eigenVectors = malloc(n * sizeof(double*)); // just the vectors
    double** eigenVecVal = malloc((n+1) * sizeof(double*)); // the first row will be the eigen values
    while(iter < 100 || prev_off - current_off > epsi){
        j = find_pivot(L, n)[1];
        i = find_pivot(L, n)[0];
        rotation = create_p(L, n);
        if(iter == 0){
            eigenVectors = rotation;
        }
        else{
            eigenVectors = matrix_multiply(eigenVectors, rotation, n);
            }
        c = rotation[i][i];
        s = rotation[i][j];
        for(int k=0; k<n;k++){
            L[i][k] = c*prev_L[i][k] - s*prev_L[j][k];
            L[j][k] = s*prev_L[i][k] + c*prev_L[j][k];
            L[k][i] = c*prev_L[k][i] - s*prev_L[k][j];
            L[k][j] = s*prev_L[k][i] + c*prev_L[k][j];
        }
        L[i][i] = c*c*prev_L[i][i] - 2*c*s*prev_L[i][j] + s*s*prev_L[j][j];
        L[j][j] = s*s*prev_L[i][i] + 2*c*s*prev_L[i][j] + c*c*prev_L[j][j];
        L[i][j] = (c*c-s*s)*L[i][j] + c*s*(prev_L[i][i] - prev_L[j][j]);
        L[j][i] = L[i][j];
        iter++;
        prev_off = current_off;
        current_off = off(L, n);
        prev_L = L;
        //free(rotation);   here free the arrays ido
    }

    for(rows = 0; rows < n + 1; rows++){
        eigenVecVal[rows] = (double*)malloc(n * sizeof(double));
        for(cols = 0; cols < n; cols++){
            if(rows == 0){
                eigenVecVal[rows][cols] = L[cols][cols];
            }
            else{
                eigenVecVal[rows][cols] = L[rows][cols];
            }
        }
    }
    // free the eigenVectors matrix ido
    // free the L matrix ido not sure though
}
    



