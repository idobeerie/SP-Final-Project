#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double** create_adj_matrix(double** centroids, int n, int d) {   // we need to remmeber to free this and also create X to be the centroids as a matrix
    double** adj_matrix = (double*) malloc(n * sizeof(double));  
    double diff; 
    for (int i = 0; i < n; i++) {
        adj_matrix[i] = (double*) malloc(n * sizeof(double));
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

double** create_degree_matrix(double** adj_matrix, int n) {    // free this also!!!
    double** degree_matrix = (double*) malloc(n * sizeof(double));
    set_matrix_to_zero(degree_matrix, n);  
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


double** create_laplacian_matrix(double** D, double** W, int n) {     // free this also!!!or not i dont have a clue rami is a benzona
    double** laplacian_matrix = malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
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


int find_number_k(int n, double ** res_d){   // find the number of clusters
    int k;
    double max_diff = 0.0;
    double* eigenVals = (double*) malloc(n * sizeof(double));
    for(int i = 0; i < n; i++){
        eigenVals[i] = res_d[i][i];
    }
    mergeSort(eigenVals, 0, n - 1);
    for(int i = 1; i < n; i++){
        if(eigenVals[i] - eigenVals[i - 1] > max_diff){
            k = eigenVals[i] - eigenVals[i - 1];
            break;
        }
    }
    
}




void merge(int arr[], int l, int m, int r) // for sorting the eigenvalues
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
  
    /* create temp arrays */
    int L[n1], R[n2];
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
  
    /* Copy the remaining elements of L[], if there
    are any */
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }
  
    /* Copy the remaining elements of R[], if there
    are any */
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}
  
void mergeSort(int arr[], int l, int r) // for sorting the eigenvalues
{
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;
  
        // Sort first and second halves
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);
  
        merge(arr, l, m, r);
    }
}


void jacobi(double** L, int n){
    int iter = 0;
    double current_off, prev_off = 0.0, epsi = 1.0*0.00001;  //global maybe?
    current_off = off(L, n);
    double** rotation = create_p(L, n);
    double** eigenVectors;
    while(iter < 100 || prev_off - current_off > epsi){
        //add the changing L
        prev_off = current_off;
        if(iter != 0){
            eigenVectors = matrix_multiply(, , n);
        }
        rotation = create_p(L, n);
        current_off = off(L, n);
        iter++;
    }
    

}
    



