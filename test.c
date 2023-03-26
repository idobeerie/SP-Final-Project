#include <stdio.h>
#include <stdlib.h>
#include <math.h>



double** wam(double** centroids, int n, int d) {   // we need to remmeber to free this and also create X to be the centroids as a matrix
    double** adj_matrix = malloc(n * sizeof(double*));  
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

double** ddg(double** adj_matrix, int n) {    
    double** degree_matrix = malloc(n * sizeof(double*));
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


double** gl(double** D, double** W, int n) {     // free this also
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

void matrix_multiply(double** A, double** B, int n){   // multiply two matrices
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
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = C[i][j];
        }
    }

    // Finally, free C and return A
    for (int i = 0; i < n; i++) {
        free(C[i]);
    }
    free(C);
    return;
}

double cmpfunc (const void * a, const void * b) {
   return (*(double*)a - *(double*)b);
}


 int find_number_of_k(double** jacobi_res, int n){   // the jacobi_res is n+1 on n
    double k;
    k = 0.0;
    for(int i = 1; i <= n / 2; i++){
        if(fabs(jacobi_res[0][i] - jacobi_res[0][i - 1]) > k){
            k = fabs(jacobi_res[0][i] - jacobi_res[0][i - 1]);
        }
    }
    return (int)k;
 }

void print_matrix(double** mat, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            printf("%f ", mat[i][j]);
        }
        printf("\n");
    }
}
// double** sort_jacobi(int n + 1, double ** res_d){   
//     double* temp = malloc((n + 1) * sizeof(double));
//     for (int i = 0; i < N; i++) {
//         for (int j = 0; j < n-i-1; j++) {
//             if (res_d[0][j] > res_d[0][j+1]) {
//                 // Swap eigenvalues
//                 temp[0] = res_d[0][j];
//                 res_d[0][j] = res_d[0][j+1];
//                 res_d[0][j+1] = temp[0];
//                 // Swap eigenvectors
//                 for (int k = 1; k < N+1; k++) {
//                     temp[k] = matrix[k][j];
//                     matrix[k][j] = matrix[k][j+1];
//                     matrix[k][j+1] = temp[k];
//                 }
//             }
//         }
//     }
//     return res_d;
// }

void free_matrix(double** mat, int n) {
    for (int i = 0; i < n; i++) {
        free(mat[i]);
    }
    free(mat);
}

double** jacobi(double** L, int n){
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
            matrix_multiply(eigenVectors, rotation, n);
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
    for(rows = 0; rows < n; rows++){
        free(rotation[rows]);
        free(eigenVectors[rows]);
    }
    free(rotation);
    free(eigenVectors);
    return eigenVecVal;
}
    
int main(int argc, char* argv){   // this will be in spkmeans.c
    FILE* fp;
    int d, int n, int sepCount;
    d = 0;
    double** laplacian;
    double** centroids;
    double** adj_matrix;
    double** degree_matrix;
    double** jacobi_res;
    char sepereator;
    char* goal;
    if(argc != 3){
        printf("not enough arguments, sory");
        return 1;
    }
    fp = fopen(argv[2], "r");
    if(fp == NULL){
        printf("file not found");
        return 1;
    }
    goal = argv[1];  
    while(sepereator = fgetc(fp) != EOF && sepereator != '\n'){  // maybe do \n differently for nova 
        if(sepereator == ','){
            d++;
        }
    }
    sepCount = d - 1;
    while (sepereator = fgetc(fp) != EOF) {
        if(sepereator == ','){
            sepCount++;
        }
    }
    n = sepCount / (d - 1);  // maybe not d-1 i dont know
    rewind(fp);
    centroids = malloc(n * sizeof(double*));
    for(int i = 0; i < n; i++){
        centroids[i] = malloc(d * sizeof(double));
        for(int j = 0; j < d; j++){
            if(fscanf(fp, "%lf", &centroids[i][j])){
                printf("error reading file");
                return 1;
            }
        }
    }
    fclose(fp);
    if(strcmp("jacobi", goal) == 0){
        jacobi_res = jacobi(centroids, n);
        print_matrix(jacobi_res, n+1, n);
        free_matrix(jacobi_res, n+1);
    }
    adj_matrix = wam(centroids, n, d);
    else if(strcmp("wam", goal) == 0){
        print_matrix(adj_matrix, n, n);
        free_matrix(adj_matrix, n); 
    }
    degree_matrix = ddg(adj_matrix, n);
    else if(strcmp("ddg", goal) == 0){
        print_matrix(degree_matrix, n, n);
        free_matrix(adj_matrix, n);
        free_matrix(degree_matrix, n);
    }
    else if(strcmp("gl", goal) == 0){
        laplacian = create_laplacian(adj_matrix, degree_matrix, n);
        print_matrix(laplacian, n, n);
        free_matrix(laplacian, n);
        free_matrix(adj_matrix, n);
        free_matrix(degree_matrix, n);
    }
    free_matrix(centroids, n);
    else{
        printf("invalid goal");
    }
    return 1;  
}


