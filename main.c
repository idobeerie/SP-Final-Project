#include "spkmeans.h"


int main(int argc, char** argv){   
    FILE* fp;
    int d, n;
    n = 0;
    d = 0;
    double** laplacian;
    double** centroids;
    double** adj_matrix;
    double** degree_matrix;
    Jacobi_output* jacobi_res;
    char sepereator;
    char* goal;
    if(argc != 3){
        printf("An Error Has Occurred");
        return 1;
    }
    fp = fopen(argv[2], "r");
    if(fp == NULL){
        printf("An Error Has Occurred");
        return 1;
    }
    goal = argv[1];  
    sepereator = fgetc(fp);
    while(sepereator != '\n' && sepereator != EOF){
        if(sepereator == ','){
            d++;
        }
        sepereator = fgetc(fp);
    }
    d++;
    rewind(fp);
    for (sepereator = getc(fp); sepereator != EOF; sepereator = getc(fp)){
        if (sepereator == '\n'){ // Increment count if this character is newline
            n = n + 1;
        }
    }
    centroids = allocateMatrix(n, d);
    if(centroids == NULL){
        printf("An Error Has Occurred");
        exit(-1);
    }
    rewind(fp);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < d; j++){
            fscanf(fp, "%lf,", &centroids[i][j]);
        
        }
    }
    fclose(fp);

    if(strcmp("jacobi", goal) == 0){
        jacobi_res = jacobi(centroids, n);
        print_list(jacobi_res->eigenValues, n);
        print_matrix(jacobi_res->V, n, n);
        free_contiguous_mat(jacobi_res->V);
        free(jacobi_res->eigenValues);
        free(jacobi_res);
    }
    adj_matrix = wam(centroids, n, d);
    degree_matrix = ddg(adj_matrix, n);
    if(strcmp("wam", goal) == 0){
        print_matrix(adj_matrix, n, n);
    }
    else if(strcmp("ddg", goal) == 0){
        print_matrix(degree_matrix, n, n);

    }
    else if(strcmp("gl", goal) == 0){
        laplacian = gl(adj_matrix, degree_matrix, n);
        print_matrix(laplacian, n, n);
        free_contiguous_mat(laplacian);
    }
    else if(strcmp("jacobi", goal) == 1){
        printf("An Error Has Occurred");
    }
    if(strcmp("gl", goal) != 0){
        free_contiguous_mat(adj_matrix);
        free_contiguous_mat(degree_matrix);
    }   
    free_contiguous_mat(centroids);
    return 1;  
}


