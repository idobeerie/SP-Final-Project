#include "spkmeans.h"


    
int main(int argc, char* argv){   // this will be in spkmeans.c
    FILE* fp;
    int d, n;
    n = 0;
    d = 0;
    double** laplacian;
    double** centroids;
    double** adj_matrix;
    double** degree_matrix;
    double** jacobi_res;
    char sepereator;
    char* goal = malloc(10 * sizeof(char));
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
    sepereator = fgetc(fp);
    while(sepereator != '\n' || sepereator != EOF){
        if(sepereator == ','){
            d++;
        }
    }
    rewind(fp);
    for (sepereator = getc(fp); sepereator != EOF; sepereator = getc(fp)){
        if (sepereator == '\n'){ // Increment count if this character is newline
            n = n + 1;
        }
    }
    centroids = malloc(n * sizeof(double*));
    for(int i = 0; i < n; i++){
        centroids[i] = malloc(d * sizeof(double));
    }
    for(int i = 0; i < n; i++){
        for(int j = 0; j < d; j++){
            fscanf(fp, "%lf", &centroids[i][j]);
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


