#ifndef SPKMEANS
#define SPKMEANS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


/* structs: */
typedef struct ROTATION_MATRIX{
    double c;
    double s;
    int i;
    int j;
} rotation_mat;

typedef struct Jacobi_output{
    double *eigenValues;
    double **V;
} Jacobi_output;

typedef struct eigenStruct{
    double *Pointer;
    int originalIndex;  
} eigenStruct;
double** alloc_nXm_matrix(int n, int m);
double** create2DfromJacobi(Jacobi_output* jacobi_output, int n);
double** create_identity_matrix(int n);
double** allocateMatrix(size_t rows, size_t cols);   // i really dont remember writing this, dont think we need it 
double** alloc_nXn_matrix(int n);
double** wam(double** centroids, int n, int d);

int* find_off_diag_max_abs_val(double** A, int n);
double** kmeans(double** dots,double **centroids, int k, int d, int n, int max_iter);
void update_centroids(double** dots,const int* dots_location, double** new_centroids, int n, int d, int k);
double calc_distance(const double* dot, const double* centroid, int d);
int check_equals_2d_list(double** list1, double** list2, int row_num, int col_num);
int find_nearest_centroid(double* dot, double** centroids, int k, int d);
void print_2d_array(double** array, int row_num, int col_num);
rotation_mat* calc_rotation_mat(rotation_mat *P, double **A, int n);
double calc_of_f_square(double** A, int n);
void set_matrix_to_zero(double** matrix, int n, int m);
void print_list(double *array, int len);
double** calc_A_tag(double** A_tag, double** A, int n, rotation_mat* P);
double** ddg(double** adj_matrix, int n);
void free_contiguous_mat(double** mat);
double** getT(double **dots, int d, int n, int *k);
double** gl(double** D, double** W, int n);
Jacobi_output* get_jacobi(double** dots, int d, int n);
void multiply_rotation_matrix(double** V, rotation_mat* P, int n);
double** transpose(double** A, int n);
double **calc_T(Jacobi_output *jacobiOutput, int n, int *k_pointer);
double **gl_py(double **D, double **W, int n);
double** deepCopy2DArray(double **A, int row_num, int col_num);
void print_matrix(double** mat, int n, int m);

void free_matrix(double** mat);
double **calc_T(Jacobi_output *jacobiOutput, int n, int *k_pointer);
int compare_eigenStruct(const void *a, const void *b);
Jacobi_output* jacobi(double** L, int n);

#endif