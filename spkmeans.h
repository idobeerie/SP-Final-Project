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

double** create2DfromJacobi(Jacobi_output* jacobi_output, int n);

double** allocateMatrix(size_t rows, size_t cols);   // i really dont remember writing this, dont think we need it 
	
double** wam(double** centroids, int n, int d);
double **kmeans(double **dots,double **centroids, int k, int d, int n, int max_iter);
void update_centroids(double **dots,const int *dots_location, double **new_centroids, int n, int d, int k);
double calc_distance(const double *dot, const double *centroid, int d);
int check_equals_2d_list(double **list1, double **list2, int row_num, int col_num);
int find_nearest_centroid(double *dot, double **centroids, int k, int d);
void print_2d_array(double **array, int row_num, int col_num);

void set_matrix_to_zero(double** matrix, int n, int m);
void print_list(double *, int);

double** ddg(double** adj_matrix, int n);
void free_contiguous_mat(double **mat);

double** gl(double** D, double** W, int n);
int* find_pivot(double** matrix, int n);

double* calculate_c_s(double ii, double jj, double ij);

double** create_p(double** A, int n);

double off(double** A, int n);

double** transpose(double** A, int n);
void matrix_multiply(double** A, double** B, int n);

int cmpfunc (const void * a, const void * b);

void print_matrix(double** mat, int n, int m);

void free_matrix(double** mat, int n);

Jacobi_output* jacobi(double** L, int n);

#endif