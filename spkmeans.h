#ifndef SPKMEANS
#define SPKMEANS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


double** allocateMatrix(size_t rows, size_t cols);   // i really dont remember writing this, dont think we need it 
	
double** wam(double** centroids, int n, int d);

void set_matrix_to_zero(double** matrix, int n) ;

double** ddg(double** adj_matrix, int n) ;

double** gl(double** D, double** W, int n);
int* find_pivot(double** matrix, int n) ;

double* calculate_c_s(double ii, double jj, double ij);

double** create_p(double** A, int n);

double off(double** A, int n);

double** transpose(double** A, int n);
void matrix_multiply(double** A, double** B, int n);

int cmpfunc (const void * a, const void * b);

void print_matrix(double** mat, int n, int m) ;

void free_matrix(double** mat, int n);

double** jacobi(double** L, int n);

#endif