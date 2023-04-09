#include "spkmeans.h"

#define MAX_JACOBI_ITERS 100
#define verifyNotNULL(var)              \
    if ((var) == NULL)                  \
    {                                   \
        printf("An Error Has Occured"); \
        exit(-1);                       \
    }
#define EPSILON 0.00001

double **allocateMatrix(size_t rows, size_t cols)
{ // i really dont remember writing this, dont think we need it
    double *data;
    double **matrix;
    size_t i;
    data = calloc(rows * cols, sizeof(double));
    matrix = calloc(rows, sizeof(double *));
    verifyNotNULL(data);
    verifyNotNULL(matrix);
    for (i = 0; i < rows; i++)
    {
        matrix[i] = data + i * cols;
    }
    set_matrix_to_zero(matrix, rows, cols);
    return matrix;
}

double **create2DfromJacobi(Jacobi_output *jacobi_output, int n)
{
    double **output = allocateMatrix(n+1, n);
    int j=0,i=0;
    output[0] = jacobi_output->eigenValues;
    for ( i= 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            output[i+1][j] = jacobi_output->V[i][j];
        }
    }
    return output;
}
double **wam(double **centroids, int n, int d)
{
    double **adj_matrix = allocateMatrix(n, n);
    double diff;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                adj_matrix[i][j] = 0;
            }
            else
            {
                double dist_squared = 0;
                for (int k = 0; k < d; k++)
                {
                    diff = centroids[i][k] - centroids[j][k];
                    dist_squared += diff * diff;
                }
                adj_matrix[i][j] = exp(-dist_squared / 2);
            }
        }
    }
    return adj_matrix;
}

void set_matrix_to_zero(double **matrix, int n, int m)
{ // helper function
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            matrix[i][j] = 0;
        }
    }
}

double **ddg(double **adj_matrix, int n)
{
    double **degree_matrix = allocateMatrix(n, n);
    double degree = 0;
    // print_matrix(adj_matrix, n, n);
    // set_matrix_to_zero(degree_matrix, n, n);
    for (int i = 0; i < n; i++)
    {
        degree = 0;
        for (int j = 0; j < n; j++)
        {
            if (adj_matrix[i][j] > 0)
            {
                degree += adj_matrix[i][j];
            }
        }
        degree_matrix[i][i] = degree;
    }
    return degree_matrix;
}

double **gl(double **D, double **W, int n)
{
    double **laplacian_matrix = allocateMatrix(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            laplacian_matrix[i][j] = (D[i][j] - W[i][j]) * (-1);
        }
    }
    free_contiguous_mat(D);
    free_contiguous_mat(W);
    return laplacian_matrix;
}

double **gl_py(double **D, double **W, int n)
{
    double **laplacian_matrix = allocateMatrix(n, n);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            laplacian_matrix[i][j] = (D[i][j] - W[i][j]);
        }
    }
    free_contiguous_mat(D);
    free_contiguous_mat(W);
    return laplacian_matrix;
}




double **alloc_nXm_matrix(int n, int m)
{
    int i;
    double **output_matrix = (double **)calloc(n, sizeof(double *));
    double *matrix_data = (double *)calloc(n * m, sizeof(double));

    verifyNotNULL(output_matrix)
        verifyNotNULL(matrix_data)

            for (i = 0; i < n; i++)
                output_matrix[i] = matrix_data + i * m;

    return output_matrix;
}

/*------------------------- allocation memory for a nXn matrix: ----------------------------*/
double **alloc_nXn_matrix(int n)
{
    return alloc_nXm_matrix(n, n);
}

/*------------------------- free contiguous memory of a nXm matrix: ---------------------------*/
void free_contiguous_mat(double **mat)
{
    free(mat[0]);
    free(mat);
}

/*------------------------------- find_nearest_centroid: -----------------------------------------*/
/* returns the index of the nearest_centroid of the dot*/
int find_nearest_centroid(double *dot, double **centroids, int k, int d){
    double min_distance = calc_distance(dot, centroids[0], d);
    double distance;
    int min_index = 0, i;
    for (i = 1 ; i < k ; ++i){
        distance = calc_distance(dot, *(centroids+i), d);
        if (distance < min_distance){
            min_distance = distance;
            min_index = i;
        }
    }
    return min_index;
}


void print_matrix(double **array, int n, int m)
{
    int i, j;
    for(i=0; i < n; i++)
    {
        for(j = 0; j < m; j++) {
            printf("%.4f", array[i][j]);
            if(j!= m-1)
                printf(",");
        }
        printf("\n");
    }
}

/*----------------------------- creating Identity matrix n: --------------------------------------*/
double **create_identity_matrix(int n)
{
    double **I = alloc_nXn_matrix(n);
    int i, j;

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            I[i][j] = i == j ? 1 : 0;
    return I;
}

void free_matrix(double **mat)
{

    free(mat[0]);
    free(mat);
}

void copy_matrix(double **mat, double **copy, int n, int m)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            copy[i][j] = mat[i][j];
        }
    }
}

/* ----------------------------- deepCopy2DArray: --------------------------------------------*/

double **deepCopy2DArray(double **A, int row_num, int col_num)
{
    int i, j;
    double **B = (double **)alloc_nXm_matrix(row_num, col_num);
    for (i = 0; i < row_num; i++)
        for (j = 0; j < col_num; j++)
            B[i][j] = A[i][j];

    return B;
}

double calc_distance(const double *dot, const double *centroid, int d){
    double sum = 0;
    int i;
    for(i=0 ; i<d ; i++)
        sum += ((dot[i]-centroid[i])*(dot[i]-centroid[i]));
    return sum;
}

void multiply_rotation_matrix(double **V, rotation_mat *P, int n)
{
    int r;
    double res_i, res_j;
    for (r = 0; r < n; r++)
    {
        /* Note: using temp vars to avoid sync issues */
        res_i = V[r][P->i] * P->c + V[r][P->j] * P->s * -1;
        res_j = V[r][P->i] * P->s + V[r][P->j] * P->c;
        V[r][P->i] = res_i;
        V[r][P->j] = res_j;
    }
}

/* ----------------------- calc_rotation_mat: ---------------------------------*/

rotation_mat *calc_rotation_mat(rotation_mat *P, double **A, int n)
{
    /* returns a struct that represent the rotation matrix of A. */
    int sign_theta, i, j;
    double Aii, Ajj, Aij, theta, t, max_val = 0, abs_curr;

    /* TODO, tests it
        int *ij = find_off_diag_max_abs_val(A, n);
        P->i = ij[0];
        P->j = ij[1];
    */
    /*if all of-diagonal elements are zeros: */
    P->i = 0;
    P->j = 1;

    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            abs_curr = fabs(A[i][j]);
            if (abs_curr > max_val)
            {
                P->i = i;
                P->j = j;
                max_val = abs_curr;
            }
        }
    }

    Aii = A[P->i][P->i];
    Ajj = A[P->j][P->j];
    Aij = A[P->i][P->j];
    if (Aij == 0)
    {
        P->c = 1;
        P->s = 0;
        return P;
    }
    theta = (double)(Ajj - Aii) / (2 * Aij);
    sign_theta = (theta == 0) ? 1 : (int)(theta / fabs(theta));
    t = (double)(sign_theta / ((fabs(theta) + sqrt((theta * theta) + 1))));

    P->c = (double)(1 / (sqrt(t * t + 1)));
    P->s = (double)(t * (P->c));
    /* todo
        free(ij);
    */
    return P;
}

double calc_of_f_square(double **A, int n)
{
    double diagonal_elementwise_sum_square = 0, matrix_elementwise_sum_square = 0, f;
    int i, j;

    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            matrix_elementwise_sum_square += pow(A[i][j], 2);

    for (i = 0; i < n; i++)
        diagonal_elementwise_sum_square += pow(A[i][i], 2);

    f = sqrt(matrix_elementwise_sum_square);

    return f - diagonal_elementwise_sum_square;
}

/* ----------------------- Calculate A': ---------------------------------*/

/* Calculate A' by the guidance in Section 6 - Relation between A and A': */
double **calc_A_tag(double **A_tag, double **A, int n, rotation_mat *P)
{
    int r, col;

    /* Copying A to A', but for each row, at indices i,j assign the formula as described in Section 6. */
    for (r = 0; r < n; r++)
    {
        for (col = 0; col < n; col++)
        {
            if (col == P->i)
            {
                A_tag[r][col] = P->c * A[r][P->i] - P->s * A[r][P->j];
                A_tag[col][r] = A_tag[r][col];
            }
            else if (col == P->j)
            {
                A_tag[r][col] = P->c * A[r][P->j] + P->s * A[r][P->i];
                A_tag[col][r] = A_tag[r][col];
            }
            else if (r != P->j && r !=
                                      P->i) /* To prevent override bc this is a symmetric matrix - and in the above if we are doing it also to the opposite indices TODO - rewrite this comment. */
                A_tag[r][col] = A[r][col];
        }
    }
    /* Overriding A_tag indices for the last 3 formulas in Section 6 (Indices ii,jj,ij,ji) */
    A_tag[P->i][P->i] = pow(P->c, 2) * A[P->i][P->i] + pow(P->s, 2) * A[P->j][P->j] - 2 * P->s * P->c * A[P->i][P->j];
    A_tag[P->j][P->j] = pow(P->s, 2) * A[P->i][P->i] + pow(P->c, 2) * A[P->j][P->j] + 2 * P->s * P->c * A[P->i][P->j];
    // A_tag[P->i][P->j] = ((pow(P->c, 2) - pow(P->s, 2)) * A[P->i][P->j]) + (P->s * P->c * (A[P->i][P->i] - A[P->j][P->j]));
    // A_tag[P->j][P->i] = A_tag[P->i][P->j];
    A_tag[P->i][P->j] = 0;    //not sure if this will solve the problem
    A_tag[P->j][P->i] = 0;
    return A_tag;
}



/*//////////////////////////////////////// Jacobi Algorithm: //////////////////////////////////////////////////////*/
/* ----------------------- Jacobi main function: ---------------------------------*/

Jacobi_output *jacobi(double **A, int n)
{
    double **A_tag, **V, **temp;
    double *eigenValues;
    Jacobi_output *output_struct;
    rotation_mat *P;
    int i = 1;

    A = deepCopy2DArray(A, n, n);
    V = create_identity_matrix(n);

    P = (rotation_mat *)malloc(sizeof(rotation_mat));
    verifyNotNULL(P)
        calc_rotation_mat(P, A, n);

    A_tag = calc_A_tag(alloc_nXn_matrix(n), A, n, P);

    multiply_rotation_matrix(V, P, n);
    while (calc_of_f_square(A, n) - calc_of_f_square(A_tag, n) > EPSILON && i <= MAX_JACOBI_ITERS)
    {
        temp = A;
        A = A_tag;
        calc_rotation_mat(P, A, n);
        A_tag = calc_A_tag(temp, A, n, P);
        multiply_rotation_matrix(V, P, n);
        i++;
    }

    output_struct = (Jacobi_output *)malloc(sizeof(Jacobi_output));
    eigenValues = (double *)calloc(n, sizeof(double));
    verifyNotNULL(output_struct)
        verifyNotNULL(eigenValues)

            for (i = 0; i < n; i++)
                eigenValues[i] = A_tag[i][i];

    output_struct->eigenValues = eigenValues;
    output_struct->V = V;

    free(P);
    free_contiguous_mat(A_tag);
    free_contiguous_mat(A);

    return output_struct;
}






int *find_off_diag_max_abs_val(double **A, int n)
{
    /*A is symmetric. n is the dimensions of A. Returns a list of [i,j]. */
    int i, j;
    double max_val = 0, abs_curr;
    int *ij_list = (int *)calloc(2, sizeof(int));
    verifyNotNULL(ij_list)

        /*if all of-diagonal elements are zeros: */
        ij_list[0] = 0;
    ij_list[1] = 1;

    for (i = 0; i < n; i++)
    {
        for (j = i + 1; j < n; j++)
        {
            abs_curr = fabs(A[i][j]);
            if (abs_curr > max_val)
            {
                ij_list[0] = i;
                ij_list[1] = j;
                max_val = abs_curr;
            }
        }
    }
    return ij_list;
}



void print_list(double *array, int len)
{
    int i;
    double num;
    for (i = 0; i < len; i++)
    {
        num = array[i];
        if (num > -0.00005 && num < 0)
            printf("%.4f", 0.0);
        else
            printf("%.4f", array[i]);
        if (i != len - 1)
            printf(",");
    }
    printf("\n");
}


/*compare between two 2 dimensions lists of the same dimensions:*/
int check_equals_2d_list(double **list1, double **list2, int row_num, int col_num){
    int i, j;
    for(i = 0; i < row_num; i++) {
        for (j = 0; j < col_num; j++) {
            if (list1[i][j] != list2[i][j])
                return 0;
        }
    }
    return 1;
}
/*------------------------------- update_centroid: -----------------------------------------*/

void update_centroids(double **dots,const int *dots_location, double **centroids, int n, int d, int k) {

    /* space allocation */
    int i, j;
    int *num_of_dots_per_cluster = (int*) malloc(k * sizeof(int));
    double **sum_of_dots_per_cluster = (double **) malloc(k*sizeof(double *));
    verifyNotNULL(num_of_dots_per_cluster)
    verifyNotNULL(sum_of_dots_per_cluster)


    for (i = 0; i < k; i++) {
        sum_of_dots_per_cluster[i] = (double *) calloc(d, sizeof(double));
        verifyNotNULL(sum_of_dots_per_cluster)
    }

    /* initializing num_of_dots_per_cluster */
    for (i = 0; i < k; i++) {
        num_of_dots_per_cluster[i] = 0;
    }

    /* assignments in the sum and num arrays */
    for (i=0; i<n; i++){
        num_of_dots_per_cluster[dots_location[i]] += 1;
        for(j=0; j<d; j++)
            sum_of_dots_per_cluster[dots_location[i]][j] += dots[i][j];
    }

    /* assignments in the centroids array */
    for (i=0; i<k; i++){
        for(j=0; j<d; j++)
            centroids[i][j] = (sum_of_dots_per_cluster[i][j])/num_of_dots_per_cluster[i];
    }
    /* free space */
    free(num_of_dots_per_cluster);
    for(i = 0; i<k; i++)
        free(sum_of_dots_per_cluster[i]);
    free(sum_of_dots_per_cluster);
}

/*---------------------------------- kmeans Algorithm: ---------------------------------------------*/

double **kmeans(double **dots,double **centroids, int k, int d, int n, int max_iter) {

    /* variables declaration: */
    int i, j, iter_counter = 0;
    int *dots_location;
    double **old_centroids;

    /* Initialize a centroids matrix, and a dots cluster location array :*/
    dots_location = (int*) calloc(n, sizeof(int));
    old_centroids = (double**) malloc(k*sizeof(double*));
    verifyNotNULL(dots_location)
    verifyNotNULL(old_centroids)

    for(i=0; i<k; i++) {
        old_centroids[i] = (double *) malloc(d * sizeof(double));
        verifyNotNULL(old_centroids[i])
    }

    for(i=0; i<n; i++)
        dots_location[i] = -1;

    /* the kmeans algorithm flow :*/
    while(iter_counter < max_iter){

        for(i=0 ; i<n ; i++) {
            dots_location[i] = find_nearest_centroid(dots[i], centroids, k, d);
        }

        for(i=0 ; i<k ; i++)
            for(j=0 ; j<d ; j++)
                old_centroids[i][j] = centroids[i][j];

        update_centroids(dots, dots_location, centroids, n, d, k);

        if(check_equals_2d_list(old_centroids, centroids, k, d))
            break;

        iter_counter+=1;
    } /* end of while */
    /* frees */
    free(dots_location);
    for(i=0 ; i<k ; i++) {
        free(old_centroids[i]);
    }
    free(old_centroids);

    return centroids;
}