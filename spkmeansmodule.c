#include "spkmeansmodule.h"

#define MAX_KMEANS_ITERS 300


/* C Api Declarations: */
static double** create2DArrayFromPyObject(PyObject *data, int n, int d);
static PyObject* create2DPyObject(double **data, int n, int d);
#define verifyNotNULL(var) if((var)==NULL) {printf("An Error Has Occured"); exit(-1);}
/*
this is the part of out last HW2
*/

// typedef struct cord
// {
//     double value;
//     struct cord *next;
// }cord;

// typedef struct vector
// {
//     struct vector *next;
//     struct cord *cords;
//     struct cord *prev_cords;
//     int member_count;
//     int choice;
// }vector;



// double distance(cord *cord1, cord *cord2){
//     double d =0;
//     cord *cptr1 = cord1;
//     cord *cptr2 = cord2;
//     while(cptr1 != NULL){
//         d += pow((*cptr1).value - (*cptr2).value, 2);
//         cptr1 = (*cptr1).next;
//         cptr2 = (*cptr2).next;
//     } 
//     d = sqrt(d);
//     return d;
// }

// void updateCord(cord *current, cord *newMember, int membersNum){
//     double sum=0;
//     cord *cptr1 = current;
//     cord *cptr2 = newMember;
//     while(cptr1 != NULL){
//         sum = (*cptr1).value * membersNum + (*cptr2).value;
//         (*cptr1).value = sum / (membersNum+1);
//         cptr1 = (*cptr1).next;
//         cptr2 = (*cptr2).next;
//     }
// }

// void freeList(cord* head)
// {
//    cord* tmp;

//    while (head != NULL)
//     {
//        tmp = head;
//        head = head->next;
//         free(tmp);
       
//     }

// }

// void freeVector(vector* head)
// {
//    vector* temp;

//    while (head != NULL)
//     {
//         freeList(head->cords);
//         freeList(head->prev_cords);
//         temp = head;
//         head = head->next;
//         free(temp);
//     }

// }


// PyObject* fit(PyObject *centroids, PyObject *data, double epsi)
// { 
//     PyObject *item, *centroid, *result_centroids, *result_centroid, *line, *sub_item;
//     vector *head_vec, *curr_vec, *head_d_vec, *curr_d_vec;
//     cord *head_cord, *curr_cord, *head_d_cord, *curr_d_cord, *prev_cord;
//     int iter = 200, count =0, small_change=0, iter_count=0;
//     double min_d = DBL_MAX, d=0.0;
//     long i, j, k = 0, line_length=0, line_count=0, min_i=0, m=0;
    

//     head_cord = malloc(sizeof(cord));
//     curr_cord = head_cord;
//     curr_cord->next = NULL;
//     curr_cord->value =0;


//     head_vec = malloc(sizeof(vector));
//     curr_vec = head_vec;
//     curr_vec->next = NULL;
//     curr_vec->cords = NULL;
//     curr_vec->prev_cords = NULL;
//     curr_vec->member_count = 0;
//     curr_vec->choice = 0;

//     head_d_cord = malloc(sizeof(cord));
//     curr_d_cord = head_d_cord;
//     curr_d_cord->next = NULL;
//     curr_d_cord->value = 0;

//     head_d_vec = malloc(sizeof(vector));
//     curr_d_vec = head_d_vec;
//     curr_d_vec->next = NULL;
//     curr_d_vec->cords = NULL;
//     curr_d_vec->prev_cords = NULL;
//     curr_d_vec->member_count = 0;
//     curr_d_vec->choice = 0;

//     k = PyList_Size(centroids);
//     line_length = PyList_Size(PyList_GetItem(centroids, 0));

//     for(i =0; i < k; i++){
//         centroid = PyList_GetItem(centroids, i);
//         for(j =0; j < line_length; j++){
//             item = PyList_GetItem(centroid, j);
//             curr_cord->value = PyFloat_AsDouble(item);
//             curr_cord->next = malloc(sizeof(cord));
//             curr_cord = curr_cord->next;
//             curr_cord->next = NULL;
//             curr_cord->value = 0;
//         }
//         curr_vec->member_count = 1;
//         curr_vec->cords = head_cord;
//         curr_vec->next = malloc(sizeof(vector));
//         curr_vec->prev_cords = malloc(sizeof(cord));
//         curr_vec->prev_cords->next = NULL;
//         curr_vec = curr_vec->next;
//         curr_vec->next = NULL;
//         curr_vec->cords = NULL;
//         curr_vec->prev_cords = NULL;
//         curr_vec->member_count = 0;
//         curr_vec->choice = 0;
//         head_cord = malloc(sizeof(cord));
//         curr_cord = head_cord;
//         curr_cord->next = NULL;
//         curr_cord->value = 0;
//         curr_vec->cords = curr_cord;
//         count++;
//     }

//     line_count = PyList_Size(data);

//     for(i =0; i < line_count; i++){
//         line = PyList_GetItem(data, i);
//         for(j =0; j < line_length; j++){
//             sub_item = PyList_GetItem(line, j);
//             curr_d_cord->value = PyFloat_AsDouble(sub_item);
//             curr_d_cord->next = malloc(sizeof(struct cord));
//             curr_d_cord = curr_d_cord->next;
//             curr_d_cord->next = NULL;
//             curr_d_cord->value = 0;
//         }
//         curr_d_vec->cords = head_d_cord;
//         curr_d_vec->next = malloc(sizeof(vector));
//         curr_d_vec = curr_d_vec->next;
//         curr_d_vec->next = NULL;
//         curr_d_vec->cords = NULL;
//         curr_d_vec->prev_cords = NULL;
//         curr_d_vec->member_count = 0;
//         curr_d_vec->choice = 0;
//         head_d_cord = malloc(sizeof(cord));
//         curr_d_cord = head_d_cord;
//         curr_d_cord->next = NULL;
//         curr_d_cord->value = 0;
//         curr_d_vec->cords = curr_d_cord;
//     }
    

//     count =0;
//     while(iter_count < iter){
//         curr_d_vec = head_d_vec;
//         curr_d_cord = curr_d_vec->cords;
//         curr_vec = head_vec;
//         curr_cord = head_vec->cords;
//         while(count < line_count){
//             min_d = DBL_MAX-1;
//             for(i=0; i < k; i++){
//                 d = distance(curr_cord, curr_d_cord);
//                 if(d < min_d){
//                     min_d = d;
//                     min_i = i;
//                 }
//                 curr_vec = curr_vec->next;
//                 curr_cord = curr_vec->cords;
//             }
//             curr_d_vec->choice = min_i;            
//             curr_vec = head_vec;
//             curr_cord = head_vec->cords;
//             curr_vec = head_vec;
//             curr_cord = curr_vec->cords;
//             curr_d_vec = curr_d_vec->next;
//             curr_d_cord = curr_d_vec->cords;
//             count++;
//         }
//         count = 0;
//         curr_vec = head_vec;
//         curr_cord = curr_vec->cords;
//         prev_cord = curr_vec->prev_cords;
//         for(i =0; i < k; i++){
//             for(m =0; m < line_length; m++){
//                 prev_cord->value = curr_cord->value;
//                 if(iter_count == 0){
//                     prev_cord->next = malloc(sizeof(cord));
//                 }
//                 prev_cord = prev_cord->next;
//                 prev_cord->value = 0;
//                 if(iter_count == 0){
//                     prev_cord->next = NULL;
//                 }
//                 curr_cord = curr_cord->next;
//             }
//             curr_vec = curr_vec->next;
//             curr_cord = curr_vec->cords;
//             prev_cord = curr_vec->prev_cords;
//         }
//         curr_d_vec = head_d_vec;
//         curr_vec = head_vec;
//         for(i=0; i < k; i++){
//             curr_cord = curr_vec->cords;
//             while(curr_cord != NULL){
//                 curr_cord->value =0;
//                 curr_cord = curr_cord->next;
//             }
//             curr_vec->member_count = 0;
//             curr_vec = curr_vec->next;
//         }
//         curr_d_vec = head_d_vec;
//         while(count < line_count){
//             curr_vec = head_vec;
//             for(i=0; i < curr_d_vec->choice; i++){
//                 curr_vec = curr_vec->next;            
//             }
//             updateCord(curr_vec->cords, curr_d_vec->cords, curr_vec->member_count);
//             curr_vec->member_count+=1;
//             curr_d_vec = curr_d_vec->next;
//             count++;       
//         }
//         count = 0;
//         curr_vec = head_vec;
//         if(iter_count > 0){
//             small_change = 1;
//             for(i =0; i < k; i++){
//                 if(distance(curr_vec->cords, curr_vec->prev_cords) >= epsi){
//                     small_change = 0;
//                 }
//                 curr_vec = curr_vec->next;
//             }
//         }
//         if(small_change == 1){
//             break;
//         }
//         iter_count++;
//     }

//     curr_vec = head_vec;
//     curr_cord = (*curr_vec).cords;
//     result_centroids = PyList_New(k);
//     for(i=0; i < k; i++){
//         curr_cord = curr_vec->cords;
//         result_centroid = PyList_New(line_length);
//         for(j=0; j < line_length; j++){
//             item = PyFloat_FromDouble(curr_cord->value);
//             PyList_SetItem(result_centroid, j, item);
//             curr_cord = curr_cord->next;
//         }
//         PyList_SetItem(result_centroids, i, result_centroid);
//         curr_vec = (*curr_vec).next;
//     }

//     freeVector(head_vec);
//     freeVector(head_d_vec);
    
//     return Py_BuildValue("O", result_centroids);   
// }


// static PyObject* kmeanspp(PyObject *self, PyObject *args)
// {
//     PyObject *listOfCentroids, *listOfPoints;
//     Py_ssize_t k, rows_num;
//     int i=0;
//     double epsi = 0.0;

//     if (!PyArg_ParseTuple(args, "OO", &listOfCentroids,&listOfPoints))
//         return NULL;

//     if (!PyList_Check(listOfCentroids) || !PyList_Check(listOfPoints)) {
//         PyErr_SetString(PyExc_RuntimeError, "Received non-list type object.");
//         return NULL;
//     }
//     k = PyList_GET_SIZE(listOfCentroids);

//     for (i = 0; i < k; i++) {
//         PyObject *listInList1 = PyList_GET_ITEM(listOfCentroids, i);

//         if (!PyList_Check(listInList1)) {
//             PyErr_SetString(PyExc_RuntimeError, "Non-list type found in list of lists.");
//             return NULL;
//         }


//     }

//     rows_num = PyList_GET_SIZE(listOfPoints);
//     for (i = 0; i < rows_num; i++) {
//         PyObject *listInList2 = PyList_GET_ITEM(listOfPoints, i);

//         if (!PyList_Check(listInList2)) {
//             PyErr_SetString(PyExc_RuntimeError, "Non-list type found in list of lists.");
//             return NULL;
//         }

//     }
//     return Py_BuildValue("O", fit(listOfCentroids, listOfPoints, epsi));
// };



//end of test
static double** create2DArrayFromPyObject(PyObject *data, int n, int d) {
    int i, j;
    double** points;
    PyObject *temp_point,*inner_item;

    points = alloc_nXm_matrix(n,d);
    verifyNotNULL(points)
    for (i = 0; i < n; i++) {
        temp_point = PyList_GetItem(data, i);
        for (j = 0; j < d; j++) {
            inner_item = PyList_GetItem(temp_point, j);
            points[i][j] = PyFloat_AsDouble(inner_item);
        }
    }
    return points;
}
static PyObject *create2DPyObject(double** matrix, int n, int d) {
    int i,j;
    PyObject *currentVector, *pyMatrix,*num;
    pyMatrix = PyList_New(n);
    verifyNotNULL(pyMatrix)
    for (i = 0; i < n; i++) {
        currentVector = PyList_New(d);
        verifyNotNULL(currentVector)
        for (j = 0; j < d; j++) {
            num = PyFloat_FromDouble(matrix[i][j]);
            verifyNotNULL(num)
            PyList_SET_ITEM(currentVector, j, num);
        }
    PyList_SET_ITEM(pyMatrix, i, currentVector);
    }
    return pyMatrix;
}



/*--------------------------------- c-api -----------------------------------------------------*/
static PyObject *wam_api(PyObject *self, PyObject *args) {
    PyObject *points;
    int n, d ;
    double **data_points;
    int print;
    if (!PyArg_ParseTuple(args, "iiOi", &n, &d, &points, &print)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }
    data_points = create2DArrayFromPyObject(points, n, d);
    double **wam_matrix = wam(data_points, n, d);
    PyObject *result = create2DPyObject(wam_matrix, n, n);
    if (print == 1){
        print_matrix(wam_matrix, n, n);}
    free_matrix(data_points);
    free_matrix(wam_matrix);
    return result;
}
static PyObject *ddg_api(PyObject *self, PyObject *args) {
    PyObject *points;
    int n, d ;
    double **data_points;
    int print;
    if (!PyArg_ParseTuple(args, "iiOi", &n, &d, &points, &print)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }
    data_points = create2DArrayFromPyObject(points, n, d);
    double **wam_matrix = wam(data_points, n, d);
    double **ddg_matrix = ddg(wam_matrix, n);
    PyObject *result = create2DPyObject(ddg_matrix, n, n);
    if (print == 1){
        print_matrix(ddg_matrix, n, n);}
    free_matrix(data_points);
    free_matrix(wam_matrix);
    free_matrix(ddg_matrix);
    return result;
}
static PyObject *gl_api(PyObject *self, PyObject *args) {
    PyObject *points;
    int n, d, p ;
    double **data_points;
    if (!PyArg_ParseTuple(args, "iiOi", &n, &d, &points, &p)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }
    data_points = create2DArrayFromPyObject(points, n, d);
    double **wam_matrix = wam(data_points, n, d);
    double **ddg_matrix = ddg(wam_matrix, n);
    double **gl_matrix = gl_py(ddg_matrix, wam_matrix, n);
    PyObject *result = create2DPyObject(gl_matrix, n, n);
    if (p == 1){
        print_matrix(gl_matrix, n, n);}
    // free_matrix(data_points);
    // free_matrix(wam_matrix);
    // free_matrix(ddg_matrix);
    // free_matrix(gl_matrix);
    return result;
}

static PyObject *jacobi_api(PyObject *self, PyObject *args) {
    PyObject *points;
    int n, d ;
    double **data_points;
    int p;
    if (!PyArg_ParseTuple(args, "iiOi", &n, &d, &points, &p)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }
    data_points = create2DArrayFromPyObject(points, n, d);
    Jacobi_output *jacobi_matrix = jacobi(data_points, n);
    double **jacobi_res = create2DfromJacobi(jacobi_matrix, n);
    if (p == 1){  
        print_matrix(jacobi_res, n+1, n);}
    PyObject *result = create2DPyObject(jacobi_res, n+1, n);
    // free_matrix(data_points);
    // free_matrix(jacobi_res);
    // free(jacobi_matrix->eigenValues);
    // free_matrix(jacobi_matrix->V);
    // free(jacobi_matrix);
    return result;
}

static PyObject *kmeans_pp_api(PyObject *self, PyObject *args) {
    PyObject *pyInitialCentroids, *pyPoints;
    int k, n, d;
    double **initialCentroids, **points;
    if (!PyArg_ParseTuple(args, "iiiOO", &k, &n, &d, &pyInitialCentroids, &pyPoints)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }
    initialCentroids = create2DArrayFromPyObject(pyInitialCentroids, d, d);
    printf("after _init");
    points = create2DArrayFromPyObject(pyPoints, n, d);
    printf("after points");
    double **finalCentroids = kmeans(points, initialCentroids, k, d, n, MAX_KMEANS_ITERS);
    PyObject *result = create2DPyObject(finalCentroids, k, d);
    free_matrix(points);
    free_matrix(initialCentroids);
    return result;
}



static PyMethodDef methods[] = {
    {"spk",(PyCFunction)kmeans_pp_api, METH_VARARGS, PyDoc_STR("takes 2 python lists, max iteration value, Convergence value")},
    {"wam",(PyCFunction)wam_api,METH_VARARGS,PyDoc_STR("takes points list(2D), returns the weight adjacency matrix")},
    {"ddg",(PyCFunction)ddg_api,METH_VARARGS,PyDoc_STR("takes points list(2D), returns the diagonal degree matrix")},
    {"gl",(PyCFunction)gl_api,METH_VARARGS,PyDoc_STR("takes points list(2D), returns the graph laplacian matrix")},
    {"jacobi",(PyCFunction)jacobi_api,METH_VARARGS,PyDoc_STR("jacobis on a symmetric matrix,second argument should be \"sorted\" for spk purposes\n, returns(values,vectors matrix,k)")},
      /*  The docstring for the function */
    {NULL, NULL, 0, NULL}     

};

static struct PyModuleDef mykmeanssp = {
    PyModuleDef_HEAD_INIT,
    "spkmeansmodule", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    methods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_mykmeanssp(void) {
    PyObject* m;
    m = PyModule_Create(&mykmeanssp);
    if (!m) {
        return NULL;
    }
    return m;
}