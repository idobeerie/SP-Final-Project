#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"



typedef struct cord
{
    double value;
    struct cord *next;
}cord;

typedef struct vector
{
    struct vector *next;
    struct cord *cords;
    struct cord *prev_cords;
    int member_count;
    int choice;
}vector;



double distance(cord *cord1, cord *cord2){
    double d =0;
    cord *cptr1 = cord1;
    cord *cptr2 = cord2;
    while(cptr1 != NULL){
        d += pow((*cptr1).value - (*cptr2).value, 2);
        cptr1 = (*cptr1).next;
        cptr2 = (*cptr2).next;
    } 
    d = sqrt(d);
    return d;
}

void updateCord(cord *current, cord *newMember, int membersNum){
    double sum=0;
    cord *cptr1 = current;
    cord *cptr2 = newMember;
    while(cptr1 != NULL){
        sum = (*cptr1).value * membersNum + (*cptr2).value;
        (*cptr1).value = sum / (membersNum+1);
        cptr1 = (*cptr1).next;
        cptr2 = (*cptr2).next;
    }
}

void freeList(cord* head)
{
   cord* tmp;

   while (head != NULL)
    {
       tmp = head;
       head = head->next;
        free(tmp);
       
    }

}

void freeVector(vector* head)
{
   vector* temp;

   while (head != NULL)
    {
        freeList(head->cords);
        freeList(head->prev_cords);
        temp = head;
        head = head->next;
        free(temp);
    }

}


PyObject* fit(PyObject *centroids, PyObject *data, double epsi)
{ 
    PyObject *item, *centroid, *result_centroids, *result_centroid, *line, *sub_item;
    vector *head_vec, *curr_vec, *head_d_vec, *curr_d_vec;
    cord *head_cord, *curr_cord, *head_d_cord, *curr_d_cord, *prev_cord;
    int iter = 200, count =0, small_change=0, iter_count=0;
    double min_d = DBL_MAX, d=0.0;
    long i, j, k = 0, line_length=0, line_count=0, min_i=0, m=0;
    

    head_cord = malloc(sizeof(cord));
    curr_cord = head_cord;
    curr_cord->next = NULL;
    curr_cord->value =0;

    head_vec = malloc(sizeof(vector));
    curr_vec = head_vec;
    curr_vec->next = NULL;
    curr_vec->cords = NULL;
    curr_vec->prev_cords = NULL;
    curr_vec->member_count = 0;
    curr_vec->choice = 0;

    head_d_cord = malloc(sizeof(cord));
    curr_d_cord = head_d_cord;
    curr_d_cord->next = NULL;
    curr_d_cord->value = 0;

    head_d_vec = malloc(sizeof(vector));
    curr_d_vec = head_d_vec;
    curr_d_vec->next = NULL;
    curr_d_vec->cords = NULL;
    curr_d_vec->prev_cords = NULL;
    curr_d_vec->member_count = 0;
    curr_d_vec->choice = 0;

    k = PyList_Size(centroids);
    line_length = PyList_Size(PyList_GetItem(centroids, 0));

    for(i =0; i < k; i++){
        centroid = PyList_GetItem(centroids, i);
        for(j =0; j < line_length; j++){
            item = PyList_GetItem(centroid, j);
            curr_cord->value = PyFloat_AsDouble(item);
            curr_cord->next = malloc(sizeof(cord));
            curr_cord = curr_cord->next;
            curr_cord->next = NULL;
            curr_cord->value = 0;
        }
        curr_vec->member_count = 1;
        curr_vec->cords = head_cord;
        curr_vec->next = malloc(sizeof(vector));
        curr_vec->prev_cords = malloc(sizeof(cord));
        curr_vec->prev_cords->next = NULL;
        curr_vec = curr_vec->next;
        curr_vec->next = NULL;
        curr_vec->cords = NULL;
        curr_vec->prev_cords = NULL;
        curr_vec->member_count = 0;
        curr_vec->choice = 0;
        head_cord = malloc(sizeof(cord));
        curr_cord = head_cord;
        curr_cord->next = NULL;
        curr_cord->value = 0;
        curr_vec->cords = curr_cord;
        count++;
    }

    line_count = PyList_Size(data);

    for(i =0; i < line_count; i++){
        line = PyList_GetItem(data, i);
        for(j =0; j < line_length; j++){
            sub_item = PyList_GetItem(line, j);
            curr_d_cord->value = PyFloat_AsDouble(sub_item);
            curr_d_cord->next = malloc(sizeof(struct cord));
            curr_d_cord = curr_d_cord->next;
            curr_d_cord->next = NULL;
            curr_d_cord->value = 0;
        }
        curr_d_vec->cords = head_d_cord;
        curr_d_vec->next = malloc(sizeof(vector));
        curr_d_vec = curr_d_vec->next;
        curr_d_vec->next = NULL;
        curr_d_vec->cords = NULL;
        curr_d_vec->prev_cords = NULL;
        curr_d_vec->member_count = 0;
        curr_d_vec->choice = 0;
        head_d_cord = malloc(sizeof(cord));
        curr_d_cord = head_d_cord;
        curr_d_cord->next = NULL;
        curr_d_cord->value = 0;
        curr_d_vec->cords = curr_d_cord;
    }
    
    count =0;
    while(iter_count < iter){
        curr_d_vec = head_d_vec;
        curr_d_cord = curr_d_vec->cords;
        curr_vec = head_vec;
        curr_cord = head_vec->cords;
        while(count < line_count){
            min_d = DBL_MAX-1;
            for(i=0; i < k; i++){
                d = distance(curr_cord, curr_d_cord);
                if(d < min_d){
                    min_d = d;
                    min_i = i;
                }
                curr_vec = curr_vec->next;
                curr_cord = curr_vec->cords;
            }
            curr_d_vec->choice = min_i;            
            curr_vec = head_vec;
            curr_cord = head_vec->cords;
            curr_vec = head_vec;
            curr_cord = curr_vec->cords;
            curr_d_vec = curr_d_vec->next;
            curr_d_cord = curr_d_vec->cords;
            count++;
        }
        count = 0;
        curr_vec = head_vec;
        curr_cord = curr_vec->cords;
        prev_cord = curr_vec->prev_cords;
        for(i =0; i < k; i++){
            for(m =0; m < line_length; m++){
                prev_cord->value = curr_cord->value;
                if(iter_count == 0){
                    prev_cord->next = malloc(sizeof(cord));
                }
                prev_cord = prev_cord->next;
                prev_cord->value = 0;
                if(iter_count == 0){
                    prev_cord->next = NULL;
                }
                curr_cord = curr_cord->next;
            }
            curr_vec = curr_vec->next;
            curr_cord = curr_vec->cords;
            prev_cord = curr_vec->prev_cords;
        }
        curr_d_vec = head_d_vec;
        curr_vec = head_vec;
        for(i=0; i < k; i++){
            curr_cord = curr_vec->cords;
            while(curr_cord != NULL){
                curr_cord->value =0;
                curr_cord = curr_cord->next;
            }
            curr_vec->member_count = 0;
            curr_vec = curr_vec->next;
        }
        curr_d_vec = head_d_vec;
        while(count < line_count){
            curr_vec = head_vec;
            for(i=0; i < curr_d_vec->choice; i++){
                curr_vec = curr_vec->next;            
            }
            updateCord(curr_vec->cords, curr_d_vec->cords, curr_vec->member_count);
            curr_vec->member_count+=1;
            curr_d_vec = curr_d_vec->next;
            count++;       
        }
        count = 0;
        curr_vec = head_vec;
        if(iter_count > 0){
            small_change = 1;
            for(i =0; i < k; i++){
                if(distance(curr_vec->cords, curr_vec->prev_cords) >= epsi){
                    small_change = 0;
                }
                curr_vec = curr_vec->next;
            }
        }
        if(small_change == 1){
            break;
        }
        iter_count++;
    }

    curr_vec = head_vec;
    curr_cord = (*curr_vec).cords;
    result_centroids = PyList_New(k);
    for(i=0; i < k; i++){
        curr_cord = curr_vec->cords;
        result_centroid = PyList_New(line_length);
        for(j=0; j < line_length; j++){
            item = PyFloat_FromDouble(curr_cord->value);
            PyList_SetItem(result_centroid, j, item);
            curr_cord = curr_cord->next;
        }
        PyList_SetItem(result_centroids, i, result_centroid);
        curr_vec = (*curr_vec).next;
    }

    freeVector(head_vec);
    freeVector(head_d_vec);
    
    return Py_BuildValue("O", result_centroids);   
}


static PyObject* kmeanspp(PyObject *self, PyObject *args)
{
    PyObject *listOfCentroids, *listOfPoints;
    Py_ssize_t k, rows_num;
    int i=0;
    double epsi;
    int iter = 300;
    if (!PyArg_ParseTuple(args, "OOd", &listOfCentroids,&listOfPoints,&epsi))  
        return NULL;

    if (!PyList_Check(listOfCentroids) || !PyList_Check(listOfPoints)) {
        PyErr_SetString(PyExc_RuntimeError, "Received non-list type object.");
        return NULL;
    }
    k = PyList_GET_SIZE(listOfCentroids);

    for (i = 0; i < k; i++) {
        PyObject *listInList1 = PyList_GET_ITEM(listOfCentroids, i);

        if (!PyList_Check(listInList1)) {
            PyErr_SetString(PyExc_RuntimeError, "Non-list type found in list of lists.");
            return NULL;
        }


    }

    rows_num = PyList_GET_SIZE(listOfPoints);
    for (i = 0; i < rows_num; i++) {
        PyObject *listInList2 = PyList_GET_ITEM(listOfPoints, i);

        if (!PyList_Check(listInList2)) {
            PyErr_SetString(PyExc_RuntimeError, "Non-list type found in list of lists.");
            return NULL;
        }

    }
    return Py_BuildValue("O", fit(listOfCentroids, listOfPoints, epsi));
}

int free_memory(double* a, double* b, double* c, double** new_cluster, double** cluster_mean, double** data) {   
    if(a) {
        free(a);
    }
    if(b) {
        free(b);
    }
    if(c) {
        free(c);
    }
    if(new_cluster) {
        free_matrix(new_cluster);
    }
    if(cluster_mean) {
        free_matrix(cluster_mean);
    }
    if(data) {
        free_matrix(data);
    }
    return 0;
}

// double euclidean_distance(double* a, double* b, size_t no_dims) {
//     double dist = 0;
//     for(size_t i = 0; i < no_dims; i++) {
//         dist += pow(a[i] - b[i], 2);
//         i++;
//     }
//     return sqrt(dist);
// }

// static PyObject* fit(PyObject* self, PyObject* args) {    //not fit maybe just call it spkcluster or something and see if this is ok with the code we did
//     PyObject* centroids;
//     PyObject* point;
//     int k, pts_num;
//     size_t iter = 300;
//     double epsilon = 0.0;
//     if(!PyArg_ParseTuple(args, "Oi", &centroids, &pts_num)) {
//         PyErr_SetString(PyExc_TypeError, "Invalid arguments");
//         return NULL;
//     }
//     point = PyList_GetItem(centroids, 0);
//     size_t no_dims = PyList_Size(point);
//     size_t no_centroids = pts_num;

//     double* a = NULL;
//     double* b = NULL;
//     double* c = NULL;
//     double** pts = NULL;
//     double** cluster_mean = NULL;
//     double** new_cluster = NULL;
//     size_t j, m, i;
//     pts = calloc(pts_num, sizeof(double*));
//     if(!pts) {
//         PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
//         free_memory(b,c, a, new_cluster, cluster_mean, pts);
//         return NULL;
//     }
//     double* p = NULL;
//     p = (double*)calloc(pts_num * no_dims, sizeof(double));
//     if(!p) {
//         PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
//         free_memory(b,c, a, new_cluster, cluster_mean, pts);
//         return NULL;
//     }
//     for(i = 0; i < pts_num; i++) {
//         pts[i] = p + i * no_dims;
//     }

//     for(i = 0; i < pts_num; i++) {
//         point = PyList_GetItem(data, i);
//         for(j = 0; j < no_dims; j++) {
//             pts[i][j] = PyFloat_AsDouble(PyList_GetItem(point, j));
//         }
//     }
//     b = calloc(no_centroids * no_dims, sizeof(double));
//     cluster_mean = calloc(no_centroids, sizeof(double*));
//     if(!b || !cluster_mean) {
//         PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
//         free_memory(b,c, a, new_cluster, cluster_mean, pts);
//         return NULL;
//     }
//     for(i = 0; i < no_centroids; i++) {
//         cluster_mean[i] = b + i * no_dims;
//     }
//     for(i = 0; i < no_centroids; i++) {
//         point = PyList_GetItem(centroids, i);
//         for(j = 0; j < no_dims; j++) {
//             cluster_mean[i][j] = PyFloat_AsDouble(PyList_GetItem(point, j));
//         }
//     }

//     double** curr_X = pts;
//     c = calloc(no_centroids * (1 + no_dims), sizeof(double));
//     new_cluster = calloc(no_centroids, sizeof(double*));
//     if(!c || !new_cluster) {
//         PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
//         free_memory(b,c, a, new_cluster, cluster_mean, pts);
//         return NULL;
//     }
//     for(i = 0; i < no_centroids; i++) {
//         new_cluster[i] = c + i * (1 + no_dims);
//     }
//     i = 0;
//     int cvg = 0;
//     while(i < iter && !cvg){
//         memset(c, 0, no_centroids * (1 + no_dims) * sizeof(double));
//         for(j = 0; j < pts_num; j++) {
//             double min_dist = DBL_MAX;
//             size_t min_index = 0;
//             for(m = 0; m < no_centroids; m++) {
//                 double dist = 0;
//                 for(size_t n = 0; n < no_dims; n++) {
//                     dist += (curr_X[j][n] - cluster_mean[m][n]) * (curr_X[j][n] - cluster_mean[m][n]);
//                 }
//                 if(dist < min_dist) {
//                     min_dist = dist;
//                     min_index = m;
//                 }
//             }
//             double* min_cluster = new_cluster[min_index];
//             for(m = 0; m < no_dims; m++) {
//                 new_cluster[min_index][m + 1] += curr_X[j][m];
//             }
//             min_cluster[no_dims] += 1;
//         }
//         for(j=0; j < no_centroids; j++) {
//             for(m = 0; m < no_dims; m++) {
//                 new_cluster[j][m] /= new_cluster[j][no_dims];
//             }
//         }

//         double max_dist = 0;
//         double current_dist = 0;
//         for(j = 0; j < no_centroids; j++) {
//             current_dist = euclidean_distance(new_cluster[j], cluster_mean[j], no_dims);
//             if(current_dist > max_dist) {
//                 max_dist = current_dist;
//             }
//         }
//         if(max_dist < epsilon) {
//             cvg = 1;
//         }
//         i++;
//         for(j = 0; j < no_centroids; j++) {
//             for(m = 0; m < no_dims; m++) {
//                 cluster_mean[j][m] = new_cluster[j][m];
//             }
//         }
//     }
//     PyObject* result = PyList_New(no_centroids);
//     for(i = 0; i < no_centroids; i++) {
//         PyObject* centroid = PyList_New(no_dims);
//         PyList_SetItem(result, i, centroid);
//         for(j = 0; j < no_dims; j++) {
//             PyList_SetItem(centroid, j, PyFloat_FromDouble(cluster_mean[i][j]));
//         }
//     }
//     free_memory(b,c, a, new_cluster, cluster_mean, data);  //change this we dont have this function !!!!
//     return result;
// }

static PyObject* wam(PyObject* self, PyObject* args) {
    PyObject* data;
    if(!PyArg_ParseTuple(args, "O", &data )) {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return NULL;
    }
    size_t no_points = PyList_Size(data);
    PyObject* point = PyList_GetItem(data, 0);
    size_t no_dims = PyObject_Length(point);
    size_t i,j;

    double** points_lst = allocateMatrix(no_points, no_dims);
    if(!points_lst) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
        return NULL;
    }
    for(i = 0; i < no_points; i++) {
        point = PyList_GetItem(data, i);
        for(j = 0; j < no_dims; j++) {
            points_lst[i][j] = PyFloat_AsDouble(PyList_GetItem(point, j));
        }
    }
    double** wam_matrix = wam(points_lst, no_points, no_dims);
    if(!wam_matrix) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
        free_matrix(points_lst, no_points);
        return NULL;
    }
    PyObject* result = PyList_New(no_points);
    for(i = 0; i < no_points; i++) {
        PyObject* row = PyList_New(no_points);
        PyList_SetItem(result, i, row);
        for(j = 0; j < no_points; j++) {
            PyList_SetItem(row, j, PyFloat_FromDouble(wam_matrix[i][j]));
        }
    }
    free_matrix(points_lst, no_points);  //this is n on d dont think it matters
    free_matrix(wam_matrix, no_points);
    return result;
}


static PyObject* ddg(PyObject* self, PyObject* args) {
    PyObject* data;
    if(!PyArg_ParseTuple(args, "O", &data )) {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return NULL;
    }
    size_t no_points = PyList_Size(data);
    PyObject* point = PyList_GetItem(data, 0);
    size_t no_dims = PyObject_Length(point);
    size_t i,j;

    double** points_lst = allocateMatrix(no_points, no_dims);
    if(!points_lst) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
        return NULL;
    }
    for(i = 0; i < no_points; i++) {
        point = PyList_GetItem(data, i);
        for(j = 0; j < no_dims; j++) {
            points_lst[i][j] = PyFloat_AsDouble(PyList_GetItem(point, j));
        }
    }
    double** wam_matrix = wam(points_lst, no_points, no_dims);
    if(!wam_matrix) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
        free_matrix(points_lst, no_points);
        return NULL;
    }
    double** ddg_matrix = ddg(wam_matrix, no_points);
    if(!ddg_matrix) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
        free_matrix(points_lst, no_points);
        return NULL;
    }
    PyObject* result = PyList_New(no_points);
    for(i = 0; i < no_points; i++) {
        PyObject* row = PyList_New(no_points);
        PyList_SetItem(result, i, row);
        for(j = 0; j < no_points; j++) {
            PyList_SetItem(row, j, PyFloat_FromDouble(ddg_matrix[i][j]));
        }
    }
    free_matrix(points_lst, no_points);
    free_matrix(ddg_matrix, no_points);
    free_matrix(wam_matrix, no_points);
    return result;
}

static PyObject* gl(PyObject* self, PyObject* args) {
    PyObject* data;
    if(!PyArg_ParseTuple(args, "O", &data )) {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return NULL;
    }
    size_t no_points = PyList_Size(data);
    PyObject* point = PyList_GetItem(data, 0);
    size_t no_dims = PyObject_Length(point);
    size_t i,j;

    double** points_lst = allocateMatrix(no_points, no_dims);
    if(!points_lst) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
        return NULL;
    }
    for(i = 0; i < no_points; i++) {
        point = PyList_GetItem(data, i);
        for(j = 0; j < no_dims; j++) {
            points_lst[i][j] = PyFloat_AsDouble(PyList_GetItem(point, j));
        }
    }
    double** wam_matrix = wam(points_lst, no_points, no_dims);
    if(!wam_matrix) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
        free_matrix(points_lst, no_points);
        free_matrix(wam_matrix, no_points);
        return NULL;
    }
    double** ddg_matrix = ddg(wam_matrix, no_points);
    if(!ddg_matrix) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
        free_matrix(points_lst, no_points);
        free_matrix(wam_matrix, no_points);
        return NULL;
    }
    double** gl_matrix = gl(ddg_matrix, wam_matrix, no_dims);  // we free the wam and ddg matrix in gl
    if(!gl_matrix) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
        free_matrix(points_lst, no_points);
        free_matrix(wam_matrix, no_points);
        free_matrix(ddg_matrix, no_points);
        return NULL;
    }
    PyObject* result = PyList_New(no_points);
    for(i = 0; i < no_points; i++) {
        PyObject* row = PyList_New(no_points);
        PyList_SetItem(result, i, row);
        for(j = 0; j < no_points; j++) {
            PyList_SetItem(row, j, PyFloat_FromDouble(gl_matrix[i][j]));
        }
    }
    free_matrix(points_lst, no_points);
    free_matrix(gl_matrix, no_points);
    return result;
}

static PyObject* jacobi(PyObject* self, PyObject* args){
    PyObject* data;
    if(!PyArg_ParseTuple(args, "O", &data )) {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return NULL;
    }
    size_t no_points = PyList_Size(data);
    PyObject* point = PyList_GetItem(data, 0);
    size_t no_dims = PyObject_Length(point);
    size_t i,j;

    double** points_lst = allocateMatrix(no_points, no_dims);  
    if(!points_lst) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
        return NULL;
    }
    for(i = 0; i < no_points; i++) {
        point = PyList_GetItem(data, i);
        for(j = 0; j < no_dims; j++) {
            points_lst[i][j] = PyFloat_AsDouble(PyList_GetItem(point, j));
        }
    }
    
    double** jacobi_matrix = jacobi(points_lst, no_dims); 
    if(!jacobi_matrix) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
        free_matrix(points_lst, no_points);
        return NULL;
    }
    PyObject* result = PyList_New(no_points);
    for(i = 0; i < no_points; i++) {
        PyObject *row = PyList_New(no_points);
        PyList_SetItem(result, i, row);
        for (j = 0; j < no_points; j++) {
            PyList_SetItem(row, j, PyFloat_FromDouble(jacobi_matrix[i][j]));
        }
    }
    free_matrix(points_lst, no_points);
    free_matrix(jacobi_matrix, no_points);   
    return result;
}


static PyMethodDef methods[] = {
    {"spk",(PyCFunction)kmeanspp, METH_VARARGS, PyDoc_STR("takes 2 python lists, max iteration value, Convergence value")},
    {"wam",(PyCFunction)wam,METH_VARARGS,PyDoc_STR("takes points list(2D), returns the weight adjacency matrix")},
    {"ddg",(PyCFunction)ddg,METH_VARARGS,PyDoc_STR("takes points list(2D), returns the diagonal degree matrix")},
    {"gl",(PyCFunction)gl,METH_VARARGS,PyDoc_STR("takes points list(2D), returns the graph laplacian matrix")},
    {"jacobi",(PyCFunction)jacobi,METH_VARARGS,PyDoc_STR("jacobis on a symmetric matrix,second argument should be \"sorted\" for spk purposes\n, returns(values,vectors matrix,k)")},
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