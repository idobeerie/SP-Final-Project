#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

int free_memory(double* a, double* b, double* c, double** new_cluster, double** cluster_mean, double** data) {   //this is shit we dont want this 
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
        free(new_cluster);
    }
    if(cluster_mean) {
        free(cluster_mean);
    }
    if(data) {
        free(data);
    }
    return 0;
}

double euclidean_distance(double* a, double* b, size_t no_dims) {
    double dist = 0;
    for(size_t i = 0; i < no_dims; i++) {
        dist += pow(a[i] - b[i], 2);
        i++;
    }
    return sqrt(dist);
}

static PyObject* fit(PyObject* self, PyObject* args) {    //not fit maybe just call it spkcluster or something and see if this is ok with the code we did
    PyObject* centroids;
    PyObject* data;
    PyObject* point;
    int k;
    size_t iter = 300;
    double epsilon = 0.0;
    if(!PyArg_ParseTuple(args, "OOi", &centroids, &data)) {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return NULL;
    }
    size_t no_points = PyList_Size(data);
    point = PyList_GetItem(data, 0);
    size_t no_dims = PyList_Size(point);
    size_t no_centroids = PyList_Size(centroids);

    double* a = NULL;
    double* b = NULL;
    double* c = NULL;
    double** data = NULL;
    double** cluster_mean = NULL;
    double** new_cluster = NULL;
    size_t j, m, i;
    data = calloc(no_points, sizeof(double*));
    if(!data) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
        free_memory(b,c, a, new_cluster, cluster_mean, data);
        return NULL;
    }
    p = NULL;
    p = (double*)calloc(no_points * no_dims, sizeof(double));
    if(!p) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
        free_memory(b,c, a, new_cluster, cluster_mean, data);
        return NULL;
    }
    for(i = 0; i < no_points; i++) {
        data[i] = p + i * no_dims;
    }

    for(i = 0; i < no_points; i++) {
        point = PyList_GetItem(data, i);
        for(j = 0; j < no_dims; j++) {
            data[i][j] = PyFloat_AsDouble(PyList_GetItem(point, j));
        }
    }
    b = calloc(no_centroids * no_dims, sizeof(double));
    cluster_mean = calloc(no_centroids, sizeof(double*));
    if(!b || !cluster_mean) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
        free_memory(b,c, a, new_cluster, cluster_mean, data);
        return NULL;
    }
    for(i = 0; i < no_centroids; i++) {
        cluster_mean[i] = b + i * no_dims;
    }
    for(i = 0; i < no_centroids; i++) {
        point = PyList_GetItem(centroids, i);
        for(j = 0; j < no_dims; j++) {
            cluster_mean[i][j] = PyFloat_AsDouble(PyList_GetItem(point, j));
        }
    }

    double** curr_X = data;
    c = calloc(no_centroids * (1 + no_dims), sizeof(double));
    new_cluster = calloc(no_centroids, sizeof(double*));
    if(!c || !new_cluster) {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory");
        free_memory(b,c, a, new_cluster, cluster_mean, data);
        return NULL;
    }
    for(i = 0; i < no_centroids; i++) {
        new_cluster[i] = c + i * (1 + no_dims);
    }
    i = 0;
    int cvg = 0;
    while(i < iter && !cvg){
        memset(c, 0, no_centroids * (1 + no_dims) * sizeof(double));
        for(j = 0; j < no_points; j++) {
            double min_dist = DBL_MAX;
            size_t min_index = 0;
            for(m = 0; m < no_centroids; m++) {
                double dist = 0;
                for(size_t n = 0; n < no_dims; n++) {
                    dist += (curr_X[j][n] - cluster_mean[m][n]) * (curr_X[j][n] - cluster_mean[m][n]);
                }
                if(dist < min_dist) {
                    min_dist = dist;
                    min_index = m;
                }
            }
            double* min_cluster = new_cluster[min_index];
            for(m = 0; m < no_dims; m++) {
                new_cluster[min_index][m + 1] += curr_X[j][m];
            }
            min_cluster[no_dims] += 1;
        }
        for(j=0; j < no_centroids; j++) {
            for(m = 0; m < no_dims; m++) {
                new_cluster[j][m] /= new_cluster[j][no_dims];
            }
        }

        double max_dist = 0;
        double current_dist = 0;
        for(j = 0; j < no_centroids; j++) {
            current_dist = euclidean_distance(new_cluster[j], cluster_mean[j], no_dims);
            if(current_dist > max_dist) {
                max_dist = current_dist;
            }
        }
        if(max_dist < epsilon) {
            cvg = 1;
        }
        i++;
        for(j = 0; j < no_centroids; j++) {
            for(m = 0; m < no_dims; m++) {
                cluster_mean[j][m] = new_cluster[j][m];
            }
        }
    }
    PyObject* result = PyList_New(no_centroids);
    for(i = 0; i < no_centroids; i++) {
        PyObject* centroid = PyList_New(no_dims);
        PyList_SetItem(result, i, centroid);
        for(j = 0; j < no_dims; j++) {
            PyList_SetItem(centroid, j, PyFloat_FromDouble(cluster_mean[i][j]));
        }
    }
    free_memory(b,c, a, new_cluster, cluster_mean, data);  //change this we dont have this function !!!!
    return result;
}

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

    double** points_lst = allocateNonSquareMatrix(no_points, no_dims);
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

    double** points_lst = allocateNonSquareMatrix(no_points, no_dims);
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

    double** points_lst = allocateNonSquareMatrix(no_points, no_dims);
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
    double** gl_matrix = gl(ddg_matrix, wam_matrix, no_dims);
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
    free_matrix(points_lst, no_points)
    free_matrix(wam_matrix, no_points);
    free_matrix(ddg_matrix, no_points);
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

    double** points_lst = allocateNonSquareMatrix(no_points, no_dims);  // this will be square dont know if this is a problem
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
    
    double** jacobi_matrix = jacobi(points_lst, no_dims); //maybe we do need to call wan ddg and lg not sure
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
    {"spk",(PyCFunction)fit, METH_VARARGS, PyDoc_STR("takes 2 python lists, max iteration value, Convergence value")},
    {"wam",(PyCFunction)wam,METH_VARARGS,PyDoc_STR("takes points list(2D), returns the weight adjacency matrix")},
    {"ddg",(PyCFunction)ddg,METH_VARARGS,PyDoc_STR("takes points list(2D), returns the diagonal degree matrix")},
    {"gl",(PyCFunction)gl,METH_VARARGS,PyDoc_STR("takes points list(2D), returns the graph laplacian matrix")},
    {"jacobi",(PyCFunction)jacobi,METH_VARARGS,PyDoc_STR("jacobis on a symmetric matrix,second argument should be \"sorted\" for spk purposes\n, returns(values,vectors matrix,k)")},
      /*  The docstring for the function */
    {NULL, NULL, 0, NULL}     
    //dont we need here also the kmeans? 
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
