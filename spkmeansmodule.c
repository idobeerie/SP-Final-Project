#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "spkmeans.h"

#define MAX_KMEANS_ITERS 300


/* C Api Declarations: */
static double** create2DArrayFromPyObject(PyObject *data, int n, int d);
static PyObject* create2DPyObject(double **data, int n, int d);
#define verifyNotNULL(var) if((var)==NULL) {printf("An Error Has Occured"); exit(-1);}

static double** create2DArrayFromPyObject(PyObject *data, int n, int d) {
    int i, j;
    double** points;
    PyObject *temp_point,*inner_item;

    points = (double **) malloc(n * sizeof(double *));
    verifyNotNULL(points)

    for (i = 0; i < n; i++) {
        double *vector = malloc(d * sizeof(double));
        verifyNotNULL(vector)


        temp_point = PyList_GetItem(data, i);
        for (j = 0; j < d; j++) {
            inner_item = PyList_GetItem(temp_point, j);
            vector[j] = PyFloat_AsDouble(inner_item);
        }
        points[i] = vector;
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
    printf("Ido");
    PyObject *points;
    int n, d ;
    double **data_points;
    if (!PyArg_ParseTuple(args, "iiO", &n, &d, &points)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }
    data_points = create2DArrayFromPyObject(points, n, d);
    double **wam_matrix = wam(data_points, n, d);
    PyObject *result = create2DPyObject(wam_matrix, n, n);
    free_matrix(data_points, n);
    free_matrix(wam_matrix, n);
    return result;
}
static PyObject *ddg_api(PyObject *self, PyObject *args) {
    PyObject *points;
    int n, d ;
    double **data_points;
    if (!PyArg_ParseTuple(args, "iiO", &n, &d, &points)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }
    data_points = create2DArrayFromPyObject(points, n, d);
    double **wam_matrix = wam(data_points, n, d);
    double **ddg_matrix = ddg(wam_matrix, n);
    PyObject *result = create2DPyObject(ddg_matrix, n, n);
    free_matrix(data_points, n);
    free_matrix(wam_matrix, n);
    free_matrix(ddg_matrix, n);
    return result;
}
static PyObject *gl_api(PyObject *self, PyObject *args) {
    PyObject *points;
    int n, d ;
    double **data_points;
    if (!PyArg_ParseTuple(args, "iiO", &n, &d, &points)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }
    data_points = create2DArrayFromPyObject(points, n, d);
    double **wam_matrix = wam(data_points, n, d);
    double **ddg_matrix = ddg(wam_matrix, n);
    double **gl_matrix = gl(ddg_matrix, wam_matrix, n);
    PyObject *result = create2DPyObject(gl_matrix, n, n);
    free_matrix(data_points, n);
    free_matrix(wam_matrix, n);
    free_matrix(ddg_matrix, n);
    free_matrix(gl_matrix, n);
    return result;
}

static PyObject *jacobi_api(PyObject *self, PyObject *args) {
    PyObject *points;
    int n, d ;
    double **data_points;
    if (!PyArg_ParseTuple(args, "iiO", &n, &d, &points)) {
        printf("An Error Has Occured");
        Py_RETURN_NONE;
    }
    data_points = create2DArrayFromPyObject(points, n, d);
    Jacobi_output *jacobi_matrix = jacobi(data_points, n);
    double **jacobi_res = create2DfromJacobi(jacobi_matrix, n);
    PyObject *result = create2DPyObject(jacobi_res, n, n);
    free_matrix(data_points, n);
    free_matrix(jacobi_res, n);
    free(jacobi_matrix->eigenValues);
    free_matrix(jacobi_matrix->V, n);
    free(jacobi_matrix);
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
    initialCentroids = create2DArrayFromPyObject(pyInitialCentroids, k, d);
    points = create2DArrayFromPyObject(pyPoints, n, d);

    double **finalCentroids = kmeans(points, initialCentroids, k, d, n, MAX_KMEANS_ITERS);

    print_2d_array(finalCentroids, k, d);
    PyObject *result = create2DPyObject(finalCentroids, k, d);
    free_matrix(points, n);
    free_matrix(initialCentroids, k);
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