# define PY_SSIZE_T_CLEAN
# include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

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
};

static PyMethodDef kmeansMethods[] = {
    {"fit",                   /* the Python method name that will be used */
      (PyCFunction) kmeanspp, /* the C-function that implements the Python function and returns static PyObject*  */
      METH_VARARGS,           /* flags indicating parameters
accepted for this function */
      PyDoc_STR("fit expects the initial centroids and datapoints")}, /*  The docstring for the function */
    {NULL, NULL, 0, NULL}     /* The last entry must be all NULL as shown to act as a
                                 sentinel. Python looks for this entry to know that all
                                 of the functions for the module have been defined. */
};


static struct PyModuleDef kmeansmodule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp", /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    kmeansMethods /* the PyMethodDef array from before containing the methods of the extension */
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&kmeansmodule);
    if (!m) {
        return NULL;
    }
    return m;
}