#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <omp.h>
#include <string.h>
#include <time.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#include <lal/LALDetectors.h>
#include <lal/LALSimulation.h>
#include <lal/TimeDelay.h>
#include <lal/LALDatatypes.h>

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>


int compute_skymap(int input) {
    int output = input;
    return output;
}

static PyObject *method_compute_skymap(PyObject *self, PyObject *args) {
    int input;
    int output = -1;

    // parse python input arguments
    if (!PyArg_ParseTuple(args, "i", &input)) return NULL;

    if (input > 10) {
        PyErr_SetString(PyExc_ValueError, "input integer value must be less than 11");
        return NULL;
    }

    output = compute_skymap(input);
    return PyLong_FromLong(output);
}

// void capsule_cleanup(PyObject *capsule) {
//     /// Perform memory cleanup using free, see: https://stackoverflow.com/a/52732077
//     void *memory = PyCapsule_GetPointer(capsule, NULL);
//     free(memory);
// }

static PyObject* method_vector_add(PyObject* self, PyObject* args) {
    PyArrayObject* array1, * array2;

    if (!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &array1, &PyArray_Type, &array2))
        return NULL;

    // if (PyArray_NDIM(array1) != 1 || PyArray_NDIM(array2) != 1 || PyArray_DESCR(array1)->type_num != PyArray_DOUBLE || PyArray_DESCR(array2)->type_num != PyArray_DOUBLE) {
    //     PyErr_SetString(PyExc_ValueError, "arrays must be one-dimensional and of type float");
    //     return NULL;
    // }

    // int n1 = array1->dimensions[0];
    // int n2 = array2->dimensions[0];

    int n1 = PyArray_NDIM(array1);
    int n2 = PyArray_NDIM(array2);
    int size = (int) PyArray_SIZE(array1);
    int type = PyArray_TYPE(array1);

    printf("vector_add | size: %d, dim1: %d, dim2: %d\n",size, n1, n2);

    if (n1 != n2) {
        PyErr_SetString(PyExc_ValueError, "arrays must have the same length");
        return NULL;
    }

    // does this copy data? or just pass pointers
    double * array1_data = (double *) PyArray_DATA(array1);
    double * array2_data = (double *) PyArray_DATA(array2);
    double * output = (double *) malloc(sizeof(double) * n1);

    for (int i = 0; i < size; i++)
        output[i] = array1_data[i] + array2_data[i];
    
    PyObject *output_array = PyArray_SimpleNewFromData(n1, dims, type, (void*)output);
    PyArray_ENABLEFLAGS((PyArrayObject*) output_array, NPY_ARRAY_OWNDATA);

    return PyArray_Return(output_array); 
}

static PyMethodDef SkymapMethods[] = {
    {
        "compute_skymap",
        method_compute_skymap,
        METH_VARARGS,
        "Python interface for SPIIR Sky Map computation method."
    },
    {
        "vector_add",
        method_vector_add,
        METH_VARARGS,
        "Adds two numpy float arrays together on the CPU."
    },
    {NULL, NULL, 0, NULL}
};


const char SkymapDocString[] = "Python interface for SPIIR C sky map module.\n";


static struct PyModuleDef skymapmodule = {
    PyModuleDef_HEAD_INIT,
    "_skymap",
    SkymapDocString,
    -1,
    SkymapMethods,
};


PyMODINIT_FUNC PyInit__skymap(void) {
    PyObject *module = NULL;
    module = PyModule_Create(&skymapmodule);
    assert(! PyErr_Occurred() );
	import_array();

	if (PyErr_Occurred()) {
        return NULL;
    }

	return module;
}