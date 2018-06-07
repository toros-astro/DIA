#include <Python.h>
#include <numpy/ndarrayobject.h>
//#include <numpy/arrayobject.h>
#include "oisdifference.h"

#if PY_MAJOR_VERSION >= 3
#define PY3
#endif

static PyObject *
oismodule_subtract(PyObject *self, PyObject *args)
{

    PyObject *py_sciobj, *py_refobj, *py_xcobj, *py_ycobj;
    int stamp_side;  // The stamp side 
    int kernel_side; // The kernel side in pixels
    int poly_deg; // The degree of the varying polynomial for the kernel

    if (!PyArg_ParseTuple(args, "OOOOiii",
        &py_sciobj, &py_refobj, &py_xcobj, &py_ycobj,
        &kernel_side, &poly_deg, &stamp_side)) {
        return NULL;
    }

    PyObject *np_sciimage = PyArray_FROM_OTF(py_sciobj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (np_sciimage == NULL) {
        Py_XDECREF(np_sciimage);
        return NULL;
    }
    PyObject *np_refimage = PyArray_FROM_OTF(py_refobj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    if (np_refimage == NULL) {
        Py_XDECREF(np_sciimage);
        Py_XDECREF(np_refimage);
        return NULL;
    }
    PyObject *np_xc = PyArray_FROM_OTF(py_xcobj, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (np_sciimage == NULL) {
        Py_XDECREF(np_sciimage);
        Py_XDECREF(np_refimage);
        Py_XDECREF(np_xc);
        return NULL;
    }
    PyObject *np_yc = PyArray_FROM_OTF(py_ycobj, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (np_refimage == NULL) {
        Py_XDECREF(np_sciimage);
        Py_XDECREF(np_refimage);
        Py_XDECREF(np_xc);
        Py_XDECREF(np_yc);
        return NULL;
    }

    double* sciimage = (double*)PyArray_DATA(np_sciimage);
    double* refimage = (double*)PyArray_DATA(np_refimage);
    int n = (int)PyArray_DIM(np_sciimage, 0);
    int m = (int)PyArray_DIM(np_sciimage, 1);

    int nstars = (int)PyArray_DIM(np_xc, 0);
    int* xc = (int *)PyArray_DATA(np_xc);
    int* yc = (int *)PyArray_DATA(np_yc);

    double *subt = (double *)malloc(n * m * sizeof(*subt));
    perform_subtraction(n, m, refimage, sciimage, kernel_side / 2, stamp_side / 2, poly_deg, nstars, xc, yc, subt);

    Py_DECREF(np_sciimage);
    Py_DECREF(np_refimage);
    Py_DECREF(np_xc);
    Py_DECREF(np_yc);

    npy_intp diff_dims[2] = {n, m};
    PyArrayObject* py_subt_img = (PyArrayObject *)PyArray_SimpleNewFromData(2, diff_dims, NPY_DOUBLE, subt);
    PyArray_ENABLEFLAGS(py_subt_img, NPY_ARRAY_OWNDATA);

    return Py_BuildValue("N", py_subt_img);
}

static PyMethodDef OISModuleMethods[] = {
    {"subtract", oismodule_subtract, METH_VARARGS, "Perform Optimal Image Subtraction"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


#ifdef PY3
static struct PyModuleDef oismodule = {
   PyModuleDef_HEAD_INIT,
   "oismodule",   /* name of module */
   NULL, /* module documentation, may be NULL */
   -1,       /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
   OISModuleMethods
};

PyMODINIT_FUNC
PyInit_oismodule(void)
{
    PyObject *m;
    m = PyModule_Create(&oismodule);
    import_array();
    return m;
}
#else
PyMODINIT_FUNC
initoismodule(void)
{
    (void) Py_InitModule("oismodule", OISModuleMethods);
    import_array();
}
#endif
