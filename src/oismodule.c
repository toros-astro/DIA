#include <Python.h>
#include <numpy/arrayobject.h>
#include "oisdifference.h"

#if PY_MAJOR_VERSION >= 3
#define PY3
#endif

static PyObject *
oismodule_subtract(PyObject *self, PyObject *args)
{
    PyArrayObject *np_sciimage, *np_refimage, *np_xc, *np_yc;
    int stamp_side;  // The stamp side 
    int kernel_side; // The kernel side in pixels
    int poly_deg; // The degree of the varying polynomial for the kernel

    if (!PyArg_ParseTuple(args, "O!O!O!O!iii",
        &PyArray_Type, &np_sciimage,
        &PyArray_Type, &np_refimage,
        &PyArray_Type, &np_xc,
        &PyArray_Type, &np_yc,
        &kernel_side, &poly_deg, &stamp_side)) {
        return NULL;
    }
    if (NULL == np_sciimage) return NULL;
    if (NULL == np_refimage) return NULL;
    if (NULL == np_xc) return NULL;
    if (NULL == np_yc) return NULL;

    double* sciimage = (double*)np_sciimage->data;
    double* refimage = (double*)np_refimage->data;
    int n = np_sciimage->dimensions[0];
    int m = np_sciimage->dimensions[1];

    int nstars = np_xc->dimensions[0];
    int* xc = (int *)np_xc->data;
    int* yc = (int *)np_yc->data;

    image img = {sciimage, n, m};
    image ref = {refimage, n, m};

    double *subt = (double *)malloc(n * m * sizeof(*subt));
    perform_subtraction(ref, img, kernel_side / 2, stamp_side / 2, poly_deg, nstars, xc, yc, subt);

    npy_intp diff_dims[2] = {n, m};
    PyObject* py_conv_img = PyArray_SimpleNewFromData(2, diff_dims, NPY_DOUBLE, subt);

    return Py_BuildValue("O", py_conv_img);
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
