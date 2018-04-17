#include <Python.h>
#include <numpy/arrayobject.h>
#include "oisdifference.h"

#if PY_MAJOR_VERSION >= 3
#define PY3
#endif

static PyObject *
oismodule_subtract(PyObject *self, PyObject *args)
{
    PyArrayObject *np_image, *np_refimage;
    int deg; // The degree of the varying polynomial for the kernel
    int stamp_side; // The stamp side 

    if (!PyArg_ParseTuple(args, "O!O!ii", &PyArray_Type, &np_image,
            &PyArray_Type, &np_refimage,
            &deg, &stamp_side)) {
        return NULL;
    }
    if (NULL == np_image) return NULL;
    if (NULL == np_refimage) return NULL;

    int n = np_image->dimensions[0];
    int m = np_image->dimensions[1];

    double* image = (double*)np_image->data;
    double* refimage = (double*)np_refimage->data;

    double* addition = (double*)malloc(n * m * sizeof(*addition));
    for(int i = 0; i < n*m; i++) {
        addition[i] = image[i] + refimage[i];
    }

    npy_intp add_dims[2] = {n, m};
    PyObject* add = PyArray_SimpleNewFromData(2, add_dims, NPY_DOUBLE, addition);

    return Py_BuildValue("Os", add, "Hello from the other side");
}

static PyMethodDef OISModuleMethods[] = {
    {"subtract", oismodule_subtract, METH_VARARGS, "Add two arrays."},
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
