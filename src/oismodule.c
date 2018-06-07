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


    int poly_dof = (poly_deg + 1) * (poly_deg + 2) / 2; // number of degree elements //
    int M_dim = pow(kernel_side, 2) * poly_dof; // size of convolution matrix//
    double *M_matrix = (double*) calloc(sizeof(double), M_dim * M_dim);
    double *b_vector = (double*) calloc(sizeof(double), M_dim);

    // Now we need to make stamps around each star to find the parameters for the kernel //
    make_matrix_system(n, m, refimage, sciimage, kernel_side / 2, stamp_side / 2,\
        poly_deg, nstars, xc, yc, M_matrix, b_vector);
    double *sol_vector = (double*) malloc(sizeof(double) * M_dim);
    // This will solve the system Mx = b and store x in sol
    solve_system(M_dim, M_matrix, b_vector, sol_vector);
    free(b_vector);
    free(M_matrix);
    int nelems = n * m;
    double *conv = (double*) calloc(sizeof(double), nelems);

    var_convolve(kernel_side / 2, poly_deg, M_dim, sol_vector, n, refimage, conv);

    double *subt = (double *)malloc(n * m * sizeof(*subt));
    // Perform the subtraction //
    for (int i = 0; i < nelems; i++) {
        subt[i] = sciimage[i] - conv[i];
    }

    Py_DECREF(np_sciimage);
    Py_DECREF(np_refimage);
    Py_DECREF(np_xc);
    Py_DECREF(np_yc);

    npy_intp diff_dims[2] = {n, m};
    PyArrayObject* py_subt_img = (PyArrayObject *)PyArray_SimpleNewFromData(2, diff_dims, NPY_DOUBLE, subt);
    PyArray_ENABLEFLAGS(py_subt_img, NPY_ARRAY_OWNDATA);

    npy_intp opt_dims[2] = {n, m};
    PyArrayObject* py_opt_img = (PyArrayObject *)PyArray_SimpleNewFromData(2, opt_dims, NPY_DOUBLE, conv);
    PyArray_ENABLEFLAGS(py_opt_img, NPY_ARRAY_OWNDATA);

    npy_intp krn_dims[3] = {kernel_side, kernel_side, poly_dof};
    PyArrayObject* py_krn = (PyArrayObject *)PyArray_SimpleNewFromData(2, krn_dims, NPY_DOUBLE, sol_vector);
    PyArray_ENABLEFLAGS(py_krn, NPY_ARRAY_OWNDATA);    

    return Py_BuildValue("NNN", py_subt_img, py_opt_img, py_krn);
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
