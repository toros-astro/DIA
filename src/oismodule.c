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

    if (!PyArg_ParseTuple(args, "O!O!0!0!iii",
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

    int n = np_sciimage->dimensions[0];
    int m = np_sciimage->dimensions[1];
    int nstars = np_xc->dimensions[0];

    double* sciimage_data = (double*)np_sciimage->data;
    double* refimage_data = (double*)np_refimage->data;
    int* xc = (int *)np_xc->data;
    int* yc = (int *)np_yc->data;

    image img = {sciimage_data, n, m};
    image ref = {refimage_data, n, m};

    int poly_dof = (poly_deg + 1) * (poly_deg + 2) / 2; // number of degrees of freedom (ind var's) for the polynomial modulation
    int total_dof = kernel_side * kernel_side * poly_dof; // size of convolution matrix
    double *M = (double*) calloc(sizeof(double), total_dof * total_dof);
    double *b = (double*) calloc(sizeof(double), total_dof);

    make_matrix_system(ref, img, kernel_side / 2, stamp_side / 2, poly_deg, nstars, xc, yc, M, b);
    free(sciimage_data);

    double *xsol = (double*) malloc(sizeof(double) * total_dof);
    // This will solve the system Mx = b and store x in xsol
    solve_system(total_dof, M, b, xsol);
    free(M);
    free(b);
    double *ruinedimg_data = (double*) calloc(sizeof(double), n * m);
    var_convolve(kernel_side / 2, poly_deg, total_dof, xsol, n, refimage_data, ruinedimg_data);
    free(refimage_data);

    npy_intp con_dims[2] = {n, m};
    PyObject* py_conv_img = PyArray_SimpleNewFromData(2, con_dims, NPY_DOUBLE, ruinedimg_data);

    return Py_BuildValue("O", py_conv_img);
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
