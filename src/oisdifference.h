#ifndef oisdifference_h
#define oisdifference_h

#include <stdlib.h>
#include <math.h>
#include <string.h>

void make_matrix_system(int nrows, int ncols, double* ref, double* sci, int w, int fwhm, int d,\
                        int nstars, int* xc, int* yc, double* C, double* D);
void solve_system(int n, double* C, double* D, double* xcs);
void var_convolve(int w, int d, int Q, double* a, int naxes, double* Ref, double* Con);
int perform_subtraction(int nrows, int ncols, double* ref, double* sci, int w, int fwhm, int d,\
                        int nstars, int* xc, int* yc, double* subtraction);

#endif /* oisdifference_h */
