#include "oisdifference.h"


int perform_subtraction(int nrows, int ncols, double* ref, double* sci, int w, int fwhm,\
                         int d, int nstars, int* xc, int* yc, double* subtraction) {
    int deg = (d + 1) * (d + 2) / 2; // number of degree elements //
    int Q = pow(2 * w + 1, 2) * deg; // size of convolution matrix//
    double *C = (double*) calloc(sizeof(double), Q * Q);
    double *D = (double*) calloc(sizeof(double), Q);

    // Now we need to make stamps around each star to find the parameters for the kernel //
    make_matrix_system(nrows, ncols, ref, sci, w, fwhm, d, nstars, xc, yc, C, D);
    double *a = (double*) malloc(sizeof(double) * Q);
    // This will solve the system Cx = D and store x in a
    solve_system(Q, C, D, a);
    free(D);
    free(C);
    int nelems = nrows * ncols;
    double *conv = (double*) calloc(sizeof(double), nelems);

    var_convolve(w, d, Q, a, nrows, ref, conv);

    // Perform the subtraction //
    for (int i = 0; i < nelems; i++) {
        subtraction[i] = sci[i] - conv[i];
    }
    free(conv);
    return EXIT_SUCCESS;
}

void make_matrix_system(int nrows, int ncols, double* ref, double* sci, \
    int w, int fwhm, int d, int nstars, int* xc, int* yc, double* C, double* D) {
    
    // Now we need to make stamps around each star to find the parameters for the kernel //
    //parameters that fall out from above//
    
    int L = 2 * w + 1;       // kernel axis //
    int nk = L * L;          // number of kernel elements //
    int stax = 2 * fwhm + 1; // size of star stamps //
    int S = stax * stax;     // number of stamp elements //
    int deg = (d + 1) * (d + 2) / 2; // number of degree elements //
    int Q = nk * deg;          //size of D matrix//
    int cent = (nk - 1) / 2; //center of the kernel//
    
    double *Rs = (double*) malloc(sizeof(double) * S);
    double *Ss = (double*) malloc(sizeof(double) * S);
    double *CRKq = (double*) malloc(sizeof(double) * S);
    double *CRKn = (double*) malloc(sizeof(double) * S);
    double *Kn = (double*) malloc(sizeof(double) * nk);
    double *Kq = (double*) malloc(sizeof(double) * nk);
    
    // now we need to solve for the kernel parameters //
    int qrs = 0; //initialize the qrs step//
    for (int q = 0; q < nk; q++) {
        //make the q kernel//
        memset(Kq, 0, nk * sizeof(*Kq));
        Kq[q] = 1.0;
        if (q != cent) Kq[cent] = -1.0;
        
        for (int r = 0; r <= d; r++){
            for (int s = 0; s <= d - r; s++){
                for (int n = 0; n < nk; n++){
                    //make the n kernel//
                    memset(Kn, 0, nk * sizeof(*Kn));
                    Kn[n] = 1.0;
                    if (n != cent) Kn[cent] = -1.0;
                    int ml = 0; //initialize the ml step//
                    for (int m = 0; m <= d; m++){
                        for (int l = 0; l <= d - m; l++){
                            D[qrs] = 0; //ensure D is only calculated once for each Q//
                            for (int k = 0; k < nstars; k++){
                                int xcent = xc[k]; //x coordinate of stamp center//
                                int ycent = yc[k]; //y coordinate of stamp center//
                                //make the star stamps//
                                for (int i = 0; i < stax; i++){
                                    for(int j = 0; j < stax; j++){
                                        Rs[i + j * stax] = ref[(i + xcent - fwhm) + (j + ycent - fwhm) * nrows];
                                        Ss[i + j * stax] = sci[(i + xcent - fwhm) + (j + ycent - fwhm) * nrows];
                                    }
                                }
                                //reinitialize the convolution matrix//
                                memset(CRKn, 0, S * sizeof(*CRKn));
                                memset(CRKq, 0, S * sizeof(*CRKq));
                                //now we do the convolution for n and q//
                                for (int i=0; i<stax;i++){
                                    for(int j=0;j<stax;j++){
                                        for (int mm = 0; mm < L; mm++){
                                            for(int nn = 0; nn < L; nn++){
                                                int ii = i + (mm - w);//index of convolution//
                                                int jj = j + (nn - w);//index of convolution//
                                                if (ii >= 0 && ii < stax && jj >= 0 && jj < stax) {
                                                    CRKn[i + j * stax] += Rs[ii + jj * stax] * Kn[mm + nn * L];
                                                    CRKq[i + j * stax] += Rs[ii + jj * stax] * Kq[mm + nn * L];
                                                }//end of if statement//
                                            }//end of nn loop//
                                        }// end of mm loop//
                                    }//end of j loop//
                                }//end of i loop//
                                
                                int mr = m + r;
                                int ls = l + s; //exponents for polynomial approximation//
                                //now we need to fill in C//
                                for (int i = 0; i < S; i++){
                                    C[n * deg + ml + qrs * Q] += pow(xcent, mr) * pow(ycent, ls) * CRKn[i] * CRKq[i];
                                }//end of C loop//
                                
                                //now we need to fill in D//
                                for (int i = 0; i < S; i++){
                                    D[qrs] += pow(xcent, r) * pow(ycent, s) * Ss[i] * CRKq[i];
                                }//end of D loop//
                                
                            }//end of k loop//
                            ml++;
                        }//end of l loop//
                    }//end of m loop//
                }//end of n loop//
                qrs++;
            }//end of s loop//
        }//end of r loop //
    }//end of q loop//
    
    //free everything//
    free(CRKn);
    free(CRKq);
    free(Kn);
    free(Kq);
    free(Ss);
    free(Rs);
}

void solve_system(int n, double* C, double* D, double* xcs) {
    int nsq = n * n;
    double *Low = (double*) calloc(sizeof(double), nsq);
    double *U = (double*) calloc(sizeof(double), nsq);
    
    
    // Now we need to do the LU decomposition
    for (int k = 0; k < n; k++) {
        Low[k + k * n] = 1.0;
        for (int i = k + 1; i < n; i++){
            Low[k + i * n] = C[k + i * n] / C[k + k * n];
            for (int j = k + 1; j < n; j++){
                C[j + i * n] = C[j + i * n] - Low[k + i * n] * C[j + k * n];
            }
        }
        for (int j = k; j < n; j++) {
            U[j + k * n] = C[k + j * n];
        }
    }
    
    // Now we will do Gaussian elimination
    // Solve for xc
    double *ycs = (double*) calloc(sizeof(double), n);
    for (int i = 0; i < (n - 1); i++) {
        for (int j = (i + 1); j < n; j++) {
            double ratio = Low[j + i * n] / Low[i + i * n];
            for (int count = i; count < n; count++) {
                Low[count + j * n] -= ratio * Low[count + i * count];
            }
            D[j] -= (ratio * D[i]);
        }
    }
    ycs[n - 1] = D[n - 1] / Low[(n - 1) + n * (n - 1)];
    for (int i = (n - 2); i >= 0; i--){
        double temp = D[i];
        for (int j = (i + 1); j < n; j++){
            temp -= (Low[j + i * n] * ycs[j]);
        }
        ycs[i] = temp / Low[i + i * n];
    }
    
    //Solve for xc
    for (int i = 0; i < (n - 1); i++){
        for (int j = (i + 1); j < n; j++){
            double ratio = U[j + i * n] / U[i + i * n];
            for (int count = i; count < n; count++){
                U[count + j * n] -= ratio * Low[count + i * count];
            }
            ycs[j] -= (ratio * ycs[i]);
        }
    }
    xcs[n - 1] = ycs[n - 1] / U[(n - 1) + n * (n - 1)];
    for (int i = (n - 2); i >= 0; i--){
        double temp = ycs[i];
        for (int j = (i + 1); j < n; j++){
            temp -= U[j + i * n] * xcs[j];
        }
        xcs[i] = temp / U[i + i * n];
    }
    
    //free everything//
    free(ycs);
    free(U);
    free(Low);
}

void var_convolve(int w, int d, int Q, double* a, int naxes, double* Ref, double* Con) {
    int L = 2 * w + 1;       // kernel axis //
    int nk = pow(2 * w + 1, 2); // number of kernel elements //
    int deg = (d + 1) * (d + 2) / 2; // number of degree elements //
    
    // Now we can do the final convolution //
    double *K = (double*) calloc(sizeof(double), Q);
    int cent = (nk - 1) / 2; //center index
    //do the convolution//
    int ml = 0;
    for (int m = 0; m <= d; m++) {
        for (int l = 0; l <= d - m; l++) {
            for (int i = 0; i < nk; i++) {
                if (i != cent) {
                    K[i + nk * ml] = a[deg * i + ml];
                    K[cent + nk * ml] -= a[deg * i + ml];
                } else {
                    K[i + nk * ml] += a[deg * i + ml];
                }
            }
            ml++;
        }
    }
    int nml = 0;
    for (int m = 0; m <= d;m++) {
        for (int l = 0; l <= d-m;l++) {
            for (int j = 0; j < naxes; j++) {
                for (int i = 0; i < naxes; i++) {
                    for(int nn = 0; nn < L; nn++) {
                        for(int mm = 0; mm < L; mm++) {
                            int ii = i + (mm - w);
                            int jj = j + (nn - w);
                            if (ii >= 0 && ii < naxes && jj >= 0 && jj < naxes) {
                                Con[i + j * naxes] += pow(i, m) * pow(j, l) * Ref[ii + jj * naxes] * K[mm + nn * L + nk * nml];
                            }
                        }
                    }
                }
            }
            nml++;
        }
    }
    free(K);
}

