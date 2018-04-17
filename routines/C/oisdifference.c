// This code is a version of the Optimal Image Subtraction method detailed in Alard and Lupton 1998// 
// and then reintroduced in Miller 2008. It uses a Delta Function Kernel to solve for the images offset //
//and subtraction. It can be used for either a constant or space-varying kernel depending on how you set the code. //
//I have tried to comment it the best I could so anyone can understand it. If you use this routine, you should cite://
//Alard & Lupton 1998, Alard 2000, Miller+2008, Oelkers+2015, Oelkers & Stassun 2018//

// Include the necessary libraries //
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"
#include <time.h>

char version_number[] = "1.0.1";

void usage(char* exec_name);
void version(char* exec_name);

typedef struct {
    double* data;
    int n;
    int m;
} image;

image fits_get_data(char* filename);
void make_matrix_system(image ref, image sci, int w, int fwhm, int d, int nstars, int* xc, int* yc, double* C, double* D);
void solve_system(int n, double* C, double* D, double* xcs);


int main (int argc, char* argv[])
{
    // Start the clock
    clock_t begin = clock();
    char *exec_name = argv[0];
    
    // Parse command line arguments:
    // Default arguments
    int nstars = -1; // The number of refernce stars in the refstar file
    int w = -1;      // The half-width of the kernel side
    int fwhm = -1;   // The full-width at half-maximum of the stars PSF
    int d = -1;      // Degree of the interpolating polynomial for the variable kernel
    char *reffile = NULL;
    char* scifile = NULL;
    char refstarsfile_default[] = "./refstars.txt";
    char *refstarsfile = refstarsfile_default;
    if ( argc < 2 ) {
        usage(exec_name);
        return EXIT_SUCCESS;
    } else {
        ++argv; // Skip the invocation program name
        --argc;
        while ( argc > 0 )
        {
            if ( !strcmp(*argv, "-nstars") )
            {
                ++argv; // Consume one word
                --argc; // Decrease word counter by 1
                nstars = atoi(*argv);
            }
            else if ( !strcmp(*argv, "-w") ) {
                ++argv; // Consume one word
                --argc; // Decrease word counter by 1
                w = atoi(*argv);
            }
            else if ( !strcmp(*argv, "-ref") ) {
                ++argv; // Consume one word
                --argc; // Decrease word counter by 1
                reffile = *argv;
            }
            else if ( !strcmp(*argv, "-sci") ) {
                ++argv; // Consume one word
                --argc; // Decrease word counter by 1
                scifile = *argv;
            }
            else if ( !strcmp(*argv, "-refstars") ) {
                ++argv; // Consume one word
                --argc; // Decrease word counter by 1
                refstarsfile = *argv;
            }
            else if ( !strcmp(*argv, "-fwhm") ) {
                ++argv; // Consume one word
                --argc; // Decrease word counter by 1
                fwhm = atoi(*argv);
            }
            else if ( !strcmp(*argv, "-d") ) {
                ++argv; // Consume one word
                --argc; // Decrease word counter by 1
                d = atoi(*argv);
            }
            else if ( !strcmp(*argv, "-h") || !strcmp(*argv, "--help") ) {
                usage(exec_name);
                return EXIT_SUCCESS;
            }
            else if ( !strcmp(*argv, "--version") ) {
                version(exec_name);
                return EXIT_SUCCESS;
            }
            else {
                printf("Unexpected Argument: %s\n", *argv);
                usage(exec_name);
                return EXIT_FAILURE;
            }
            ++argv;
            --argc;
        }
    }
    // Check here which variables were not set
    if (nstars == -1) {
        printf("Undefined value for nstars. Exiting.\n");
        usage(exec_name);
        return EXIT_FAILURE;
    }
    if (d == -1) {
        printf("Undefined value for d. Exiting.\n");
        usage(exec_name);
        return EXIT_FAILURE;
    }
    if (w == -1) {
        printf("Undefined value for w. Exiting.\n");
        usage(exec_name);
        return EXIT_FAILURE;
    }
    if (fwhm == -1) {
        printf("Undefined value for fwhm. Exiting.\n");
        usage(exec_name);
        return EXIT_FAILURE;
    }

    // read in the star list //
    int xc[nstars], yc[nstars];
    FILE *fp = fopen(refstarsfile, "r");
    if (fp == NULL) {
        printf("Unable to open refstars.txt file. Exiting.\n");
        return EXIT_FAILURE;
    }
    for (int i = 0; i < nstars; i++) fscanf(fp, "%d %d", xc + i, yc + i);
    fclose(fp);

    
    // Open and read fits files
    image refimg = fits_get_data(reffile);
    double* Ref = refimg.data;
    image sciimg = fits_get_data(scifile);
    double *Sci = sciimg.data;
    int N = (int)(sciimg.n * sciimg.m); // The total number of pixels in the images.
    int naxes = sciimg.n; // This is to fix future references to naxes
    
    // Put a guard in case images are of different shape
    if (refimg.n != sciimg.n || refimg.m != sciimg.m) {
        printf("ERROR: Reference and Science images have different dimensions.\n");
        return EXIT_FAILURE;
    }
 
    int deg = (d + 1) * (d + 2) / 2; // number of degree elements //
    int Q = pow(2 * w + 1, 2) * deg; // size of convolution matrix//
    double *C = (double*) calloc(sizeof(double), Q * Q);
    double *D = (double*) calloc(sizeof(double), Q);
    
    // Now we need to make stamps around each star to find the parameters for the kernel //
    make_matrix_system(refimg, sciimg, w, fwhm, d, nstars, xc, yc, C, D);
    
    double *a = (double*) malloc(sizeof(double) * Q);
    // This will solve the system Cx = D and store x in a
    solve_system(Q, C, D, a);
    free(D);
    free(C);
    
// Now we can do the final convolution // 
    double *Con, *K;
    int nml;
    Con = (double*) malloc(sizeof(double)*N);
    K = (double*) malloc(sizeof(double)*Q);
    
    int L = 2 * w + 1;       // kernel axis //
    int nk = L * L;          // number of kernel elements //

    int cent = (nk - 1)/2;//center index
    for (int i = 0; i < N; i++){
        Con[i]=0;}
    for (int i = 0; i < Q; i++){
        K[i] = 0.0;}
    //do the convolution//
    int ml=0;
    for (int m = 0; m <= d;m++){
        for (int l = 0;l <= d-m;l++){
            for (int i = 0; i < nk; i++){
                if (i != cent){
                    K[i+nk*ml] = a[deg*i+ml];
                    K[cent+nk*ml] =  K[cent+nk*ml] - a[deg*i+ml];
                }
                if (i == cent) {
                    K[i+nk*ml] = K[i+nk*ml]+a[deg*i+ml];
                }
            }
            
            ml++;
        }
    }
    
    nml=0;
    for (int m = 0; m <= d;m++){
        for (int l = 0;l <= d-m;l++){
            for (int j = 0;j < naxes; j++){
                for (int i = 0;i<naxes;i++){
                    for(int nn=0;nn<L; nn++){
                        for(int mm=0;mm<L;mm++){
                            int ii=i+(mm-w);
                            int jj=j+(nn-w);
                            if (ii >=0 && ii < naxes && jj >=0 && jj < naxes){
                                Con[i+j*naxes] = Con[i+j*naxes] + pow(i,m)*pow(j,l)*Ref[ii+jj*naxes]*K[mm+nn*L+nk*nml];
                            }// end of if statment //
                            
                        }// end of nn loop //
                    }// end of mm loop //
                }// end of i loop //
            }// end of j loop //
            nml++;
        }// end of l loop //
    }// end of m loop //

    //free everything//
    free(K); 
    
    //Now we can do the subtraction//
    double *Diff;
    
    Diff = (double*) malloc(sizeof(double)*N);
    
    for (int i = 0; i < N; i++){
        Diff[i] = 0;
        Diff[i] = Sci[i]-Con[i]; // the difference //
    }
    // free everything//
    free(Sci); free(Con);
    
    // Now we need to make a fits file for the differenced image //
    fitsfile *fptd;
    long fpixel, nelements, nax[2];
    double **array;
    int bpix = DOUBLE_IMG;
    char *dfilename, *sciname;
    //printf("working\n");
    sciname = scifile;
    int naxis = 2;
    nax[0] = naxes; nax[1] = naxes;
    
    array = malloc(naxes*sizeof(double*));
    
    for (int i = 0; i< naxes; i++){
        array[i] = (double*) malloc(sizeof(double)*N);}
    for (int i = 1; i < naxes; i++){
        array[i] = array[i-1]+naxes;}
    dfilename = "dimg.fits";
    
    remove(dfilename);
    int status = 0;
    
    fits_create_file(&fptd, dfilename, &status);
    fits_create_img(fptd, bpix, naxis, nax, &status);
    
    for (int j = 0; j < naxes; j++){
        for(int i = 0; i < naxes; i++){
            array[j][i] = Diff[i+j*naxes];}}
    
    fpixel = 1;
    nelements = nax[0]*nax[1];
    
    fits_write_img(fptd, TDOUBLE, fpixel, nelements, array[0], &status);
    
    //free everything//
    free(Diff); free(array);
    
    // Set the Header on the new image
    fitsfile *infptr;      /* pointer to the FITS file, defined in fitsio.h */
    
    char *infilename;
    infilename = scifile;  /* name for existing FITS file   */
    
    char card[FLEN_CARD];
    //printf("%s\n", filename);
    int nkeys1;
    
    status = 0;
    
    /*open the existing FITS file */
    fits_open_file(&infptr, infilename, READWRITE, &status);
    
    //copy the header data
    fits_get_hdrspace(infptr, &nkeys1, NULL, &status);
    
    
    for (int i = 11; i < nkeys1; i++){
        fits_read_record(infptr, i, card, &status);
        //printf("%s\n", card);
        fits_write_record(fptd, card, &status);}
    
    // close everything //
    fits_close_file(infptr, &status) ;
    fits_close_file(fptd, &status);
    
    clock_t end = clock();
    printf("The difference took %f seconds\n", (double) (end-begin)/CLOCKS_PER_SEC);

    free(Ref);

} // end of main file //


image fits_get_data(char* filename) {
    fitsfile *fp;
    int status = 0;
    int bitpix, naxis;
    long naxes[2];
    fits_open_file(&fp, filename, READONLY, &status);
    fits_get_img_param(fp, 2, &bitpix, &naxis, naxes, &status);
    int N = (int)(naxes[0] * naxes[1]); // size of the image //
    long initial_pixel[] = {1, 1};
    double* data = (double*) calloc(N, sizeof(double));
    fits_read_pix(fp, TDOUBLE, initial_pixel, N, 0, data, 0, &status);
    fits_close_file(fp, &status);
    image img = {data, (int)naxes[0], (int)naxes[1]};
    return img;
}

void make_matrix_system(image ref, image sci, int w, int fwhm, int d, int nstars, int* xc, int* yc, double* C, double* D) {
    
    // Now we need to make stamps around each star to find the parameters for the kernel //
    //parameters that fall out from above//
    int naxes = ref.n;
    double* Ref = ref.data;
    double* Sci = sci.data;
    
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
        memset(Kq, 0, nk);
        Kq[q] = 1.0;
        if (q != cent) Kq[cent] = -1.0;
        
        for (int r = 0; r <= d; r++){
            for (int s = 0; s <= d - r; s++){
                for (int n = 0; n < nk; n++){
                    //make the n kernel//
                    memset(Kn, 0, nk);
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
                                        Rs[i + j * stax] = Ref[(i + xcent - fwhm) + (j + ycent - fwhm) * naxes];
                                        Ss[i + j * stax] = Sci[(i + xcent - fwhm) + (j + ycent - fwhm) * naxes];
                                    }
                                }
                                //reinitialize the convolution matrix//
                                memset(CRKn, 0, S);
                                memset(CRKq, 0, S);
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

void usage(char *exec_name) {
    char *exec_basename = strrchr(exec_name, '/') + 1;
    if (exec_basename == NULL) exec_basename = exec_name;
    printf("%s\nAuthor: Ryan Oelkers (c)\n", exec_basename);
    printf("------------------------\n\n");
    printf("usage: %s -fwhm <int> -w <int> -d <int> \
-ref <filename> -sci <filename> [-refstars <filename>] [-h, --help] [--version]\n\n", exec_basename);
    printf("Arguments:\n");
    printf("\t-fwhm: The approximate full-width at half-maximum of the stars PSF in pixels (integer).\n");
    printf("\t-w: the half-width of the kernel to calculate the optimal difference.\n");
    printf("\t-d: Degree of the interpolating polynomial for the variable kernel.\n");
    printf("\t-ref: The reference image path.\n");
    printf("\t-sci: The science image path.\n");
    printf("\t-refstars [optional]: The path to the file with the x, y values \
of the reference stars to estimate the convolution kernel.\n");
    printf("\t\tDefault value is \"refstars.txt\".\n");
    printf("\t-h, --help: Print this help and exit.\n");
    printf("\t--version: Print version information and exit.\n");
    printf("\n");
}


void version(char *exec_name) {
    char *exec_basename = strrchr(exec_name, '/') + 1;
    if (exec_basename == NULL) exec_basename = exec_name;
    printf("%s %s\n", exec_basename, version_number);
}
