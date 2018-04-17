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

// Start the program to output a float //
int main (int argc, char* argv[])
{
    // Start the clock
    clock_t begin = clock();
    char *exec_name = argv[0];
    int t = 0;
    
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
    
    char listr[1][30];
    char sname[30];
    fitsfile *fptr;
    int status, bitpix, naxis;
    long naxes;// the size of the axes //

    reffile = listr[0]; // the current reference file //
    status = 0;
    
    // get the image dimensions //
    fits_open_file(&fptr, reffile, READONLY, &status);
    fits_get_img_param(fptr, 2, &bitpix, &naxis, &naxes, &status);

    int N; 
    N = naxes*naxes; // size of the image //
    
    // now we can read in the pixel values from the reference //
    float *pixr;
    long fpixr[2];
    double *Ref;
    
    Ref = (double*) malloc(N*sizeof(double));
    pixr = (float*) malloc(N*sizeof(float));
    
    for (int i = 0; i < N; i++){
        Ref[i] = 0.0;} // initialize the reference matrix //
    
    fits_get_img_dim(fptr, &naxis, &status);
    
    fpixr[0] = 1.0;fpixr[1]=1.0;
    int j = 0;
    
    // read in the pixel values //
    for (fpixr[1] = 1; fpixr[1] <= naxes; fpixr[1]++){
        fits_read_pix(fptr, TFLOAT, fpixr, naxes, 0, pixr, 0, &status);
        for (int i = 0; i < naxes; i++){
            Ref[i+j*naxes] = (double) pixr[i];}
        j++;}
    fits_close_file(fptr, &status); free(pixr);

    char listn[1][30];
        
    sprintf(sname, "./img.txt");
    fp = fopen(sname,"r");
        
    for (int i = 0; i < 1; i++){
        fscanf(fp, "%s\n", listn[i]);} // read in file names only //
    fclose(fp);
    
    // now we can read in each file to difference against the reference frame //
    fitsfile *fpts;
    float *pixs;
    long fpixs[2];
    double *Sci;
    
    Sci = (double*) malloc(N*sizeof(double));
    pixs = (float*) malloc(N*sizeof(float));
    
    for (int i = 0; i < N; i++){
        Sci[i] = 0.0;}
    status = 0;
    
    scifile = listn[0];
    
    fits_open_file(&fpts, scifile, READONLY, &status);
    fits_get_img_dim(fpts, &naxis, &status);
    
    fpixs[0] = 1.0;fpixs[1] = 1.0;
    j = 0;
    
    // read in the pixel values //
    for (fpixs[1] = 1.0; fpixs[1] <= naxes; fpixs[1]++){
        fits_read_pix(fpts, TFLOAT, fpixs, naxes, 0, pixs, 0, &status);
        for (int i = 0; i < naxes; i++){
            Sci[i+j*naxes] = (double)pixs[i];
            //printf("%f\n", Sci[i+j*naxes]);
        }
        j++;}
    fits_close_file(fpts, &status); free(pixs);
    
    // Now we need to make stamps around each star to find the parameters for the kernel //
    
    double *Rs, *Ss, *C, *D, *CRKn, *CRKq, *Kn, *Kq;
    int n, q, mm, nn, l, r, s, m, ii, jj, xcent, ycent;
    int stax, L, mr, ls, ml,qrs, nk, S, deg, P, Q, Qtwo, cent;
    
    //parameters that fall out from above//
    L = 2*w + 1; // kernel axis //
    nk = L*L; // number of kernel elements //
    stax = 2*fwhm + 1; // size of star stamps //
    S = stax*stax; // number of stamp elements //
    deg = (0.5)*(d+1)*(d+2); // number of degree elements //
    Q = nk*deg;//size of D matrix//
    Qtwo = Q*Q;//size of C matrix//
    P = nstars; // number of star stamps to use //
    cent = (nk-1)/2;//center of the kernel//
    
    printf("The kernel size is %d x %d, the polynomial degree is %d, and %d stars were used.\n", L, L, d, P);
    
    Rs = (double*) malloc(sizeof(double)*S);
    Ss = (double*) malloc(sizeof(double)*S);
    
    CRKq = (double*) malloc(sizeof(double)*S);
    CRKn = (double*) malloc(sizeof(double)*S);

    C = (double*) malloc(sizeof(double)*Qtwo);
    D = (double*) malloc(sizeof(double)*Q);

    Kn = (double*) malloc(sizeof(double)*nk);
    Kq = (double*) malloc(sizeof(double)*nk);
    
    for (int i = 0; i< Qtwo; i++){
        C[i] = 0;}
    
    // now we need to solve for the kernel paramters //
    
    qrs=0;//initialize the qrs step//
    for (q=0; q<nk; q++){
        
        //make the q kernel//
        for (int i=0;i<nk;i++){
            Kq[i]=0;}
        Kq[q]=1.0;
        if (q !=cent){
            Kq[cent]=-1.0;
        }
        
        for (r = 0; r <= d; r++){
            for (s = 0; s <= d-r;s++){
                
                for (n = 0; n < nk;n++){
                    
                    //make the n kernel//
                    for (int i=0; i < nk;i++){
                        Kn[i]=0;}
                    Kn[n]=1.0;
                    if (n !=cent){
                        Kn[cent]=-1.0;
                    }
                    ml=0;//initialize the ml step//
                    for (m = 0; m <= d;m++){
                        for (l = 0; l <= d-m; l++){
                            D[qrs]=0;//ensure D is only calculated once for each Q//
                            
                            for (int k = 0; k < P; k++){
                                
                                xcent = xc[k];//x coordinate of stamp center//
                                ycent = yc[k];//y coordinate of stamp center//
                                
                                //make the star stamps//
                                for (int i=0;i<stax;i++){
                                    for(j=0;j<stax;j++){
                                        Rs[i+j*stax] = Ref[(i+xcent-fwhm)+(j+ycent-fwhm)*naxes];
                                        Ss[i+j*stax] = Sci[(i+xcent-fwhm)+(j+ycent-fwhm)*naxes];
                                        //printf("Ref value is %f, Sci value is %f\n", Rs[i+j*stax], Ss[i+j*stax]);
                                    }//end of i loop//
                                }//end of j loop//
                               
                                //reinitialize the convolution matrix//
                                for (int i = 0; i < S; i++){
                                    CRKn[i]=0;CRKq[i]=0;}
                                
                                //now we do the convolution for n and q//
                                for (int i=0; i<stax;i++){
                                    for(int j=0;j<stax;j++){
                                        for (mm=0;mm<L;mm++){
                                            for(nn=0;nn<L;nn++){
                                                ii=i+(mm-w);//index of convolution//
                                                jj=j+(nn-w);//index of convolution//
                                                if (ii>=0 && ii<stax && jj>=0 && jj < stax){
                                                    CRKn[i+j*stax]=CRKn[i+j*stax]+Rs[ii+jj*stax]*Kn[mm+nn*L];
                                                    CRKq[i+j*stax]=CRKq[i+j*stax]+Rs[ii+jj*stax]*Kq[mm+nn*L];
                                                }//end of if statement//
                                            }//end of nn loop//
                                        }// end of mm loop//
                                    }//end of j loop//
                                }//end of i loop//
                                
                                mr = m+r; ls = l+s;//exponents for polynomial approximation//
                                
                                //now we need to fill in C//
                                for (int i=0;i<S;i++){
                                    C[n*deg+ml+qrs*Q]=C[n*deg+ml+qrs*Q]+pow(xcent,mr)*pow(ycent,ls)*CRKn[i]*CRKq[i];
                                    //printf("C is %f\n", C[]);
                                }//end of C loop//
                                
                                //now we need to fill in D//
                                for (int i = 0; i < S;i++){
                                    D[qrs]=D[qrs]+pow(xcent,r)*pow(ycent,s)*Ss[i]*CRKq[i];
                                    //printf("D is %f\n", D[qrs]);
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
    free(CRKn); free(CRKq); free(Kn); free(Kq);free(Ss); free(Rs);
    //for (i=0; i < Qtwo; i++){
      //  printf("%f\n", C[i]);}
    
    double *Low, *U, *xcs, *ycs, *a;
    int count;
    double ratio, temp;
    
    Low = (double*) malloc(sizeof(double)*Qtwo);
    U =  (double*) malloc(sizeof(double)*Qtwo);
    xcs = (double*) malloc(sizeof(double)*Q);
    ycs = (double*) malloc(sizeof(double)*Q);
    a = (double*) malloc(sizeof(double)*Q);
    
   // Now we need to do the LU decomposition
   
    for (int i = 0; i < Qtwo; i++){
        Low[i] = 0.0; U[i] = 0.0;}
        
    for (int k = 0; k < Q; k++){
        Low[k+k*Q] = 1.0;
        for (int i = k + 1; i < Q; i++){
            Low[k+i*Q] = C[k+i*Q]/C[k+k*Q];
            for (j = k + 1; j < Q; j++){
                C[j+i*Q] = C[j+i*Q] - Low[k+i*Q]*C[j+k*Q];
            }
        }
        
        for (j = k; j < Q; j++){
            U[j+k*Q] = C[k+j*Q];
        }
    }
        
    // Now we will do Gaussian elimination
    // Solve for yc
    ratio = 0; temp = 0;
    for(int i=0; i<Q;i++){
        ycs[i]=0;xcs[i]=0;}
    for (int i = 0; i < (Q-1); i++){
        for (j = (i+1); j<Q; j++){
            ratio = Low[j+i*Q] / Low[i+i*Q];
            for (count = i; count < Q; count++){
                Low[count+j*Q] -= (ratio*Low[count+i*count]);}
            D[j] -= (ratio*D[i]);}}

    ycs[Q-1] = D[Q-1] / Low[(Q-1)+Q*(Q-1)]; 

    for (int i = (Q-2); i >= 0; i--){
        temp = D[i]; 
        for (j = (i+1); j < Q; j++){
            temp -= (Low[j+i*Q]*ycs[j]);}
        ycs[i] = temp / Low[i+i*Q];}

        
    //Solve for xc
    for (int i = 0; i < (Q-1); i++){
        for (j = (i+1); j<Q; j++){
            ratio = U[j+i*Q] / U[i+i*Q];
            for (count = i; count < Q; count++){
                U[count+j*Q] -= (ratio*Low[count+i*count]);}
            ycs[j] -= (ratio*ycs[i]);}}
    
    xcs[Q-1] = ycs[Q-1] / U[(Q-1)+Q*(Q-1)]; 
    
    for (int i = (Q-2); i >= 0; i--){
        temp = ycs[i];
        for (j = (i+1); j < Q; j++){
            temp -= (U[j+i*Q]*xcs[j]);}
        xcs[i] = temp / U[i+i*Q];}

    for (int i = 0; i < Q; i++){
        a[i] = 0;
        a[i] =  xcs[i];
    }
    
    //free everything//
    free(xcs); free(ycs); free(U); free(Low); free(D); free(C); 
    
// Now we can do the final convolution // 
    double *Con, *K;
    int nml;
    Con = (double*) malloc(sizeof(double)*N);
    K = (double*) malloc(sizeof(double)*Q);
    
    cent = (nk-1)/2;//center index
    for (int i = 0; i < N; i++){
        Con[i]=0;}
    for (int i = 0; i < Q; i++){
        K[i] = 0.0;}
    //do the convolution//
    ml=0;
    for (m = 0; m <= d;m++){
        for (l = 0;l <= d-m;l++){
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
    
    //printf("%f\n", a[n+nk*ml], m,l,ml);
    nml=0;
    for (m = 0; m <= d;m++){
        for (l = 0;l <= d-m;l++){
            for (int j = 0;j < naxes; j++){
                for (int i = 0;i<naxes;i++){
                    for(nn=0;nn<L; nn++){
                        for(mm=0;mm<L;mm++){
                            ii=i+(mm-w);
                            jj=j+(nn-w);
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
    sciname = listn[t];
    naxis = 2;
    nax[0] = naxes; nax[1] = naxes;
    
    array = malloc(naxes*sizeof(double*));
    
    for (int i = 0; i< naxes; i++){
        array[i] = (double*) malloc(sizeof(double)*N);}
    for (int i = 1; i < naxes; i++){
        array[i] = array[i-1]+naxes;}
    dfilename = "dimg.fits";
    
    //sprintf(dfilename, "./d%s",sciname);
    remove(dfilename);
    status = 0;
    
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
    infilename = listn[t];  /* name for existing FITS file   */
    
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
