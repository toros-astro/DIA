/***********************************************************************************************************************
 *                                                                                                                     *
 *  This code is a version of the Optimal Image Subtraction method detailed in Alard and Lupton 1998                   *
 *  and then reintroduced in Miller 2008. It uses a Delta Function Kernel to solve for the images offset               *
 *  and subtraction. It can be used for either a constant or space-varying kernel depending on how you set the code.   *
 *  I have tried to comment it the best I could so anyone can understand it. If you use this routine, you should cite: *
 *  Alard & Lupton 1998, Alard 2000, Miller+2008, Oelkers+2015, Oelkers & Stassun 2018                                 *
 *                                                                                                                     *
 ***********************************************************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"
#include <time.h>
#include "fitshelper.h"
#include "oisdifference.h"

char version_number[] = "1.0.2";

void usage(char* exec_name);
void version(char* exec_name);

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
    char outputfile_default[] = "diff_img.fits";
    char *outputfile =outputfile_default;
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
            else if ( !strcmp(*argv, "-o") ) {
                ++argv; // Consume one word
                --argc; // Decrease word counter by 1
                outputfile = *argv;
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
    
    // Put a guard in case images are of different shape
    if (refimg.n != sciimg.n || refimg.m != sciimg.m) {
        printf("ERROR: Reference and Science images have different dimensions.\n");
        return EXIT_FAILURE;
    }
    
    // Put a guard in case images are not square
    // This could be relaxed in the future.
    if (refimg.n != refimg.m || sciimg.n != sciimg.m) {
        printf("ERROR: Reference or Science images have not square dimensions.\n");
        return EXIT_FAILURE;
    }
    int naxes = sciimg.n; // This is to fix future references to naxes
    
    //Here do subtraction
    double *Diff = (double*) malloc(sizeof(double) * sciimg.n * sciimg.m);
    perform_subtraction(refimg.n, refimg.m, Ref, Sci, w, fwhm, d, nstars, xc, yc, Diff);
    free(Ref);
    free(Sci);
    
    image diffimg = {Diff, naxes, naxes};
    int success = fits_write_to(outputfile, diffimg, scifile);
    if (success == EXIT_FAILURE) {
        printf("Problem writing diff FITS file.\n");
        return EXIT_FAILURE;
    }
    free(Diff);
    
    clock_t end = clock();
    printf("The difference took %f seconds\n", (double) (end-begin)/CLOCKS_PER_SEC);
    
}

void usage(char *exec_name) {
    char *exec_basename = strrchr(exec_name, '/') + 1;
    if (exec_basename == NULL) exec_basename = exec_name;
    printf("%s\nAuthor: Ryan Oelkers (c)\n", exec_basename);
    printf("------------------------\n\n");
    printf("usage: %s -fwhm <int> -w <int> -d <int> \
-ref <filename> -sci <filename> -nstars <int> [-refstars <filename>] [-o <filename>] [-h, --help] [--version]\n\n", exec_basename);
    printf("Arguments:\n");
    printf("\t-fwhm: The approximate full-width at half-maximum of the stars PSF in pixels (integer).\n");
    printf("\t-w: the half-width of the kernel to calculate the optimal difference.\n");
    printf("\t-d: Degree of the interpolating polynomial for the variable kernel.\n");
    printf("\t-ref: The reference image path.\n");
    printf("\t-sci: The science image path.\n");
    printf("\t-nstars: The number of reference stars.\n");
    printf("\t-refstars [optional]: The path to the file with the x, y values \
of the reference stars to estimate the convolution kernel.\n");
    printf("\t\tDefault value is \"refstars.txt\".\n");
    printf("\t-o [optional]: The path to the subtraction fits file.\n");
    printf("\t\tDefault value is \"diff_img.fits\".\n");
    printf("\t-h, --help: Print this help and exit.\n");
    printf("\t--version: Print version information and exit.\n");
    printf("\n");
}

void version(char *exec_name) {
    char *exec_basename = strrchr(exec_name, '/') + 1;
    if (exec_basename == NULL) exec_basename = exec_name;
    printf("%s %s\n", exec_basename, version_number);
}
