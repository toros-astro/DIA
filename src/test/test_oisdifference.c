#include <stdio.h>
#include <stdlib.h>
#include "oisdifference.h"

int main(void) {
    FILE* fp;
    int h = 100, w = 100;
    double *sci = (double *)malloc(h * w * sizeof(*sci));
    fp = fopen("sample_img_100.txt", "r");
    if (fp == NULL) {
        printf("Problem opening file sample_img_100.txt\n");
        return EXIT_FAILURE;
    }
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            fscanf(fp, "%lf", sci + i * w + j);
        }
    }
    fclose(fp);
    double *ref = (double *)malloc(h * w * sizeof(*ref));
    fp = fopen("sample_ref_100.txt", "r");
    if (fp == NULL) {
        printf("Problem opening file sample_ref_100.txt\n");
        return EXIT_FAILURE;
    }
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            fscanf(fp, "%lf", ref + i * w + j);
        }
    }
    fclose(fp);
    int kw = 3;
    int fwhm = 11;
    int d = 0;
    int nstars = 10;
    int xc[nstars], yc[nstars];
    double *subt = (double *)malloc(h * w * sizeof(*subt));
    
    fp = fopen("refstars.txt", "r");
    if (fp == NULL) {
        printf("Problem opening file refstars.txt\n");
        return EXIT_FAILURE;
    }
    for (int i = 0; i < nstars; i++) {
        fscanf(fp, "%d %d", xc + i, yc + i);
    }
    fclose(fp);
    
    image sci_img = {sci, h, w};
    image ref_img = {ref, h, w};
    perform_subtraction(ref_img, sci_img, kw, fwhm, d, nstars, xc, yc, subt);
    
    double norm = 0;
    for (int i = 0; i < h * w; i++) norm += subt[i] * subt[i];
    norm = sqrt(norm);
    
    double refnorm = 0;
    for (int i = 0; i < h * w; i++) refnorm += ref[i] * ref[i];
    refnorm = sqrt(refnorm);
    
    free(ref);
    free(sci);
    free(subt);
    
    double perc_error = norm / refnorm;
    if (perc_error > 0.1) return EXIT_FAILURE;
    return EXIT_SUCCESS;
}

