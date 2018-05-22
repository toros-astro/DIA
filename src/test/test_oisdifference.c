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
    
    int kw = 1; //3, 5
    int fwhm = 0;
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
    
    printf("Norm: %f\n", norm);
    
    return EXIT_SUCCESS;
}

