#include <stdio.h>
#include <stdlib.h>
#include "oisdifference.h"

int main(void) {
    FILE* fp;
    int h = 100, w = 100;
    double *img = (double *)malloc(h * w * sizeof(*img));
    fp = fopen("sample_img_100.txt", "r");
    if (fp == NULL) {
        printf("Problem opening file sample_img_100.txt\n");
        return EXIT_FAILURE;
    }
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            fscanf(fp, "%lf", img + i * w + j);
        }
    }
    fclose(fp);
    double *ref = (double *)malloc(h * w * sizeof(*img));
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
    
    return EXIT_SUCCESS;
}

