#include <stdio.h>
#include "common_functions.h"
#include <stdlib.h>

int main(int argc, char *argv[]){

    // Task 1
    int fc = 400;
    int fs = 1600;
    int N = 39;
    float *out = lpf(fc, fs, N);
    printf("\nLow Pass Filter (wc = pi/2) Output: \n");
    for(int loop = 0; loop < N; loop++){
        printf("%f, ", out[loop]);
    }

    // Task 2
    fc = 400;
    fs = 3200;
    N=39;
    float *out2 = lpf(fc, fs, N);
    printf("\n\nLow Pass Filter (wc = pi/4) Output: \n");
    for(int loop = 0; loop < N; loop++){
        printf("%f, ", out2[loop]);
    }

    // Task 3
    int fc1 = 500;
    int fc2 = 1200;
    fs = 6000;
    N=39;
    float *out3 = bpf(fc1, fc2, fs, N);
    printf("\n\nBand Pass Filter (fs = 6000) Output: \n");
    for(int loop = 0; loop < N; loop++){
        printf("%f, ", out3[loop]);
    }
}