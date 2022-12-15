#include <stdio.h>
#include "common_functions.h"
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]){

    // input
    int f0 = 100;
    int f1 = 200;
    int f2 = 300;
    int fs = 2400;
    int sampSize = 80;
    const double pi = 22.0 / 7.0;

    float *sig = (float*) malloc(sampSize * sizeof(float));

    for(int n=0; n<sampSize; n++){
        sig[n] = sin((2*pi*f0*n)/fs) + 0.5*sin((2*pi*f1*n)/fs) + 0.6*sin((2*pi*f2*n)/fs);
    }

    printf("\n\nSignal Vector: \n");
    for(int loop = 0; loop < sampSize; loop++){
        printf("%f, ", sig[loop]);
    }

    // Task 1
    int fc = 600;
    int N = 51;
    int M = 2;
    int L = 2;

    float *deciOut = decimator(sig, sampSize, M, fc, fs, N);

    // printf("\nDecimater (M=2) Output: \n");
    // for(int loop = 0; loop < sampSize/M; loop++){
    //     printf("%f, ", deciOut[loop]);
    // }
    
    float *yOut = interpolater(deciOut, sampSize/M, L, fc, fs, N);

    printf("\n\nOutput (M=L=2): \n");
    for(int loop = 0; loop < sampSize; loop++){
        printf("%f, ", yOut[loop]);
    }

    float *errors = errorVec(sig, yOut, sampSize);

    printf("\n\nError Vector (M=L=2): \n");
    float errSum = 0;
    for(int loop = 0; loop < sampSize; loop++){
        printf("%f, ", errors[loop]);
        errSum += errors[loop];
    }

    printf("\n\nAvg Error (M=L=2): %f", errSum/sampSize);

    // Task 2
    fc = 300;
    N = 51;
    M = 4;
    L = 4;

    float *deciOut4 = decimator(sig, sampSize, M, fc, fs, N);

    // printf("\nDecimater (M=4) Output: \n");
    // for(int loop = 0; loop < sampSize/M; loop++){
    //     printf("%f, ", deciOut[loop]);
    // }
    
    float *yOut4 = interpolater(deciOut4, sampSize/M, L, fc, fs, N);

    printf("\n\nOutput (M=L=4): \n");
    for(int loop = 0; loop < sampSize; loop++){
        printf("%f, ", yOut4[loop]);
    }

    float *errors4 = errorVec(sig, yOut4, sampSize);

    printf("\n\nError Vector (M=L=4): \n");
    errSum = 0;
    for(int loop = 0; loop < sampSize; loop++){
        printf("%f, ", errors4[loop]);
        errSum += errors4[loop];
    }

    printf("\n\nAvg Error (M=L=4): %f", errSum/sampSize);

    return 0;
}