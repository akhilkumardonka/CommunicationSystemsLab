#include "common_functions.h"
#include <stdlib.h>
#include <math.h>

const double pi = 22.0 / 7.0;

float* hammingWindow(int N){
    float *out = (float*) malloc(N * sizeof(float));
    for(int i=0; i<N; i++){
        out[i] = 0.54 - 0.46*cos((2*pi*i)/(N-1));
    }
    return out;
}

float* lpf(int fc, int fs, int N){

    // getting hamming window
    float *window = hammingWindow(N);

    // designing Low Pass Filter 
    float wc = (2*pi*fc)/fs;
    float *out = (float*) malloc(N * sizeof(float));
    for(int i=0; i<=N; i++){
        int n = i - (N-1)/2;
        if(n==0){
            out[i] = (wc/pi) * window[i];
        }else{
            out[i] = (sin(wc*n)/(pi*n)) * window[i];
        }
    }

    return out;
}

float* bpf(int fc1, int fc2, int fs, int N){

    // getting hamming window
    float *window = hammingWindow(N);

    // designing Band Pass Filter 
    float wc1 = (2*pi*fc1)/fs;
    float wc2 = (2*pi*fc2)/fs;
    float *out = (float*) malloc(N * sizeof(float));
    for(int i=0; i<=N; i++){
        int n = i - (N-1)/2;
        if(n==0){
            out[i] = ((wc2 - wc1)/pi) * window[i];
        }else{
            out[i] = (sin(wc2*n)/(pi*n) - sin(wc1*n)/(pi*n)) * window[i];
        }
    }
    return out;
}