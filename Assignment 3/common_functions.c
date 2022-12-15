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

float* convolve(float h[], float x[], int lenH, int lenX, int* lenY){ 
  int nconv = lenH+lenX-1;
  (*lenY) = nconv;
  float *y = (float*) malloc(nconv * sizeof(float));

  for(int n=0; n<nconv; n++){
    y[n] = 0;
    for(int k=0; k<lenX; k++){
      if((n-k)>=0 && (n-k)<lenH){
        y[n] += x[k]*h[n-k];
      } 
    }
  }
  return y;
}

float* correlate(float x[], float y[], int lenX, int lenY, int* lenR){
  int ncorr = lenY+lenX-1;
  (*lenR) = ncorr;
  float *r = (float*) malloc(ncorr * sizeof(float));
  float h[lenY];

  for(int i=0; i<lenY; i++){
    h[i] = y[lenY-1-i];
  }

  for(int n=0; n<ncorr; n++){
    r[n] = 0;
    for(int k=0; k<lenX; k++){
      if((n-k)>=0 && (n-k)<lenY){
        r[n] += x[k]*h[n-k];
      } 
    }
  }

  return r;
}

float* downsample(float x[], int lenX, int M){
  float *xd = (float*) malloc((lenX/M) * sizeof(float));
  for(int i=0; i<lenX; i++){
    if(i%M == 0){
      xd[i/M] = x[i];
    }
  }
  return xd;
}

float* upsample(float x[], int lenX, int L){
  float *xu = (float*) malloc((lenX*L) * sizeof(float));
  for(int i=0; i<lenX; i++){
    xu[i*L] = x[i];
    for(int j=1; j<L; j++){
      xu[i*L + j] = 0;
    } 
  }
  return xu;
}

float* decimator(float x[], int xLen, int M, int fc, int fs, int N){
  float *h = lpf(fc, fs, N);
  int yLen;
  float *xf = convolve(h, x, N, xLen, &yLen);
  float *select = (float*) malloc(xLen * sizeof(float));
  for(int i=0; i<xLen; i++){
    select[i] = xf[(N-1)/2 + i];
  }
  float *dx = downsample(select, xLen, M);
  return dx;
}

float* interpolater(float x[], int xLen, int L, int fc, int fs, int N){
  float *dx = upsample(x, xLen, L);
  float *h = lpf(fc, fs, N);
  for(int i=0; i<N; i++){
    h[i] = h[i]*L;
  }
  int yLen;
  float *y = convolve(h, dx, N, xLen*L, &yLen);
  float *out = (float*) malloc(xLen*L * sizeof(float));
  for(int i=0; i<xLen*L; i++){
    out[i] = y[(N-1)/2 + i];
  }
  return out;
}

float* errorVec(float x[], float y[], float len){
  float *out = (float*) malloc(len * sizeof(float));
  for(int i=0; i<len; i++){
    if(x[i] < y[i]){
      out[i] = -1*(x[i] - y[i]);
    }else{
      out[i] = x[i] - y[i];
    }
  }
  return out;
}