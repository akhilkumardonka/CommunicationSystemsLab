#include "common_functions.h"
#include <stdlib.h>

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