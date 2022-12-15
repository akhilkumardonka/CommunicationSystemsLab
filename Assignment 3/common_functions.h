#include <stdio.h>

float* lpf(int fc, int fs, int N); 

float* bpf(int fc1, int fc2, int fs, int N); 

float* convolve(float h[], float x[], int lenH, int lenX, int* lenY);   

float* correlate(float x[], float y[], int lenX, int lenY, int* lenR);

float* downsample(float x[], int lenX, int M);

float* upsample(float x[], int lenX, int L);

float* decimator(float x[], int xLen, int M, int fc, int fs, int N);

float* interpolater(float x[], int xLen, int L, int fc, int fs, int N);

float* errorVec(float x[], float y[], float len);