#include <stdio.h>

extern float* convolve(float h[], float x[], int lenH, int lenX, int* lenY);   

extern float* correlate(float x[], float y[], int lenX, int lenY, int* lenR);

extern float* downsample(float x[], int lenX, int M);

extern float* upsample(float x[], int lenX, int L); 

