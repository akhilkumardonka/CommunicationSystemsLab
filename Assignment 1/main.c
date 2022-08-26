#include <stdio.h>
#include "common_functions.h"
#include <stdlib.h>

int main(int argc, char *argv[])
{
  float h[] = {0.6715, -1.2075, 0.7172, 1.6302, 0.4889, 1.0347, 0.7269, -0.3034, 0.2939, -0.7873, 0.8884, -1.1471, -1.0689, -0.8095, -2.9443};
  float x[] = {0.5377, 1.8339, -2.2588, 0.8622, 0.3188, -1.3077, -0.4336, 0.3426, 3.5784, 2.7694, -1.3499, 3.0349, 0.7254, -0.0631, 0.7147, -0.2050, -0.1241, 1.4897, 1.4090, 1.4172};
  
  int xLen = sizeof(x)/sizeof(x[0]);
  int hLen = sizeof(h)/sizeof(h[0]);
  int yLen;

  // Convolution
  float *y = convolve(h, x, hLen, xLen, &yLen);
  printf("Convolution Output: \n");
  for(int loop = 0; loop < yLen; loop++){
      printf("%f, ", y[loop]);
  }
  printf("\n\n");
  free(y);

  // Correlation
  int rLen;
  float *r = correlate(x, h, xLen, hLen, &rLen);
  printf("Correlation Output: \n");
  for(int loop = 0; loop < rLen; loop++){
      printf("%f, ", r[loop]);
  }
  printf("\n\n");
  free(r);

  // Downsampling
  float d_x[] = {0.3252, -0.7549, 1.3703, -1.7115, -0.1022, -0.2414, 0.3192, 0.3129, -0.8649, -0.0301, -0.1649, 0.6277, 1.0933, 1.1093, -0.8637, 0.0774, -1.2141, -1.1135, -0.0068, 1.5326, -0.7697, 0.3714, -0.2256, 1.1174, -1.0891, 0.0326, 0.5525, 1.1006, 1.5442, 0.0859, -1.4916, -0.7423, -1.0616, 2.3505, -0.6156, 0.7481};
  int lenDx = sizeof(d_x)/sizeof(d_x[0]);
  int M = 2;
  float *dx = downsample(d_x, lenDx, M);
  printf("Downsampled (M = %d) Output: \n", M);
  for(int loop = 0; loop < (lenDx/M); loop++){
      printf("%f, ", dx[loop]);
  }
  printf("\n\n");
  free(dx);

  M = 3;
  float *dx3 = downsample(d_x, lenDx, M);
  printf("Downsampled (M = %d) Output: \n", M);
  for(int loop = 0; loop < (lenDx/M); loop++){
      printf("%f, ", dx3[loop]);
  }
  printf("\n\n");
  free(dx3);

  // Upsampling
  float u_x[] = {0.3252, 1.3703, -0.1022, 0.3192, -0.8649, -0.1649, 1.0933, -0.8637, -1.2141, -0.0068, -0.7697, -0.2256, -1.0891, 0.5525, 1.5442, -1.4916, -1.0616, -0.6156};
  int lenUx = sizeof(u_x)/sizeof(u_x[0]);
  int L = 2;
  float *ux = upsample(u_x, lenUx, L);
  printf("Upsampled (L = %d) Output: \n", L);
  for(int loop = 0; loop < (lenUx*L); loop++){
      printf("%f, ", ux[loop]);
  }
  printf("\n\n");
  free(ux);

  L = 3;
  float *ux3 = upsample(u_x, lenUx, L);
  printf("Upsampled (L = %d) Output: \n", L);
  for(int loop = 0; loop < (lenUx*L); loop++){
      printf("%f, ", ux3[loop]);
  }
  printf("\n\n");
  free(ux3);

  return 0;
}