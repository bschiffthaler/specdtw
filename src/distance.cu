
#include <stdio.h>

#include "distance.cuh"
#include "common.h"

/*
  Calculate:
    sqrt(a^2 + b^2 + 2ab)
  for two vectors, which equates the euclidean distance of the outer
  indices. In R: `abs(outer(a, b, "-"))`
*/
__global__ void cuda_kernel_euclidean(uint64_t N, num_t *sp_a, num_t *sp_b,
                                      num_t *out, int64_t n, int64_t m) {
  int64_t row = blockIdx.x * blockDim.x + threadIdx.x;
  int64_t col = blockIdx.y * blockDim.y + threadIdx.y;

  // Shift out index one row downwards.
  int64_t ix = (col * (n + 1) + row + 1);

  if (row < n && col < m) {
    //printf("[%ld, %ld, %ld, %d, %d, %d]\n", row, col, ix, blockIdx.x, blockDim.x, threadIdx.x);
    out[ix] = fabs(sp_a[row] - sp_b[col]);
    // if (isnan(out[ix])) {
    //   printf("%e, %e, %ld, %ld, %ld\n", sp_a[row], sp_b[col], row, col, ix);
    // }
  }
}

void euclidean(num_t *sp_a, num_t *sp_b, num_t *out, int64_t n, int64_t m) {
  
  num_t *d_sp_a;
  num_t *d_sp_b;
  num_t *d_out;
  num_t *d_pow_a;
  num_t *d_pow_b;
  cudaMalloc(&d_sp_a, sizeof(num_t) * n);
  cudaMalloc(&d_sp_b, sizeof(num_t) * m);
  cudaMalloc(&d_out, sizeof(num_t) * (n + 1) * m);
  cudaMalloc(&d_pow_a, sizeof(num_t) * n);
  cudaMalloc(&d_pow_b, sizeof(num_t) * m);

  cudaMemcpy(d_sp_a, sp_a, sizeof(num_t) * n, cudaMemcpyHostToDevice);
  cudaMemcpy(d_sp_b, sp_b, sizeof(num_t) * m, cudaMemcpyHostToDevice);

  // TODO: Think about launchparams
  dim3 block(16, 16);
  dim3 grid(((n + 15) / 16), ((m + 15) / 16)); 

  // Run kernel
  cuda_kernel_euclidean<<<grid, block>>>(
      n * m, d_sp_a, d_sp_b, d_out, n, m);

  cudaDeviceSynchronize();

  cudaMemcpy(out, d_out, sizeof(num_t) * (n + 1) * m, cudaMemcpyDeviceToHost);

  cudaFree(d_sp_a);
  cudaFree(d_sp_b);
  cudaFree(d_out);
}