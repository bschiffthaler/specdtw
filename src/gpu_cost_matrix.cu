extern "C" {

#include <math.h>
#include <omp.h>
#include <stdlib.h>

#include "common.h"
#include "cost_matrix.h"
#include "gpu_cost_matrix.cuh"
#include "step_pattern.h"

__device__ int64_t gpu_argmin(num_t const *clist, int64_t len) {
  int64_t ii = -1;
  num_t vv = INFINITY;
  for (int64_t i = 0; i < len; i++) {
    if (clist[i] < vv) {
      ii = i;
      vv = clist[i];
    }
  }
  return ii;
}

/*
g[i,j] = min(
     g[i-2,j-1] +     d[i-1,j  ] +     d[i  ,j  ] ,
     g[i-1,j-1] +     d[i  ,j  ] ,
     g[i-1,j-2] +     d[i  ,j-1] + 0 * d[i  ,j  ] ,
  )
*/
__global__ void cuda_kernel_cm_rjIIc(num_t const *lm, num_t *cm, int64_t *sm,
                                     int64_t n, int64_t m, int64_t j) {
  int64_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < n) {
    num_t clist[3];  // 3 patterns
    clist[0] = NAN;
    clist[1] = NAN;
    clist[2] = NAN;
    int64_t ix = IXM(i, j, n);
    num_t dij = lm[ix];
    if ((i - 2) >= 0 && (j - 1) >= 0) {
      clist[0] = cm[IXM(i - 2, j - 1, n)] + lm[IXM(i - 1, j, n)] + dij;
    }
    if ((i - 1) >= 0 && (j - 1) >= 0) {
      clist[1] = cm[IXM(i - 1, j - 1, n)] + dij;
    }
    if ((i - 1) >= 0 && (j - 2) >= 0) {
      num_t c = cm[IXM(i - 1, j - 2, n)] + lm[IXM(i, j - 1, n)];
    }
    int64_t minc = gpu_argmin(clist, 3);
    if (minc > -1) {
      cm[ix] = clist[minc];
      sm[ix] = minc;
    }
  }
}

void gpu_cost_matrix(num_t const *lm, step_pattern const *p, num_t *cm,
                     int64_t *sm, int64_t n, int64_t m) {
  int64_t npats = p->pattern[p->np - 1] + 1;
  int64_t nsteps = p->np;
  num_t *clist = (num_t *)malloc(sizeof(num_t) * npats);

  // First two rows always are CPU
  for (int64_t j = 0; j < 2; j++) {
    for (int64_t i = 0; i < n; i++) {
      int64_t ix = IX(i, j);

      if (!isnan(cm[ix])) {
        continue;
      }

      for (int64_t _c = 0; _c < npats; _c++) {
        clist[_c] = NAN;
      }

      for (int64_t s = 0; s < nsteps; s++) {
        int64_t _p = p->pattern[s];
        // printf("%ld, %ld, %ld, %ld\n", i, j, s, ix);
        int64_t ii = i - p->di[s];
        int64_t jj = j - p->dj[s];
        if (ii >= 0 && jj >= 0) { /* address ok? C convention */
          num_t cc = p->cost[s];
          if (near(cc, -1.0, SFLT_EPSILON)) {
            clist[_p] = cm[IX(ii, jj)];
          } else { /* we rely on NAN to propagate */
            clist[_p] += cc * lm[IX(ii, jj)];
          }
        }
      }

      int64_t minc = argmin(clist, npats);
      if (minc > -1) {
        // fprintf(stderr, "%s: %ld\n", "hit", minc);
        cm[ix] = clist[minc];
        sm[ix] = minc;
      }
    }
  }

  // Now this is pod racing (i.e. where the GPU code starts)
  // Multiple GPUs?
  int devcount;
  cudaGetDeviceCount(&devcount);
  int thread = omp_get_thread_num();
  cudaSetDevice(thread % devcount);

  num_t *d_cm;
  int64_t *d_sm;
  num_t *d_lm;

  cudaMalloc(&d_cm, sizeof(num_t) * n * m);
  cudaMalloc(&d_sm, sizeof(int64_t) * n * m);
  cudaMalloc(&d_lm, sizeof(num_t) * n * m);

  cudaMemcpy(d_cm, cm, sizeof(num_t) * n * m, cudaMemcpyHostToDevice);
  cudaMemcpy(d_sm, sm, sizeof(int64_t) * n * m, cudaMemcpyHostToDevice);
  cudaMemcpy(d_lm, lm, sizeof(num_t) * n * m, cudaMemcpyHostToDevice);

  for (int64_t j = 0; j < m; j++) {
    dim3 block(256);
    dim3 grid(((n + 255) / 256));

    cuda_kernel_cm_rjIIc<<<grid, block>>>(d_lm, d_cm, d_sm, n, m, j);
    cudaStreamSynchronize(0);
  }

  cudaMemcpy(cm, d_cm, sizeof(num_t) * n * m, cudaMemcpyDeviceToHost);
  cudaMemcpy(sm, d_sm, sizeof(int64_t) * n * m, cudaMemcpyDeviceToHost);

  cudaFree(d_cm);
  cudaFree(d_sm);
  cudaFree(d_lm);
  free(clist);
}
}