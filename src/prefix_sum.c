#include <stdint.h>
#include <stdlib.h>

#define MAX_BLOCK_SZ 1024
#define NUM_BANKS 32
#define LOG_NUM_BANKS 5

#ifdef ZERO_BANK_CONFLICTS
#define CONFLICT_FREE_OFFSET(n) ((n) >> NUM_BANKS + (n) >> (2 * LOG_NUM_BANKS))
#else
#define CONFLICT_FREE_OFFSET(n) ((n) >> LOG_NUM_BANKS)
#endif

#include "common.h"

__global__
void gpu_add_block_sums(num_t*  d_out,
  num_t const*  d_in,
  num_t*  d_block_sums,
  uint64_t const num_elems)
{

  uint32_t d_block_sum_val = d_block_sums[blockIdx.x];

  // Simple implementation's performance is not significantly (if at all)
  //  better than previous verbose implementation
  uint32_t cpy_idx = 2 * blockIdx.x * blockDim.x + threadIdx.x;
  if (cpy_idx < num_elems)
  {
    d_out[cpy_idx] = d_in[cpy_idx] + d_block_sum_val;
    if (cpy_idx + blockDim.x < num_elems)
      d_out[cpy_idx + blockDim.x] = d_in[cpy_idx + blockDim.x] + d_block_sum_val;
  }
}

__global__ void gpu_prescan(num_t* d_out, num_t const* d_in,
                            num_t* d_block_sums, uint32_t const len,
                            uint32_t const shmem_sz,
                            uint32_t const max_elems_per_block) {
  // Allocated on invocation
  extern __shared__ unsigned int s_out[];

  int thid = threadIdx.x;
  int ai = thid;
  int bi = thid + blockDim.x;

  // Zero out the shared memory
  // Helpful especially when input size is not power of two
  s_out[thid] = 0;
  s_out[thid + blockDim.x] = 0;
  // If CONFLICT_FREE_OFFSET is used, shared memory
  //  must be a few more than 2 * blockDim.x
  if (thid + max_elems_per_block < shmem_sz)
    s_out[thid + max_elems_per_block] = 0;

  __syncthreads();

  // Copy d_in to shared memory
  // Note that d_in's elements are scattered into shared memory
  //  in light of avoiding bank conflicts
  unsigned int cpy_idx = max_elems_per_block * blockIdx.x + threadIdx.x;
  if (cpy_idx < len) {
    s_out[ai + CONFLICT_FREE_OFFSET(ai)] = d_in[cpy_idx];
    if (cpy_idx + blockDim.x < len)
      s_out[bi + CONFLICT_FREE_OFFSET(bi)] = d_in[cpy_idx + blockDim.x];
  }

  // For both upsweep and downsweep:
  // Sequential indices with conflict free padding
  //  Amount of padding = target index / num banks
  //  This "shifts" the target indices by one every multiple
  //   of the num banks
  // offset controls the stride and starting index of
  //  target elems at every iteration
  // d just controls which threads are active
  // Sweeps are pivoted on the last element of shared memory

  // Upsweep/Reduce step
  int offset = 1;
  for (int d = max_elems_per_block >> 1; d > 0; d >>= 1) {
    __syncthreads();

    if (thid < d) {
      int ai = offset * ((thid << 1) + 1) - 1;
      int bi = offset * ((thid << 1) + 2) - 1;
      ai += CONFLICT_FREE_OFFSET(ai);
      bi += CONFLICT_FREE_OFFSET(bi);

      s_out[bi] += s_out[ai];
    }
    offset <<= 1;
  }

  // Save the total sum on the global block sums array
  // Then clear the last element on the shared memory
  if (thid == 0) {
    d_block_sums[blockIdx.x] =
        s_out[max_elems_per_block - 1 +
              CONFLICT_FREE_OFFSET(max_elems_per_block - 1)];
    s_out[max_elems_per_block - 1 +
          CONFLICT_FREE_OFFSET(max_elems_per_block - 1)] = 0;
  }

  // Downsweep step
  for (int d = 1; d < max_elems_per_block; d <<= 1) {
    offset >>= 1;
    __syncthreads();

    if (thid < d) {
      int ai = offset * ((thid << 1) + 1) - 1;
      int bi = offset * ((thid << 1) + 2) - 1;
      ai += CONFLICT_FREE_OFFSET(ai);
      bi += CONFLICT_FREE_OFFSET(bi);

      unsigned int temp = s_out[ai];
      s_out[ai] = s_out[bi];
      s_out[bi] += temp;
    }
  }
  __syncthreads();

  // Copy contents of shared memory to global memory
  if (cpy_idx < len) {
    d_out[cpy_idx] = s_out[ai + CONFLICT_FREE_OFFSET(ai)];
    if (cpy_idx + blockDim.x < len)
      d_out[cpy_idx + blockDim.x] = s_out[bi + CONFLICT_FREE_OFFSET(bi)];
  }
}

void prefix_sum(num_t* const d_out, const num_t* d_in,
                const uint64_t num_elems) {
  cudaMemset(d_out, 0, num_elems * sizeof(num_t));

  uint32_t block_sz = MAX_BLOCK_SZ / 2;
  uint32_t max_elems_per_block =
      2 * block_sz;  // due to binary tree nature of algorithm

  uint32_t grid_sz = num_elems / max_elems_per_block;

  if (num_elems % max_elems_per_block != 0) {
    grid_sz += 1;
  }

  // Conflict free padding requires that shared memory be more than 2 * block_sz
  uint32_t shmem_sz =
      max_elems_per_block + ((max_elems_per_block - 1) >> LOG_NUM_BANKS);

  // Allocate memory for array of total sums produced by each block
  // Array length must be the same as number of blocks
  num_t* d_block_sums;
  cudaMalloc(&d_block_sums, sizeof(num_t) * grid_sz);
  cudaMemset(d_block_sums, 0, sizeof(num_t) * grid_sz);

  // Sum scan data allocated to each block
  // gpu_sum_scan_blelloch<<<grid_sz, block_sz, sizeof(unsigned int) *
  // max_elems_per_block >>>(d_out, d_in, d_block_sums, numElems);
  gpu_prescan<<<grid_sz, block_sz, sizeof(num_t) * shmem_sz>>>(
      d_out, d_in, d_block_sums, num_elems, shmem_sz, max_elems_per_block);

  // Sum scan total sums produced by each block
  // Use basic implementation if number of total sums is <= 2 * block_sz
  //  (This requires only one block to do the scan)
  if (grid_sz <= max_elems_per_block) {
    num_t* d_dummy_blocks_sums;
    cudaMalloc(&d_dummy_blocks_sums, sizeof(num_t));
    cudaMemset(d_dummy_blocks_sums, 0, sizeof(num_t));
    // gpu_sum_scan_blelloch<<<1, block_sz, sizeof(unsigned int) *
    // max_elems_per_block>>>(d_block_sums, d_block_sums, d_dummy_blocks_sums,
    // grid_sz);
    gpu_prescan<<<1, block_sz, sizeof(num_t) * shmem_sz>>>(
        d_block_sums, d_block_sums, d_dummy_blocks_sums, grid_sz, shmem_sz,
        max_elems_per_block);
    cudaFree(d_dummy_blocks_sums);
  }
  // Else, recurse on this same function as you'll need the full-blown scan
  //  for the block sums
  else {
    num_t* d_in_block_sums;

    cudaMalloc(&d_in_block_sums, sizeof(unsigned int) * grid_sz);
    cudaMemcpy(d_in_block_sums, d_block_sums, sizeof(unsigned int) * grid_sz,
               cudaMemcpyDeviceToDevice);
    prefix_sum(d_block_sums, d_in_block_sums, grid_sz);
    cudaFree(d_in_block_sums);
  }

  gpu_add_block_sums<<<grid_sz, block_sz>>>(d_out, d_out, d_block_sums,
                                            num_elems);

  cudaFree(d_block_sums);
}