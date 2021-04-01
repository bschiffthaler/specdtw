#pragma once

#include <stdint.h>

#include "common.h"

__global__ void cuda_kernel_euclidean(uint64_t N, num_t *sp_a, num_t *sp_b,
                                      num_t *out, int64_t n, int64_t m);

void euclidean(num_t const *sp_a, num_t const *sp_b, num_t *out, int64_t n, int64_t m);