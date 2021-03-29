#pragma once

#include "common.h"
#include "step_pattern.h"

bool almost_equal(num_t A, num_t B);
int64_t argmin(num_t const *clist, int64_t len);

void cost_matrix(num_t const * lm, step_pattern * p, num_t * out,
  int64_t * sm, int64_t n, int64_t m);