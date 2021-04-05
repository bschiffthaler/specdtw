#pragma once

#include "common.h"
#include "step_pattern.h"

void cost_matrix(num_t const * lm, step_pattern const * p, num_t * out,
  int64_t * sm, int64_t n, int64_t m);