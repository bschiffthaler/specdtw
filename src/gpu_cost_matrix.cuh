#pragma once

#include <stdint.h>
#include "common.h"
#include "step_pattern.h"

void gpu_cost_matrix(num_t const *lm, step_pattern const *p, num_t *cm,
                     int64_t *sm, int64_t n, int64_t m);