#pragma once

#include <stdint.h>
#include "common.h"

enum norm_t
{
  NA,
  NM,
  N,
  M
};

typedef struct {
  int64_t * pattern;
  int64_t * di;
  int64_t * dj;
  num_t * cost;
  int64_t np;
  norm_t norm;
} step_pattern;


step_pattern rjIIc(void);

void free_step_pattern(step_pattern * pat);