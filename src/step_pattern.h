#pragma once

#include <stdint.h>
#include "common.h"

typedef enum
{
  NA,
  NM,
  N,
  M
} norm_t;

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
void norm_row(step_pattern const *p, int64_t const n, int64_t const m,
              num_t *row);