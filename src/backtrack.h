#pragma once

#include "common.h"
#include "step_pattern.h"

typedef struct {
  int64_t* index1;
  int64_t* index2;
  int64_t* index1s;
  int64_t* index2s;
  int64_t* steps;
  int64_t ii;
  int64_t jj;
  int64_t ss;
} bt_index_t;

void backtrack(step_pattern const* p, num_t const* cm, int64_t const* sm,
               int64_t const n, int64_t const m,  int64_t const jmin,
               bt_index_t * path);