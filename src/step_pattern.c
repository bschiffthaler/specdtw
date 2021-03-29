#include <stdlib.h>
#include <stdint.h>

#include "step_pattern.h"
#include "common.h"

step_pattern rjIIc(void) {

  step_pattern _pat;

  _pat.np = 8;
  _pat.norm = norm_t::N;

  // Straight up index patterns 0-based since we don't need R compatibility
  _pat.pattern = (int64_t *)malloc(sizeof(int64_t) * _pat.np);
  _pat.pattern[0] = 0; _pat.pattern[1] = 0; _pat.pattern[2] = 0;
  _pat.pattern[3] = 1; _pat.pattern[4] = 1; _pat.pattern[5] = 2;
  _pat.pattern[6] = 2; _pat.pattern[7] = 2;

  _pat.di = (int64_t *)malloc(sizeof(int64_t) * _pat.np);
  _pat.di[0] = 2; _pat.di[1] = 1; _pat.di[2] = 0;
  _pat.di[3] = 1; _pat.di[4] = 0; _pat.di[5] = 1;
  _pat.di[6] = 0; _pat.di[7] = 0;

  _pat.dj = (int64_t *)malloc(sizeof(int64_t) * _pat.np);
  _pat.dj[0] = 1; _pat.dj[1] = 0; _pat.dj[2] = 0;
  _pat.dj[3] = 1; _pat.dj[4] = 0; _pat.dj[5] = 2;
  _pat.dj[6] = 1; _pat.dj[7] = 0;

  _pat.cost = (num_t *)malloc(sizeof(num_t) * _pat.np);
  _pat.cost[0] = -1.; _pat.cost[1] = 1.; _pat.cost[2] = 1.;
  _pat.cost[3] = -1.; _pat.cost[4] = 1.; _pat.cost[5] = -1.;
  _pat.cost[6] = 1.; _pat.cost[7] = 0.;

  return _pat;
}

void free_step_pattern(step_pattern * pat) {
  free(pat->pattern);
  free(pat->di);
  free(pat->dj);
  free(pat->cost);
}