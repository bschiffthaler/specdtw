#include "backtrack.h"

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "step_pattern.h"

void backtrack(step_pattern const* p, int64_t const* sm, int64_t const n,
               int64_t const jmin, bt_index_t* path) {
  int64_t npat = p->pattern[p->np - 1] + 1;

  int64_t i = n - 1;
  int64_t j = jmin;

  /*
    What follows here is a really convoluted way to precompute the steps
    in reverse row order. Essentially all it does is store the corresponding
    indices of non-null deltas in an array.
  */

  // Make an array that counts how many non-null deltas each pattern has
  int64_t* pat_cnt = (int64_t*)malloc(sizeof(int64_t) * npat);
  for (int64_t _i = 0; _i < npat; _i++) {
    pat_cnt[_i] = 0;
  }
  // Now actually count them
  for (int64_t _i = 0; _i < p->np; _i++) {
    if (p->di[_i] != 0 || p->dj[_i] != 0) {
      int64_t pp = p->pattern[_i];
      pat_cnt[pp]++;
    }
  }

  // Now make an array of pointers. Each will hold the indices of the
  // non null deltas
  int64_t** pat_map = (int64_t**)malloc(sizeof(int64_t**) * npat);
  for (int64_t _i = 0; _i < npat; _i++) {
    pat_map[_i] = (int64_t*)malloc(sizeof(int64_t) * pat_cnt[_i]);
  }
  // Go through the steps
  for (int64_t _i = 0; _i < p->np; _i++) {
    // Check if not null
    if (p->di[_i] != 0 || p->dj[_i] != 0) {
      int64_t pp = p->pattern[_i];
      // Add index to array starting from the back
      pat_map[pp][pat_cnt[pp] - 1] = _i;
      pat_cnt[pp]--;
    }
  }
  // Fill count again
  for (int64_t _i = 0; _i < p->np; _i++) {
    if (p->di[_i] != 0 || p->dj[_i] != 0) {
      int64_t pp = p->pattern[_i];
      // fprintf(stderr, "Inc pp: %ld\n", pp);
      pat_cnt[pp]++;
    }
  }

  // for (int64_t _i = 0; _i < npat; _i++) {
  //   for (int64_t _j = 0; _j < pat_cnt[_i]; _j++) {
  //     fprintf(stderr, "P: %ld, I: %ld, N: %ld\n", _i, pat_map[_i][_j],
  //     pat_cnt[_i]);
  //   }
  // }

  // Hold end positions for arrays
  uint64_t ss = 0;
  uint64_t ii = 0;
  uint64_t jj = 0;

  // add endpoints
  path->index1[ii] = i;
  path->index1s[jj] = i;
  path->index2[ii++] = j;
  path->index2s[jj++] = j;

  while (true) {
    if (i == 0 && j == 0) {
      // fprintf(stderr, "Hit break @ [%ld, %ld]\n", i, j);
      break;
    }

    int64_t s = sm[IXM(i, j, n)];

    if (s == LONG_MIN) {
      // fprintf(stderr, "Hit break @ [%ld, %ld] for s=%ld\n", i, j, s);
      break;
    }

    path->steps[ss++] = s;
    int64_t ns = pat_cnt[s];

    for (int64_t k = 0; k < ns; k++) {
      path->index1[ii] = i - p->di[pat_map[s][k]];
      path->index2[ii] = j - p->dj[pat_map[s][k]];
      ii++;
    }

    i = i - p->di[pat_map[s][ns - 1]];
    j = j - p->dj[pat_map[s][ns - 1]];

    path->index1s[jj] = i;
    path->index2s[jj] = j;
    jj++;
  }
  path->ii = ii;
  path->jj = jj;
  path->ss = ss;

  for (int64_t _i = 0; _i < npat; _i++) {
    free(pat_map[_i]);
  }
  free(pat_cnt);
  free(pat_map);
}

void free_bt_index(bt_index_t * f) {
  free(f->index1);
  free(f->index1s);
  free(f->index2);
  free(f->index2s);
  free(f->steps);
}