#include <math.h>
#include <stdlib.h>

#include "common.h"
#include "cost_matrix.h"
#include "step_pattern.h"


void cost_matrix(num_t const *lm, step_pattern const *p, num_t *cm, int64_t *sm,
                 int64_t n, int64_t m) {
  int64_t npats = p->pattern[p->np - 1] + 1;
  int64_t nsteps = p->np;
  num_t *clist = (num_t *)malloc(sizeof(num_t) * npats);

  for (int64_t j = 0; j < m; j++) {
    for (int64_t i = 0; i < n; i++) {
      int64_t ix = IX(i, j);

      if (!isnan(cm[ix])) {
        continue;
      }

      for (int64_t _c = 0; _c < npats; _c++) {
        clist[_c] = NAN;
      }

      for (int64_t s = 0; s < nsteps; s++) {
        int64_t _p = p->pattern[s];
        //printf("%ld, %ld, %ld, %ld\n", i, j, s, ix);
        int64_t ii = i - p->di[s];
        int64_t jj = j - p->dj[s];
        if (ii >= 0 && jj >= 0) { /* address ok? C convention */
          num_t cc = p->cost[s];
          if (near(cc, -1.0, SFLT_EPSILON)) {
            clist[_p] = cm[IX(ii, jj)];
          } else { /* we rely on NAN to propagate */
            clist[_p] += cc * lm[IX(ii, jj)];
          }
        }
      }

      int64_t minc = argmin(clist, npats);
      if (minc > -1) {
        // fprintf(stderr, "%s: %ld\n", "hit", minc);
        cm[ix] = clist[minc];
        sm[ix] = minc;
      }
    }
  }
  free(clist);
}