#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "backtrack.h"
#include "common.h"
#include "cost_matrix.h"
#include "gpu_cost_matrix.cuh"
#include "distance.cuh"
#include "step_pattern.h"

int main(int argc, char const *argv[]) {
  clock_t program_start = clock();

  if (argc != 6) {
    fprintf(stderr,
      "usage: %s <query> <reference> <index_out> <indexs_out> <dist_out>\n",
     argv[0]);
    return 1;
  }


  FILE *index_out = fopen(argv[3], "w");
  FILE *indexs_out = fopen(argv[4], "w");
  FILE *dist_out = fopen(argv[5], "w");

  if (! index_out) {
    fprintf(stderr, "Error opening output: %s\n", argv[3]);
    return 1;
  }
  if (! indexs_out) {
    fprintf(stderr, "Error opening output: %s\n", argv[4]);
    return 1;
  }
  if (! dist_out) {
    fprintf(stderr, "Error opening output: %s\n", argv[5]);
    return 1;
  }

  fprintf(stderr, "%s\n", "Reading input spectra...");
  clock_t begin_read = clock();
  fvec spec_a = make_fvec(argv[1]);
  fvec spec_b = make_fvec(argv[2]);
  if (spec_a.size == 0 || spec_b.size == 0) {
    return 1;
  }
  clock_t end_read = clock();
  fprintf(stderr, "Time to read: " num_fmt "\n",
          (num_t)(end_read - begin_read) / CLOCKS_PER_SEC);

  // Check input dims
  if ((spec_a.size + 1) * spec_b.size > LONG_MAX) {
    fprintf(stderr, "%s\n", "Input dimensions too large");
  }
  int64_t n = spec_a.size;
  int64_t m = spec_b.size;

  fprintf(stderr, "Matrix dims: [%ld, %ld]\n", n, m);

  fprintf(stderr, "%s\n", "Calculating distance matrix lm...");
  clock_t begin_lm = clock();
  // Extra for null row
  num_t *lm = (num_t *)malloc(sizeof(num_t) * (n + 1) * m);
  euclidean(spec_a.data, spec_b.data, lm, n, m);
  for (int64_t i = 0; i < m; i++) {
    lm[(i * (n + 1))] = 0.f;
  }
  clock_t end_lm = clock();
  fprintf(stderr, "Time to calculate lm: " num_fmt "\n",
          (num_t)(end_lm - begin_lm) / CLOCKS_PER_SEC);


  // Allocate cost matrix. Since we always do open.begin we allocate an
  // additional null row
  num_t *cm = (num_t *)malloc(sizeof(num_t) * (n + 1) * m);
  for (int64_t i = 0; i < ((n + 1) * m); i++) {
    cm[i] = NAN;
  }
  for (int64_t i = 0; i < m; i++) {
    cm[(i * (n + 1))] = 0.f;
  }

  // Allocate step matrix. Initialize to INTEGER_MIN
  int64_t *sm = (int64_t *)malloc(sizeof(int64_t) * (n + 1) * m);
  for (int64_t i = 0; i < ((n + 1) * m); i++) {
    sm[i] = LONG_MIN;
  }

  step_pattern p = rjIIc();

  fprintf(stderr, "%s\n", "Calculating cost matrix cm...");
  clock_t begin_cm = clock();
#ifdef CPU_ONLY
  cost_matrix(lm, &p, cm, sm, n + 1, m);
#else
  gpu_cost_matrix(lm, &p, cm, sm, n + 1, m);
#endif
  clock_t end_cm = clock();
  fprintf(stderr, "Time to calculate cm: " num_fmt "\n",
          (num_t)(end_cm - begin_cm) / CLOCKS_PER_SEC);

  num_t *last_row = (num_t *)malloc(sizeof(num_t) * m);
  for (int64_t i = 0; i < m; i++) {
    last_row[i] = cm[IXM(n, i, n + 1)];
  }

  norm_row(&p, n, m, last_row);

  // open.end=TRUE
  int64_t jmin = argmin(last_row, m);
  num_t distance = cm[IXM(n, jmin, n + 1)];
  num_t normalized_distance = last_row[jmin];

  fprintf(stderr, "Distance:" num_fmt "\n", distance);
  fprintf(stderr, "Normalized distance:" num_fmt "\n", normalized_distance);

  fprintf(dist_out, num_fmt " " num_fmt "\n", distance, normalized_distance);
  fclose(dist_out);

  bt_index_t path;
  path.index1 = (int64_t *)malloc(sizeof(int64_t) * (n + m + 1));
  path.index1s = (int64_t *)malloc(sizeof(int64_t) * (n + m + 1));
  path.index2 = (int64_t *)malloc(sizeof(int64_t) * (n + m + 1));
  path.index2s = (int64_t *)malloc(sizeof(int64_t) * (n + m + 1));
  path.steps = (int64_t *)malloc(sizeof(int64_t) * (n + m + 1));

  fprintf(stderr, "%s\n", "Backtracking...");
  clock_t begin_bt = clock();
  backtrack(&p, sm, n + 1, jmin, &path);
  clock_t end_bt = clock();
  fprintf(stderr, "Time to backtrack: " num_fmt "\n",
          (num_t)(end_bt - begin_bt) / CLOCKS_PER_SEC);

  fprintf(stderr, "%s\n", "Writing files...");
  clock_t begin_write = clock();
  for (int64_t _i = (path.ii - 2); _i > -1; _i--) {
    fprintf(index_out, "%ld %ld\n", path.index1[_i] - 1, path.index2[_i]);
  }
  fclose(index_out);
  for (int64_t _i = (path.jj - 2); _i > -1; _i--) {
    fprintf(indexs_out, "%ld %ld\n", path.index1s[_i] - 1, path.index2s[_i]);
  }
  fclose(indexs_out);
  clock_t end_write = clock();
  fprintf(stderr, "Time to write: " num_fmt "\n",
          (num_t)(end_write - begin_write) / CLOCKS_PER_SEC);

  // Free spectra matrices
  free_fvec(&spec_a);
  free_fvec(&spec_b);
  // Free distance, cost, pattern, step matrices
  free(lm);
  free(cm);
  free(sm);
  free(last_row);
  // Free the remapped paths
  free_bt_index(&path);
  // Free the step pattern data
  free_step_pattern(&p);

  clock_t program_end = clock();
  fprintf(stderr, "Total runtime: " num_fmt "\n",
          (num_t)(program_end - program_start) / CLOCKS_PER_SEC);
  return 0;
}