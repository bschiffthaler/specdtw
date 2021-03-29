#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "backtrack.h"
#include "common.h"
#include "cost_matrix.cuh"
#include "distance.cuh"
#include "step_pattern.h"

int main(int argc, char const *argv[]) {
  clock_t program_start = clock();

  FILE *spec_a_handle = fopen(argv[1], "r");
  FILE *spec_b_handle = fopen(argv[2], "r");

  fprintf(stderr, "%s\n", "Reading input spectra...");
  clock_t begin_read = clock();
  fvec spec_a = make_fvec(spec_a_handle);
  fvec spec_b = make_fvec(spec_b_handle);
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
  cost_matrix(lm, &p, cm, sm, n + 1, m);
  clock_t end_cm = clock();
  fprintf(stderr, "Time to calculate cm: " num_fmt "\n",
          (num_t)(end_cm - begin_cm) / CLOCKS_PER_SEC);

  num_t *last_row = (num_t *)malloc(sizeof(num_t) * m);
  for (int64_t i = 0; i < m; i++) {
    last_row[i] = cm[IXM(n, i, n + 1)];
  }

  // norm the last row
  switch (p.norm) {
    case norm_t::NA:
      break;
    case norm_t::NM:
      for (int64_t i = 0; i < m; i++) {
        last_row[i] /= ((num_t)n + (num_t)(i + 1));
      }
      break;
    case norm_t::N:
      for (int64_t i = 0; i < m; i++) {
        last_row[i] /= ((num_t)n);
      }
      break;
    case norm_t::M:
      for (int64_t i = 0; i < m; i++) {
        last_row[i] /= ((num_t)(i + 1));
      }
      break;
  }

  // open.end=TRUE
  int64_t jmin = argmin(last_row, m);
  num_t distance = cm[IXM(n, jmin, n + 1)];
  num_t normalized_distance = last_row[jmin];

  fprintf(stderr, "Distance:" num_fmt "\n", distance);
  fprintf(stderr, "Normalized distance:" num_fmt "\n", normalized_distance);

  bt_index_t path;
  path.index1 = (int64_t *)malloc(sizeof(int64_t) * (n + m + 1));
  path.index1s = (int64_t *)malloc(sizeof(int64_t) * (n + m + 1));
  path.index2 = (int64_t *)malloc(sizeof(int64_t) * (n + m + 1));
  path.index2s = (int64_t *)malloc(sizeof(int64_t) * (n + m + 1));
  path.steps = (int64_t *)malloc(sizeof(int64_t) * (n + m + 1));

  fprintf(stderr, "%s\n", "Backtracking...");
  clock_t begin_bt = clock();
  backtrack(&p, cm, sm, n + 1, m, jmin, &path);
  clock_t end_bt = clock();
  fprintf(stderr, "Time to backtrack: " num_fmt "\n",
          (num_t)(end_bt - begin_bt) / CLOCKS_PER_SEC);

  fprintf(stderr, "%s\n", "Writing files...");
  clock_t begin_write = clock();
  FILE *outf = fopen("index1_2.txt", "w");
  for (int64_t _i = (path.ii - 2); _i > -1; _i--) {
    fprintf(outf, "%ld %ld\n", path.index1[_i] - 1, path.index2[_i]);
  }
  fclose(outf);
  outf = fopen("index1_2s.txt", "w");
  for (int64_t _i = (path.jj - 2); _i > -1; _i--) {
    fprintf(outf, "%ld %ld\n", path.index1s[_i] - 1, path.index2s[_i]);
  }
  fclose(outf);
  clock_t end_write = clock();
  fprintf(stderr, "Time to write: " num_fmt "\n",
          (num_t)(end_write - begin_write) / CLOCKS_PER_SEC);

  // Close files
  fclose(spec_a_handle);
  fclose(spec_b_handle);
  // Free spectra matrices
  free_fvec(&spec_a);
  free_fvec(&spec_b);
  // Free distance, cost, pattern, step matrices
  free(lm);
  free(cm);
  free(sm);
  free(last_row);
  // Free the remapped paths
  free(path.index1);
  free(path.index1s);
  free(path.index2);
  free(path.index2s);
  free(path.steps);
  // Free the step pattern data
  free_step_pattern(&p);

  clock_t program_end = clock();
  fprintf(stderr, "Total runtime: " num_fmt "\n",
          (num_t)(program_end - program_start) / CLOCKS_PER_SEC);
  return 0;
}