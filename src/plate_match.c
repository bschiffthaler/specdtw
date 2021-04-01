#include "plate_match.h"

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

group_t make_group(char const *filepath) {
  FILE *stream = fopen(filepath, "r");
  group_t buf;

  buf.size = 0;

  if (!stream) {
    fprintf(stderr, "Error opening input: %s\n", filepath);
    return buf;
  }

  buf.capacity = BUFSIZE;
  buf.index = (int64_t *)malloc(sizeof(int64_t) * BUFSIZE);
  buf.ppm = (num_t *)malloc(sizeof(num_t) * BUFSIZE);
  buf.value = (num_t *)malloc(sizeof(num_t) * BUFSIZE);
  buf.snr = (num_t *)malloc(sizeof(num_t) * BUFSIZE);
  buf.scale = (num_t *)malloc(sizeof(num_t) * BUFSIZE);
  buf.sample = (uint64_t *)malloc(sizeof(uint64_t) * BUFSIZE);

  uint64_t i = 0;  // running line counter
  char line[GROUP_FILE_MAX_LINE_LENGTH];

  while (!feof(stream)) {
    if (fgets(line, GROUP_FILE_MAX_LINE_LENGTH, stream) != NULL) {
      buf.sample[i] = 0;
      buf.index[i] = 0;
      sscanf(line, "%ld" num_fmt num_fmt num_fmt num_fmt "%ld", &buf.index[i],
             &buf.ppm[i], &buf.value[i], &buf.snr[i], &buf.scale[i],
             &buf.sample[i]);

      // All parse failures. Probably a header line
      if (buf.sample[i] == 0 && buf.index[i] == 0) {
        fprintf(stderr, "[WARNING]: Parse failure in '%s'. Skip line: %s",
                filepath, line);
        continue;
      }
      // Convert to 0 based
      buf.sample[i]--;
      buf.index[i]--;

      i++;
      buf.size = i;
      if (i == buf.capacity) {
        buf.index =
            (int64_t *)realloc(buf.index, sizeof(int64_t) * 2 * buf.capacity);
        buf.ppm = (num_t *)realloc(buf.ppm, sizeof(num_t) * 2 * buf.capacity);
        buf.value =
            (num_t *)realloc(buf.value, sizeof(num_t) * 2 * buf.capacity);
        buf.snr = (num_t *)realloc(buf.snr, sizeof(num_t) * 2 * buf.capacity);
        buf.scale =
            (num_t *)realloc(buf.scale, sizeof(num_t) * 2 * buf.capacity);
        buf.sample = (uint64_t *)realloc(buf.sample,
                                         sizeof(uint64_t) * 2 * buf.capacity);
        buf.capacity = 2 * buf.capacity;
      }
    }
  }
  --buf.size;
  fclose(stream);
  return buf;
}

void free_group(group_t *data) {
  free(data->index);
  free(data->ppm);
  free(data->value);
  free(data->snr);
  free(data->scale);
  free(data->sample);
}

void dtw(num_t const *sp1, num_t const *sp2, int64_t const n,
         step_pattern const *p, bt_index_t *path, num_t *dist,
         num_t *normdist) {
  // Euclidean
  num_t *lm = (num_t *)malloc(sizeof(num_t) * (n + 1) * n);  // n==m
  euclidean(sp1, sp2, lm, n, n);
  for (int64_t i = 0; i < n; i++) {
    lm[(i * (n + 1))] = 0.f;
  }
  // Cost matrix
  num_t *cm = (num_t *)malloc(sizeof(num_t) * (n + 1) * n);
  for (int64_t i = 0; i < ((n + 1) * n); i++) {
    cm[i] = NAN;
  }
  for (int64_t i = 0; i < n; i++) {
    cm[(i * (n + 1))] = 0.f;
  }
  // Step matrix
  int64_t *sm = (int64_t *)malloc(sizeof(int64_t) * (n + 1) * n);
  for (int64_t i = 0; i < ((n + 1) * n); i++) {
    sm[i] = LONG_MIN;
  }
  cost_matrix(lm, p, cm, sm, n + 1, n);
  num_t *last_row = (num_t *)malloc(sizeof(num_t) * n);
  for (int64_t i = 0; i < n; i++) {
    last_row[i] = cm[IXM(n, i, n + 1)];
  }
}

int main(int argc, char const *argv[]) {
  if (argc != 6) {
    fprintf(stderr, "usage: %s <g1> <g2> <time> <p1> <p2>\n", argv[0]);
    return 1;
  }

  group_t g1 = make_group(argv[1]);
  group_t g2 = make_group(argv[2]);
  fvec time = make_fvec(argv[3]);
  mat_t p1 = make_matrix_rowwise(argv[4], time.size);
  mat_t p2 = make_matrix_rowwise(argv[5], time.size);

  // Error reading files
  if (g1.size == 0 || g2.size == 0 || time.size == 0 || p1.nrow == 0 ||
      p2.nrow == 0) {
    return 1;
  }

  fprintf(stderr, "Have [ %ld | %ld ] peaks in [ G1 | G2 ]\n", g1.size,
          g2.size);
  fprintf(stderr, "Have %ld time points\n", time.size);
  fprintf(stderr, "Matrix dims: [[%ld, %ld], [%ld, %ld]] in [P1, P2]\n",
          p1.nrow, p1.ncol, p2.nrow, p2.ncol);

  free(p1.data);
  free(p2.data);
  free_group(&g1);
  free_group(&g2);
  free_fvec(&time);
  return 0;
}