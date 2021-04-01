#include "common.h"

#include <stdint.h>
#include <stdlib.h>
#include <math.h>

bool near(num_t A, num_t B) {
  // Calculate the difference.
  num_t diff = fabs(A - B);
  A = fabs(A);
  B = fabs(B);
  // Find the largest
  num_t largest = (B > A) ? B : A;

  return (diff <= largest * SFLT_EPSILON);
}

int64_t argmin(num_t const *clist, int64_t len) {
  int64_t ii = -1;
  num_t vv = INFINITY;
  for (int64_t i = 0; i < len; i++) {
    if (clist[i] < vv) {
      ii = i;
      vv = clist[i];
    }
  }
  return ii;
}

fvec make_fvec(char const *filepath) {
  // Alloc initial vector
  FILE *stream = fopen(filepath, "r");
  fvec buf;

  buf.size = 0;

  if (!stream) {
    fprintf(stderr, "[ERROR]: opening input: %s\n", filepath);
    return buf;
  }

  buf.data = (num_t *)malloc(sizeof(num_t) * BUFSIZE);
  buf.capacity = BUFSIZE;

  uint64_t i = 0;
  while (!feof(stream)) {
    fscanf(stream, num_fmt, &buf.data[i]);
    i++;
    buf.size = i;
    // Check if capacity reached
    if (i == buf.capacity) {
      buf.data = (num_t *)realloc(buf.data, sizeof(num_t) * 2 * buf.capacity);
      buf.capacity = 2 * buf.capacity;
    }
  }
  --buf.size;
  fclose(stream);
  return buf;
}

void free_fvec(fvec *v) { free(v->data); }

void print_matrix(num_t *data, int64_t n, int64_t m) {
  for (int64_t row = 0; row < n; row++) {
    for (int64_t col = 0; col < m; col++) {
      printf("%e", data[IX(row, col)]);

      if (col == (m - 1)) {
        printf("%s", "\n");
      } else {
        printf("%s", " ");
      }
    }
  }
}

void print_imatrix(int64_t *data, int64_t n, int64_t m) {
  for (int64_t row = 0; row < n; row++) {
    for (int64_t col = 0; col < m; col++) {
      printf("%ld", data[IX(row, col)]);

      if (col == (m - 1)) {
        printf("%s", "\n");
      } else {
        printf("%s", " ");
      }
    }
  }
}

mat_t make_matrix_rowwise(char const *filepath, uint64_t ncol) {
  mat_t buf;
  buf.storage = R;
  FILE *stream = fopen(filepath, "r");

  buf.nrow = 0;
  buf.ncol = ncol;

  if (!stream) {
    fprintf(stderr, "[ERROR]: opening input: %s\n", filepath);
    return buf;
  }

  // Initial alloc for 2 rows
  buf.data = (num_t *)malloc(sizeof(num_t) * ncol * 2);
  buf.capacity = 2;

  uint64_t col_ctr = 0;
  uint64_t row_ctr = 0;
  uint64_t i = 0;

  num_t tmp = 0;
  while (!feof(stream)) {
    if (!fscanf(stream, num_fmt, &tmp)) {
      // parse failure
      fprintf(stderr, "[WARNING]: Parse failure at index %ld for %s\n", i,
              filepath);
      continue;
    }
    buf.data[i] = tmp;
    i++;
    col_ctr++;
    if (col_ctr == ncol) {
      row_ctr++;
      col_ctr = 0;  // reset
      // Check if we need to increase size
      if (row_ctr == buf.capacity) {
        buf.data =
            (num_t *)realloc(buf.data, sizeof(num_t) * 2 * ncol * buf.capacity);
        buf.capacity = 2 * buf.capacity;
      }
    }
  }
  if (col_ctr != 1) {
    fprintf(stderr,
            "[WARNING]: Data in %s was not evenly divisble by "
            "column number %ld\n",
            filepath, ncol);
  }
  buf.nrow = row_ctr;

  fclose(stream);
  return buf;
}