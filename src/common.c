#include "common.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

bool near(num_t A, num_t B, num_t tol) {
  // Calculate the difference.
  num_t diff = fabs(A - B);
  A = fabs(A);
  B = fabs(B);
  // Find the largest
  num_t largest = (B > A) ? B : A;

  return (diff <= largest * tol);
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
  char line[SPECTRUM_MAX_LINE_LENGTH];
  while (!feof(stream)) {
    if (fgets(line, SPECTRUM_MAX_LINE_LENGTH, stream) != NULL) {
      double tmp = NAN;
      sscanf(line, num_fmt, &tmp);
      if (isnan(tmp)) {
        fprintf(stderr, "[WARNING]: Parse failure in '%s'. Skip line: %s",
                filepath, line);
        continue;
      }
      buf.data[i] = tmp;
      i++;
      buf.size = i;
      // Check if capacity reached
      if (i == buf.capacity) {
        buf.data = (num_t *)realloc(buf.data, sizeof(num_t) * 2 * buf.capacity);
        buf.capacity = 2 * buf.capacity;
      }
    }
  }
  fclose(stream);
  return buf;
}

void free_fvec(fvec *v) { free(v->data); }
void free_uvec(uvec *v) { free(v->data); }
void free_ivec(ivec *v) { free(v->data); }

void push_back_f(fvec *v, num_t *d) {
  if (v->size == v->capacity) {
    v->data = (num_t *)realloc(v->data, sizeof(num_t) * v->capacity * 2);
    v->capacity *= 2;
  }
  v->data[v->size] = *d;
  v->size++;
}
void push_back_u(uvec *v, uint64_t *d) {
  if (v->size == v->capacity) {
    v->data = (uint64_t *)realloc(v->data, sizeof(uint64_t) * v->capacity * 2);
    v->capacity *= 2;
  }
  v->data[v->size] = *d;
  v->size++;
}
void push_back_i(ivec *v, int64_t *d) {
  if (v->size == v->capacity) {
    v->data = (int64_t *)realloc(v->data, sizeof(int64_t) * v->capacity * 2);
    v->capacity *= 2;
  }
  v->data[v->size] = *d;
  v->size++;
}

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

uint64_t lower_bound(num_t const *arr, num_t const *v, uint64_t const l,
                     num_t const *tol) {
  uint64_t hi = l - 1;
  uint64_t lo = 0;
  uint64_t mi = hi / 2;

  while (!near(arr[mi], *v, *tol)) {
    if (arr[mi] > *v) {
      hi = mi;
    } else {
      lo = mi;
    }
    if (hi == lo) {
      mi = hi;
      break;
    }

    mi = lo + ((hi - lo) / 2);

    if (mi == lo || mi == hi) {
      break;
    }
  }
  if (!near(arr[mi], *v, *tol)) {  // Not found
    return l;
  }
  if (mi == 0) {  // Found and already at the lowest value
    return 0;
  }
  while (mi != 0 && near(arr[mi - 1], *v, *tol)) {  // Check if we can go lower
    mi--;
  }
  return mi;
}

uint64_t abs_diff(uint64_t const * x, uint64_t const * y) {
  if ((*x) > (*y)) {
    return (*x) - (*y);
  }
  return (*y) - (*x);
}