#include "common.h"

#include <stdint.h>
#include <stdlib.h>

fvec make_fvec(FILE *stream) {
  // Alloc initial vector
  fvec buf;
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

void print_matrix(int64_t *data, int64_t n, int64_t m) {
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
