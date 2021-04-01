#pragma once

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#include "common.h"

#define BUFSIZE 1024
#define GROUP_FILE_MAX_LINE_LENGTH 65536
#define SFLT_EPSILON 0.000000001
#define IX(I, J) ((J) * (n) + (I))
#define IXM(I, J, K) ((J) * (K) + (I))

#define num_t double
#define num_fmt "%lf"

typedef struct {
  num_t *data;
  uint64_t capacity;
  uint64_t size;
} fvec;

typedef enum {
  R,
  C
} mat_storage;

typedef struct {
  num_t * data;
  uint64_t nrow;
  uint64_t ncol;
  uint64_t capacity;
  mat_storage storage;
} mat_t;

fvec make_fvec(char const *stream);
void free_fvec(fvec *v);
void print_matrix(num_t *data, int64_t n, int64_t m);
void print_imatrix(int64_t *data, int64_t n, int64_t m);
bool near(num_t A, num_t B);
int64_t argmin(num_t const *clist, int64_t len);
mat_t make_matrix_rowwise(char const *filepath, uint64_t ncol);
