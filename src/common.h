#pragma once

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#include "common.h"

#define BUFSIZE 1024
#define SPECTRUM_MAX_LINE_LENGTH 1024
#define GROUP_FILE_MAX_LINE_LENGTH 65536
#define SFLT_EPSILON 0.000000001
#define IX(I, J) ((J) * (n) + (I))
#define IXM(I, J, K) ((J) * (K) + (I))

#define FVEC(__X__, ___Y___)         \
  fvec __X__;               \
  __X__.capacity = ___Y___; \
  __X__.size = 0;           \
  __X__.data = (num_t *)malloc(sizeof(num_t) * ___Y___);
#define UVEC(__X__, ___Y___)         \
  uvec __X__;               \
  __X__.capacity = ___Y___; \
  __X__.size = 0;           \
  __X__.data = (uint64_t *)malloc(sizeof(uint64_t) * ___Y___);
#define IVEC(__X__, ___Y___)         \
  ivec __X__;               \
  __X__.capacity = ___Y___; \
  __X__.size = 0;           \
  __X__.data = (int64_t *)malloc(sizeof(int64_t) * ___Y___);

#define num_t double
#define num_fmt "%lf"

typedef struct {
  num_t *data;
  uint64_t capacity;
  uint64_t size;
} fvec;

typedef struct {
  uint64_t *data;
  uint64_t capacity;
  uint64_t size;
} uvec;

typedef struct {
  int64_t *data;
  uint64_t capacity;
  uint64_t size;
} ivec;

typedef enum { R, C } mat_storage;

typedef struct {
  num_t *data;
  uint64_t nrow;
  uint64_t ncol;
  uint64_t capacity;
  mat_storage storage;
} mat_t;

fvec make_fvec(char const *stream);
void free_fvec(fvec *v);
void free_uvec(uvec *v);
void free_ivec(ivec *v);

void push_back_f(fvec *v, num_t *d);
void push_back_u(uvec *v, uint64_t *d);
void push_back_i(ivec *v, int64_t *d);

void print_matrix(num_t *data, int64_t n, int64_t m);
void print_imatrix(int64_t *data, int64_t n, int64_t m);
bool near(num_t A, num_t B, num_t tol);
int64_t argmin(num_t const *clist, int64_t len);
uint64_t abs_diff(uint64_t const *x, uint64_t const *y);
mat_t make_matrix_rowwise(char const *filepath, uint64_t ncol);

uint64_t lower_bound(num_t const *arr, num_t const *v, uint64_t const l,
                     num_t const *tol);