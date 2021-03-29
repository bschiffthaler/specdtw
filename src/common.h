#pragma once

#include <stdio.h>
#include <stdint.h>

#include "common.h"

#define BUFSIZE 1024
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

fvec make_fvec(FILE *stream);
void free_fvec(fvec *v);
void print_matrix(num_t *data, int64_t n, int64_t m);
void print_matrix(int64_t *data, int64_t n, int64_t m);
