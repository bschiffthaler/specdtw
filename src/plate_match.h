#pragma once

#include <stdio.h>
#include <stdint.h>

#include "common.h"

#define PSIVEC(__X__, ___Y___)         \
  psi_vec __X__;               \
  __X__.capacity = ___Y___; \
  __X__.size = 0;           \
  __X__.data = (psi_peak *)malloc(sizeof(psi_peak) * ___Y___);

typedef struct {
  int64_t *index;
  num_t *ppm;
  num_t *value;
  num_t *snr;
  num_t *scale;
  uint64_t *sample;
  uint64_t capacity;
  uint64_t size;
} group_t;

typedef struct {
  uint64_t *count;
  uint64_t *offset;
  uint64_t *index;
  uint64_t size;
} group_index;

typedef struct
{
  uint64_t pos;
  int64_t idx;
} psi_peak;

typedef struct
{
  psi_peak * data;
  uint64_t size;
  uint64_t capacity;
} psi_vec;

void push_back_psi(psi_vec *v, psi_peak *d);
int comp_peak(void const *x, void const *y);
int comp_i64(void const *x, void const *y);

group_t make_group(char const *filepath);
group_index gindex(group_t *gm, uint64_t ngroup);

void free_group(group_t *data);
void free_group_index(group_index *data);