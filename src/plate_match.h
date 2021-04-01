#pragma once

#include <stdio.h>
#include <stdint.h>

#include "common.h"


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

group_t make_group(FILE *stream);

void free_group(group_t *data);