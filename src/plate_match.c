#include "plate_match.h"

#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "backtrack.h"
#include "common.h"
#include "cost_matrix.h"
#include "distance.cuh"
#include "step_pattern.h"

void drop_na(char *s, uint64_t len) {
  for (uint64_t i = 0; i < (len - 1); i++) {
    if (!s[i]) {
      break;
    }
    if (s[i] == 'N' && s[i + 1] == 'A') {
      s[i] = '0';
      s[i + 1] = '.';
    }
  }
}

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
      drop_na(line, GROUP_FILE_MAX_LINE_LENGTH);
      buf.sample[i] = 0;
      buf.index[i] = 0;
      sscanf(line, "%ld" num_fmt num_fmt num_fmt num_fmt "%lu", &buf.index[i],
             &buf.ppm[i], &buf.value[i], &buf.snr[i], &buf.scale[i],
             &buf.sample[i]);

      // All parse failures. Probably a header line
      if (buf.sample[i] == 0 || buf.index[i] == 0) {
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
  norm_row(p, n, n, last_row);
  // Get distance results
  int64_t jmin = argmin(last_row, n);
  // fprintf(stderr, "%ld\n", jmin);
  *dist = cm[IXM(n, jmin, n + 1)];
  *normdist = last_row[jmin];

  // Backtrack least cost
  path->index1 = (int64_t *)malloc(sizeof(int64_t) * (n + n + 1));
  path->index1s = (int64_t *)malloc(sizeof(int64_t) * (n + n + 1));
  path->index2 = (int64_t *)malloc(sizeof(int64_t) * (n + n + 1));
  path->index2s = (int64_t *)malloc(sizeof(int64_t) * (n + n + 1));
  path->steps = (int64_t *)malloc(sizeof(int64_t) * (n + n + 1));

  backtrack(p, sm, n + 1, jmin, path);

  free(cm);
  free(sm);
  free(lm);
  free(last_row);
}

void free_group_index(group_index *data) {
  free(data->count);
  free(data->offset);
  free(data->index);
  data->size = 0;
}

group_index gindex(group_t *g, uint64_t ngroup) {
  group_index idx;
  idx.size = 0;
  idx.count = (uint64_t *)calloc(ngroup, sizeof(uint64_t));
  idx.offset = (uint64_t *)calloc(ngroup, sizeof(uint64_t));
  idx.index = (uint64_t *)malloc(sizeof(uint64_t) * g->size);

  uint64_t *tmp_storage = (uint64_t *)calloc(ngroup, sizeof(uint64_t));

  idx.size = ngroup;

  for (uint64_t i = 0; i < g->size; i++) {
    idx.count[g->sample[i]]++;
  }
  uint64_t cumsum = 0;
  for (uint64_t i = 0; i < ngroup; i++) {
    idx.offset[i] = cumsum;
    cumsum += idx.count[i];
  }
  for (uint64_t i = 0; i < g->size; i++) {
    uint64_t s = g->sample[i];
    uint64_t pos = idx.offset[s] + tmp_storage[s];
    idx.index[pos] = i;
    tmp_storage[s]++;
  }

  free(tmp_storage);

  return idx;
}

void push_back_psi(psi_vec *v, psi_peak *d) {
  if (v->size == v->capacity) {
    v->data = (psi_peak *)realloc(v->data, sizeof(psi_peak) * v->capacity * 2);
    v->capacity *= 2;
  }
  v->data[v->size] = *d;
  v->size++;
}
int comp_peak(void const *x, void const *y) {
  psi_peak *A = (psi_peak *)x;
  psi_peak *B = (psi_peak *)y;
  if (A->pos == B->pos) {
    return 0;
  }
  if (A->pos < B->pos) {
    return -1;
  }
  return 1;
}
int comp_i64(void const *x, void const *y) {
  int64_t A = *(int64_t *)x;
  int64_t B = *(int64_t *)y;
  if (A == B) {
    return 0;
  }
  if (A < B) {
    return -1;
  }
  return 1;
}

void update_diff_mat(int64_t const *rowindex, int64_t const *colindex,
                     int64_t *c_peak, uint64_t *c_diff, int64_t *b_peak,
                     uint64_t *b_diff, int64_t *n_peak, num_t const *normdist,
                     int64_t const *nrow, num_t *dist_matrix,
                     num_t const *maxii) {
  int64_t row = rowindex[*c_peak];
  int64_t col = colindex[*b_peak];
// Update the peak matrix with the difference in peak positions weighted by
// the spectral difference from DTW. A higher spectral difference incurs
// a higher penslty
#pragma omp critical(mat_update)
  {
    // Reward has two factors.
    //   normdist: Similarity of spectra according to DTW
    //   b_diff  : Distance between peaks on pseudotime
    num_t val = *maxii - ((num_t)*b_diff * (1.0 - *normdist));
    // fprintf(stderr, "U: [%ld, %ld]\n", row, col);
    if (isnan(val)) {
      fprintf(stderr,
              "[WARNING]: NAN introduced in dist matrix here: [%ld, %ld]\n",
              row, col);
    }
    int64_t index = IXM(row, col, *nrow);
    /*if (isnan(dist_matrix[index])) {
      dist_matrix[index] = 0.0;
    }*/
    dist_matrix[index] += val;
  }
  *c_peak = *n_peak;
  *b_peak = 0;
  *b_diff = ULONG_MAX;
  *c_diff = ULONG_MAX;
}

void pseudomatch(num_t const *sp1, num_t const *sp2, group_t const *g1,
                 group_t const *g2, fvec const *time, step_pattern const *p,
                 uint64_t const xi, uint64_t const xj, group_index const *g1i,
                 group_index const *g2i, int64_t const *rowindex,
                 int64_t const *colindex, num_t *dist_matrix,
                 int64_t const nrow, uint64_t *penalty) {
  int64_t n = time->size;
  num_t dist = 0;
  num_t normdist = 0;
  bt_index_t path;

  // Calculate the dynamic time warping
  dtw(sp1, sp2, n, p, &path, &dist, &normdist);

  // No warp function found
  if (isnan(normdist)) {
    return;
  }

  num_t maxii = (num_t)(path.ii - 2);

  #pragma omp critical(incr_penalty)
  { // Write to shared mem. Penalty is multiplied by 2 since we have both a 
    // forward and reverse algorithm. Max reward is therefore (2 * L - 0)
    *penalty += 2 * (path.ii - 2);
  }

  // Now let's convert the DTW indices to pseudo-time
  num_t *psiA = (num_t *)malloc(sizeof(num_t) * (path.ii - 1));
  num_t *psiB = (num_t *)malloc(sizeof(num_t) * (path.ii - 1));
  uint64_t _j = 0;
  for (int64_t i = (path.ii - 2); i > -1; i--) {
    psiA[_j] = time->data[path.index1[i] - 1];
    psiB[_j] = time->data[path.index2[i]];
    _j++;
  }
  // END convert DTW to pseudo-time

  // Match pseudotime to known peaks
  //  Tolerance for finding an index on the time axis is half the length of
  //  a time step. This allows to match more fuzzily imputed peak ppms.
  num_t tol = (time->data[1] - time->data[0]) / 2.0;

  // Vectors for peak ppm position and peak identity
  PSIVEC(pi_a, BUFSIZE);
  PSIVEC(pi_b, BUFSIZE);

  // Now let's find the indices
  for (uint64_t i = 0; i < g1i->count[xi]; i++) {
    uint64_t pos = g1i->index[g1i->offset[xi] + i];
    num_t *v = &g1->ppm[pos];
    uint64_t bound = lower_bound(psiA, v, path.ii - 1, &tol);
    int64_t idx = g1->index[pos];
    while (bound < ((uint64_t)path.ii - 1) && near(psiA[bound], *v, tol)) {
      // fprintf(stderr, "%lf, %ld, %ld, %lf\n", g1->ppm[pos], idx, bound,
      //         psiA[bound]);
      psi_peak peak;
      peak.pos = bound;
      peak.idx = idx;
      push_back_psi(&pi_a, &peak);
      bound++;
    }
  }
  for (uint64_t i = 0; i < g2i->count[xj]; i++) {
    uint64_t pos = g2i->index[g2i->offset[xj] + i];
    num_t *v = &g2->ppm[pos];
    uint64_t bound = lower_bound(psiB, v, path.ii - 1, &tol);
    int64_t idx = g2->index[pos];
    while (bound < ((uint64_t)path.ii - 1) && near(psiB[bound], *v, tol)) {
      // fprintf(stderr, "%lf, %ld, %ld, %lf\n", g2->ppm[pos], idx, bound,
      //         psiB[bound]);
      psi_peak peak;
      peak.pos = bound;
      peak.idx = idx;
      push_back_psi(&pi_b, &peak);
      bound++;
    }
  }

  // Need sorted input for next step
  qsort(pi_a.data, pi_a.size, sizeof(psi_peak), &comp_peak);
  qsort(pi_b.data, pi_b.size, sizeof(psi_peak), &comp_peak);
  // END Match pseudotime to known peaks

  // Now we calculate a best match for each peak index
  //   Get first position and identity for both spectra
  int64_t c_peak = pi_a.data[0].idx;
  uint64_t c_diff = abs_diff(&pi_a.data[0].pos, &pi_b.data[0].pos);

  int64_t b_peak = 0;
  uint64_t b_diff = ULONG_MAX;

  uint64_t ii = 0;
  uint64_t jj = 0;

  int64_t n_peak;
  uint64_t n_diff;

  while (ii < pi_a.size) {
    n_peak = pi_a.data[ii].idx;
    n_diff = abs_diff(&pi_a.data[ii].pos, &pi_b.data[jj].pos);

    if (n_peak != c_peak) {
      update_diff_mat(rowindex, colindex, &c_peak, &c_diff, &b_peak, &b_diff,
                      &n_peak, &normdist, &nrow, dist_matrix, &maxii);
    }

    // Forward algorithm: search future timepoints for better matching peaks.
    // Stop if going further into future doesn't yield same or better distance
    while (jj < pi_b.size && n_diff <= c_diff) {
      if (n_diff < b_diff) {
        b_peak = pi_b.data[jj].idx;
        b_diff = n_diff;
        c_diff = n_diff;
      }
      jj++;
      n_diff = abs_diff(&pi_a.data[ii].pos, &pi_b.data[jj].pos);
    }
    // Search backward
    while (n_diff <= c_diff) {
      if (n_diff < b_diff) {
        b_peak = pi_b.data[jj].idx;
        b_diff = n_diff;
        c_diff = n_diff;
      }
      jj--;
      n_diff = abs_diff(&pi_a.data[ii].pos, &pi_b.data[jj].pos);
      if (jj == 0) {            // Stop if we hit 0
        if (n_diff < b_diff) {  // Stil check if jj=0 is better than jj=1
          b_peak = pi_b.data[jj].idx;
          b_diff = n_diff;
        }
        break;
      }
    }
    ii++;
  }
  update_diff_mat(rowindex, colindex, &c_peak, &c_diff, &b_peak, &b_diff,
                  &n_peak, &normdist, &nrow, dist_matrix, &maxii);

  // Now the reciprocal
  c_peak = pi_b.data[0].idx;
  c_diff = abs_diff(&pi_b.data[0].pos, &pi_a.data[0].pos);

  b_peak = 0;
  b_diff = ULONG_MAX;

  ii = 0;
  jj = 0;

  while (ii < pi_b.size) {
    n_peak = pi_b.data[ii].idx;
    n_diff = abs_diff(&pi_b.data[ii].pos, &pi_a.data[jj].pos);

    if (n_peak != c_peak) {
      update_diff_mat(rowindex, colindex, &c_peak, &c_diff, &b_peak, &b_diff,
                      &n_peak, &normdist, &nrow, dist_matrix, &maxii);
    }

    // Forward algorithm: search future timepoints for better matching peaks.
    // Stop if going further into future doesn't yield same or better distance
    while (jj < pi_a.size && n_diff <= c_diff) {
      if (n_diff < b_diff) {
        b_peak = pi_a.data[jj].idx;
        b_diff = n_diff;
        c_diff = n_diff;
      }
      jj++;
      n_diff = abs_diff(&pi_b.data[ii].pos, &pi_a.data[jj].pos);
    }
    // Search backward
    while (n_diff <= c_diff) {
      if (n_diff < b_diff) {
        b_peak = pi_a.data[jj].idx;
        b_diff = n_diff;
        c_diff = n_diff;
      }
      jj--;
      n_diff = abs_diff(&pi_b.data[ii].pos, &pi_a.data[jj].pos);
      if (jj == 0) {            // Stop if we hit 0
        if (n_diff < b_diff) {  // Stil check if jj=0 is better than jj=1
          b_peak = pi_a.data[jj].idx;
          b_diff = n_diff;
        }
        break;
      }
    }
    ii++;
  }
  update_diff_mat(rowindex, colindex, &c_peak, &c_diff, &b_peak, &b_diff,
                  &n_peak, &normdist, &nrow, dist_matrix, &maxii);

  // END Now we calculate a best match for each peak index

  free(psiA);
  free(psiB);
  free(pi_a.data);
  free(pi_b.data);
  free_bt_index(&path);
}

int main(int argc, char const *argv[]) {
  if (argc != 9) {
    fprintf(stderr,
            "usage: %s <g1> <g2> <time> <p1> <p2> <dist_out> <rows_out> "
            "<cols_out>\n",
            argv[0]);
    return 1;
  }

  group_t g1 = make_group(argv[1]);
  group_t g2 = make_group(argv[2]);
  fvec time = make_fvec(argv[3]);
  mat_t p1 = make_matrix_rowwise(argv[4], time.size);
  mat_t p2 = make_matrix_rowwise(argv[5], time.size);

  FILE *dist_out = fopen(argv[6], "w");
  FILE *rows_out = fopen(argv[7], "w");
  FILE *cols_out = fopen(argv[8], "w");

  if (!dist_out || !rows_out || !cols_out) {
    fprintf(stderr, "Error opening output file\n");
    return 1;
  }

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

  step_pattern p = rjIIc();

  group_index g1i = gindex(&g1, p1.nrow);
  group_index g2i = gindex(&g2, p2.nrow);

  IVEC(all_peaks_a, BUFSIZE);
  IVEC(all_peaks_b, BUFSIZE);

  for (uint64_t i = 0; i < g1.size; i++) {
    push_back_i(&all_peaks_a, &g1.index[i]);
  }
  for (uint64_t i = 0; i < g2.size; i++) {
    push_back_i(&all_peaks_b, &g2.index[i]);
  }
  qsort(all_peaks_a.data, all_peaks_a.size, sizeof(int64_t), &comp_i64);
  qsort(all_peaks_b.data, all_peaks_b.size, sizeof(int64_t), &comp_i64);

  // Get unique peak IDs
  IVEC(rowids, BUFSIZE);
  IVEC(colids, BUFSIZE);

  push_back_i(&rowids, &all_peaks_a.data[0]);
  for (uint64_t i = 1; i < all_peaks_a.size; i++) {
    if (rowids.data[rowids.size - 1] != all_peaks_a.data[i]) {
      push_back_i(&rowids, &all_peaks_a.data[i]);
    }
  }
  push_back_i(&colids, &all_peaks_b.data[0]);
  for (uint64_t i = 1; i < all_peaks_b.size; i++) {
    if (colids.data[colids.size - 1] != all_peaks_b.data[i]) {
      push_back_i(&colids, &all_peaks_b.data[i]);
    }
  }
  fprintf(stderr, "Have matrix of [%lu, %lu] unique peak IDs\n", rowids.size,
          colids.size);
  // Now make a O(1) vector to map peak IDs to their respective row and column
  int64_t max_row = rowids.data[rowids.size - 1];
  int64_t max_col = colids.data[colids.size - 1];
  IVEC(rowindex, (max_row + 1));
  IVEC(colindex, (max_col + 1));

  for (int64_t i = 0; i < max_row + 1; i++) {
    rowindex.data[i] = -1;
  }
  for (int64_t i = 0; i < max_col + 1; i++) {
    colindex.data[i] = -1;
  }

  rowindex.size = max_row + 1;
  colindex.size = max_col + 1;
  rowindex.capacity = max_row + 1;
  colindex.capacity = max_col + 1;

  for (int64_t i = 0; i < (int64_t)rowids.size; i++) {
    rowindex.data[rowids.data[i]] = i;
  }
  for (int64_t i = 0; i < (int64_t)colids.size; i++) {
    colindex.data[colids.data[i]] = i;
  }

  num_t *dist_matrix =
      (num_t *)calloc(rowids.size * colids.size, sizeof(num_t));

  uint64_t total = (p1.nrow * p2.nrow);
  num_t d_total = (num_t)total;
  for (uint64_t i = 0; i < (rowids.size * colids.size); i++) {
    dist_matrix[i] = 0;
  }
  uint64_t ctr = 0;
  double begin_match = omp_get_wtime();
  uint64_t penalty = 0;
// Match only lower triangle
#pragma omp parallel for schedule(dynamic)
  for (uint64_t i = 0; i < p1.nrow; i++) {  // p1.nrow
    for (uint64_t j = 0; j < p2.nrow; j++) {
      double begin_submatch = omp_get_wtime();
      pseudomatch(&p1.data[p1.ncol * i], &p2.data[p2.ncol * j], &g1, &g2, &time,
                  &p, i, j, &g1i, &g2i, rowindex.data, colindex.data,
                  dist_matrix, (int64_t)rowids.size, &penalty);
      double end_match = omp_get_wtime();
      num_t t_over = end_match - begin_match;
      num_t t_sub = end_match - begin_submatch;
#pragma omp critical(pm_fp)
      {
        ctr++;
        double perc = ((double)ctr / (double)total) * 100;
        double mean = t_over / (double)ctr * omp_get_num_threads();
        double eta = ((100.0 / perc) * t_over) - t_over;
        fprintf(stderr,
                "Progress: %3.4lf%%, Overall: %.2lfs, Lap: %.2lfs, Mean: "
                "%.2lfs, ETA: %.2f, i: %lu, j: %lu\r",
                perc, t_over, t_sub, mean, eta, i, j);
      }
    }
  }

  // Panlize and norm by comparisons
  for (uint64_t i = 0; i < (rowids.size * colids.size); i++) {
    num_t nn = (dist_matrix[i] - penalty) / d_total;
    dist_matrix[i] = nn;
  }

  for (uint64_t i = 0; i < rowids.size; i++) {
    fprintf(rows_out, "%ld\n", rowids.data[i]);
    for (uint64_t j = 0; j < colids.size; j++) {
      if (i == 0) {
        fprintf(cols_out, "%ld\n", colids.data[j]);
      }
      num_t d = dist_matrix[IXM(i, j, rowids.size)];
      fprintf(dist_out, "%lf", d);
      if (j == (colids.size - 1)) {
        fprintf(dist_out, "\n");
      } else {
        fprintf(dist_out, " ");
      }
    }
  }

  fclose(dist_out);
  fclose(rows_out);
  fclose(cols_out);

  free_step_pattern(&p);

  free(p1.data);
  free(p2.data);
  free_group(&g1);
  free_group_index(&g1i);
  free_group_index(&g2i);
  free_group(&g2);
  free_fvec(&time);
  free_ivec(&all_peaks_a);
  free_ivec(&all_peaks_b);
  free_ivec(&rowids);
  free_ivec(&colids);
  free_ivec(&rowindex);
  free_ivec(&colindex);
  free(dist_matrix);
  return 0;
}