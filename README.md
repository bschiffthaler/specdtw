# Spec-DTW - Efficient dynamic time warping of spectra

![spectrum](https://raw.githubusercontent.com/bschiffthaler/specdtw/main/doc/img/align.png)

This is an implementation of a single subset Toni Giorgino's `dtw` algorithm to 
push throughput a bit higher and to not have R as a dependency.

Specifically, we implement this call:

```r
dtw(
  x,
  y,
  step.pattern=rabinerJuangStepPattern(2, "c"),
  open.begin=TRUE,
  open.end=TRUE,
  window.type="none"
)
```

The distance matrix is done in CUDA, albeit there's no huge gain here. The rest is plain C. `specdtw` is ~4x faster compared to running through R.


## Compiling

Edit `Makefile` to point to your `nvcc` location.

```bash
make release
```

## Usage

```bash
specdtw <query> <reference> <index_out> <indexs_out> <dist_out>
```

You can try by runing the example:

```bash
./bin/specdtw \
  example/spec1.tsv \
  example/spec2.tsv \
  example/index1_2.txt \
  example/index1_2s.txt \
  example/dist.txt
```

## Output

Currently hardcoded in `$PWD`:

* `index1_2.txt` - A 2-column matrix with the **0-based** indices analogous to R's `dtw$index1` and `dtw$index2`
* `index1_2s.txt` - A 2-column matrix with the **0-based** indices analogous to R's `dtw$index1s` and `dtw$index2s`
