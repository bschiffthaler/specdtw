# Spec-DTW - Efficient dynamic time warping of spectra

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
specdtw <query> <reference>
```

You can try by runing the example:

```bash
./bin/specdtw example/spec1.tsv example/spec2.tsv
```
