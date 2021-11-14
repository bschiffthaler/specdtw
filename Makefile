LIBS=-lm
CC=gcc
NVCC=nvcc
CFLAGS=-Wall -Wextra -fopenmp -pthread
CFLAGS+=$(CFLAGS_EXTRA)

debug: FLAGS += --default-stream per-thread -g -O0 --compiler-options '$(CFLAGS)' 
release: FLAGS += --default-stream per-thread -g -O3 --compiler-options '$(CFLAGS)' 

default: debug

debug: all
release: all

src/backtrack.o: src/backtrack.c
	$(NVCC) $(FLAGS) -c -o src/backtrack.o src/backtrack.c

src/common.o: src/common.c
	$(NVCC) $(FLAGS) -c -o src/common.o src/common.c

src/cost_matrix.o: src/cost_matrix.c
	$(NVCC) $(FLAGS) -c -o src/cost_matrix.o src/cost_matrix.c

src/gpu_cost_matrix.o: src/gpu_cost_matrix.cu
	$(NVCC) $(FLAGS) -c -o src/gpu_cost_matrix.o src/gpu_cost_matrix.cu

src/distance.o: src/distance.cu
	$(NVCC) $(FLAGS) -c -o src/distance.o src/distance.cu

src/step_pattern.o: src/step_pattern.c
	$(NVCC) $(FLAGS) -c -o src/step_pattern.o src/step_pattern.c

src/specdtw.o: src/specdtw.c
	$(NVCC) $(FLAGS) -c -o src/specdtw.o src/specdtw.c

src/plate_match.o: src/plate_match.c
	$(NVCC) $(FLAGS) -c -o src/plate_match.o src/plate_match.c

specdtw: src/backtrack.o src/common.o src/distance.o src/specdtw.o src/step_pattern.o src/cost_matrix.o src/gpu_cost_matrix.o
	mkdir -p bin
	$(NVCC) $(FLAGS) -o bin/specdtw src/backtrack.o src/common.o \
	  src/cost_matrix.o src/gpu_cost_matrix.o src/distance.o src/specdtw.o src/step_pattern.o $(LIBS)

plate_match: src/backtrack.o src/common.o src/distance.o src/step_pattern.o src/cost_matrix.o src/gpu_cost_matrix.o src/plate_match.o
	mkdir -p bin
	$(NVCC) $(FLAGS) -o bin/plate_match src/plate_match.o src/backtrack.o src/common.o \
	  src/cost_matrix.o src/gpu_cost_matrix.o src/distance.o src/step_pattern.o $(LIBS)

all: plate_match specdtw

clean:
	rm -f src/*.o
	rm -f bin/specdtw
	rm -f bin/plate_match
