LIBS=-lm
CC=gcc
NVCC=/usr/local/cuda/bin/nvcc

debug:	CFLAGS += -g -O0 --compiler-options '-Wall -Wextra'
release:	CFLAGS += -O2

default: release

debug: all
release: all

backtrack.o: src/backtrack.c
	$(NVCC) $(CFLAGS) --x cu -c -o src/backtrack.o src/backtrack.c

common.o: src/common.c
	$(NVCC) $(CFLAGS) --x cu -c -o src/common.o src/common.c

cost_matrix.o: src/cost_matrix.cu
	$(NVCC) $(CFLAGS) --x cu -c -o src/cost_matrix.o src/cost_matrix.cu

distance.o: src/distance.cu
	$(NVCC) $(CFLAGS) --x cu -c -o src/distance.o src/distance.cu

step_pattern.o: src/step_pattern.c
	$(NVCC) $(CFLAGS) --x cu -c -o src/step_pattern.o src/step_pattern.c

main.o: src/main.cu
	$(NVCC) $(CFLAGS) --x cu -c -o src/main.o src/main.cu

all: backtrack.o common.o distance.o main.o step_pattern.o cost_matrix.o
	mkdir -p bin
	$(NVCC) $(CFLAGS) -o bin/specdtw src/backtrack.o src/common.o \
	  src/cost_matrix.o src/distance.o src/main.o src/step_pattern.o $(LIBS)

clean:
	rm -f src/*.o
	rm -f bin/specdtw