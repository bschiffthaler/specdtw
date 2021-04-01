LIBS=-lm
CC=gcc
NVCC=/usr/local/cuda/bin/nvcc

debug:	FLAGS += -g -O0 --compiler-options '-Wall -Wextra'
release:	FLAGS += -O2

default: debug

debug: all
release: all

src/backtrack.o: src/backtrack.c
	$(NVCC) $(FLAGS) --x cu -c -o src/backtrack.o src/backtrack.c

src/common.o: src/common.c
	$(NVCC) $(FLAGS) --x cu -c -o src/common.o src/common.c

src/cost_matrix.o: src/cost_matrix.cu
	$(NVCC) $(FLAGS) --x cu -c -o src/cost_matrix.o src/cost_matrix.cu

src/distance.o: src/distance.cu
	$(NVCC) $(FLAGS) --x cu -c -o src/distance.o src/distance.cu

src/step_pattern.o: src/step_pattern.c
	$(NVCC) $(FLAGS) --x cu -c -o src/step_pattern.o src/step_pattern.c

src/specdtw.o: src/specdtw.cu
	$(NVCC) $(FLAGS) --x cu -c -o src/specdtw.o src/specdtw.cu

src/plate_match.o: src/plate_match.c
	$(NVCC) $(FLAGS) --x cu -c -o src/plate_match.o src/plate_match.c

specdtw: src/backtrack.o src/common.o src/distance.o src/specdtw.o src/step_pattern.o src/cost_matrix.o
	mkdir -p bin
	$(NVCC) $(FLAGS) -o bin/specdtw src/backtrack.o src/common.o \
	  src/cost_matrix.o src/distance.o src/specdtw.o src/step_pattern.o $(LIBS)

plate_match: src/backtrack.o src/common.o src/distance.o src/step_pattern.o src/cost_matrix.o src/plate_match.o
	mkdir -p bin
	$(NVCC) $(FLAGS) -o bin/plate_match src/backtrack.o src/common.o \
	  src/cost_matrix.o src/distance.o src/plate_match.o src/step_pattern.o $(LIBS)

all: plate_match specdtw

clean:
	rm -f src/*.o
	rm -f bin/specdtw
	rm -f bin/plate_match