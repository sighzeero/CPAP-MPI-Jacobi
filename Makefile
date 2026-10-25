MPICC = mpicc
PYTHON ?= python
CFLAGS = -O3 -std=gnu99 -Wall -Wextra
CPPFLAGS =
LDLIBS = -lm
SOURCES = src/jacobi.c src/jacobi_parallel_mpi.c src/matrix.c src/timer.c src/results.c
HEADERS = src/jacobi.h src/matrix.h src/timer.h src/results.h
.PHONY: all clean test
all: src/jacobi_mpi_bin
src/jacobi_mpi_bin: $(SOURCES) $(HEADERS) Makefile
	mkdir -p results
	$(MPICC) $(CPPFLAGS) $(CFLAGS) $(SOURCES) -o $@ $(LDLIBS)
test: all
	$(PYTHON) tests/check_mpi.py
clean:
	rm -f src/jacobi_mpi_bin
