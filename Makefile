# Makefile for Conservation Spectral SDK (Chapel)

CHPL = chpl
CHPL_FLAGS = --fast -M src
SRC = $(wildcard src/*.chpl)

.PHONY: all test bench clean

all:
	@echo "Build individual targets: make test, make bench"

test: test_chords
	./test_chords

test_chords: $(SRC) tests/test_chords.chpl
	$(CHPL) $(CHPL_FLAGS) -o $@ tests/test_chords.chpl

bench: bench_scale
	./bench_scale

bench_scale: $(SRC) benchmarks/bench_scale.chpl
	$(CHPL) $(CHPL_FLAGS) -o $@ benchmarks/bench_scale.chpl

clean:
	rm -f test_chords bench_scale *.dat

# Multi-locale build (requires GASNet-enabled Chapel)
bench_dist: $(SRC) benchmarks/bench_scale.chpl
	$(CHPL) $(CHPL_FLAGS) -nl 4 -o bench_dist benchmarks/bench_scale.chpl
