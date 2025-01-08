all: runplot

# TODO 1: add OpenMP flags
main: main.cpp Makefile
	g++ -O3 -std=c++11 -fopenmp main.cpp -o main

run: main
	./main 1000 1

plot:
	./plot

runplot: run
	./plot

clean:
	rm -rf main out experiment hist_0.dat hist_1.dat

.PHONY: all clean run plot runplot
