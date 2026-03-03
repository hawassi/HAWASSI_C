#!/bin/sh
g++ --verbose main.c -O3 -I/home/arnidallatifah/soft/include -L/home/arnidallatifah/soft/lib/ -fopenmp -lfftw3_omp -lfftw3 -lgsl -lgslcblas -lm -Wno-unused-result -ggdb3 -g -o Hawassi2D
