#!/bin/sh
g++ main.cc -O2 -I/home/usr/soft/include -L/home/usr/soft/lib -fopenmp -lfftw3_omp -lfftw3 -lgsl -lgslcblas -lm -Wno-format-overflow -Wno-unused-result -o Hawassi

# For MacOS user:
#  g++-14 main.cc -O2 -I/opt/homebrew/include -L/opt/homebrew/lib -fopenmp -lfftw3_omp -lfftw3 -lgsl -lgslcblas -lm -Wno-format-overflow -Wno-unused-result -o Hawassi
