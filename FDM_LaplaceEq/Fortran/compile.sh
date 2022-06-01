#! /bin/bash

rm *.txt
gfortran -g -c -Wall -lm -O2 $1.f90 $2.f90
gfortran -g -Wall -lm -O2 $1.o $2.o -o $1.out
