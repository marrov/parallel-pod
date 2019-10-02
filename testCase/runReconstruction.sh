#!/bin/bash

clear

mkdir -p output/reconstruct

cd ../cppSource/
g++ -I ./depends/eigen/ -I ./depends/ezOptionParser/ -O3 -std=c++11 -fopenmp -o reconstructField.out reconstructField.cpp
cd ../testCase/

time ../cppSource/reconstructField.out -r output/reconstruct -c output/chronos -m output/sortedModes/U -p 381600 -v 3 -nm 3 -s 200 -np 6 -nr 10