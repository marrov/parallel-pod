#!/bin/bash

clear

mkdir -p output/sorted

cd ../cppSource/
g++ -I ./depends/eigen/ -I ./depends/ezOptionParser/ -O3 -std=c++11 -fopenmp -o sortModes.out sortModes.cpp
cd ../testCase/

time ../cppSource/sortModes.out -s output/sorted -m output/mode -r ./ -p 28 -v 3 -nm 5 -np 4