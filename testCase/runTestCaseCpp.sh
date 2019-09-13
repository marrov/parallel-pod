#!/bin/bash

printf "Making directories if not present...\n\n"
mkdir -p chronos/ mode/ VTK/

cppbin="./3DPOD.out"
if [ -e "$cppbin" ]
then
    printf "Running 3D POD script...\n\n"
    time $cppbin
    printf "Done \n"
else
    printf "Compiling and running 3D POD script...\n\n"
    g++ -std=c++11 -O3 -I ./eigen 3DPOD.cpp -o 3DPOD.out; time ./3DPOD.out
    printf "Done \n"
fi
