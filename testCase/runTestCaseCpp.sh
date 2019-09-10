#!/bin/bash

printf "Making directories if not present...\n\n"
mkdir -p chronos/ mode/ VTK/

printf "Running 3D POD script...\n\n"
g++ -O3 -I ./eigen 3DPOD.cpp -o 3DPOD.out; time ./3DPOD.out
printf "Done \n"