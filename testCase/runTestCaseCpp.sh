#!/bin/bash

printf "Making directories if not present...\n\n"
mkdir -p chronos/ mode/ VTK/

cppbin="3DPOD.out"
if [ -e ./"$cppbin" ]
then
    printf "Running 3D POD script...\n\n"
    time ./$cppbin
    printf "Done \n"
else
    printf "$cppbin not found. Please compile and install.\n"
fi
