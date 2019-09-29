#!/bin/bash

if [[ ${OS} =~ "Windows_NT" ]]
then
    uname_S="Windows"
else
    uname_S="$(uname -s)"
fi

if [[ "${uname_S}" =~ "Windows" ]]
then
    echo "Running on Windows"
    cppbin="3DPOD.exe"
else
    echo "Running on Linux"
    cppbin="../cppSource/P3DPOD.out"
fi

printf "Making directories if not present...\n\n"
mkdir -p output
cd output
mkdir -p chronos/ mode/ VTK/
cd ..

if [ -e ./"$cppbin" ]
then
    printf "Running parallel 3D POD script...\n\n"
    time ./$cppbin
    printf "Done \n"
    exit 0
else
    printf "$cppbin not found. Please compile and install.\n"
    exit 1
fi

