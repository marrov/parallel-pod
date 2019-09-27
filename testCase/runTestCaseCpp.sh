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
    cppbin="3DPOD.out"
fi

printf "Making directories if not present...\n\n"
mkdir -p chronos/ mode/ VTK/

if [ -e ./"$cppbin" ]
then
    printf "Running 3D POD script...\n\n"
    time ./$cppbin -i ./input/U100 -c ./chronos -m ./mode -p 381600 -v 3 -nm 5 -s 10
    printf "Done \n"
    exit 0
else
    printf "$cppbin not found. Please compile and install.\n"
    exit 1
fi

