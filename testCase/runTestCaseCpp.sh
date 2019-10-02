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
    printf "\nRunning on Linux \n\n"
    cppbin="../cppSource/3DPOD_s.out"
fi

printf "Making directories if not present...\n\n"
mkdir -p output/chronos/ output/mode/ output/VTK/

if [ -e ./"$cppbin" ]
then
    printf "Running 3D POD script...\n\n"
    time ./$cppbin -i ./input -c output/chronos -m output/mode -p 381600 -v 1 -nm 5 -s 200 -np 6

    # Matlab sort modes
    #printf "\nRunning Matlab sorting script for modes...\n"
    #matlab -nodisplay -nosplash -nodesktop -r "run('sortModes.m');exit;"

    # Python sorted modes to VTK
    #printf "\nRunning Python mode to VTK script...\n\n"
    #python2.7 modeToVTK_U.py
    
    exit 0
else
    printf "$cppbin not found. Please compile and install.\n"
    exit 1
fi

