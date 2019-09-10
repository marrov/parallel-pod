#!/bin/bash

printf "Making directories if not present...\n\n"
mkdir -p chronos/ mode/ VTK/

printf "Running 3D POD script...\n\n"
python2.7 3DPOD_U.py
printf "Done \n"