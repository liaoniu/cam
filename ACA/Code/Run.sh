#!/bin/bash

# Give permission: chmod +x Run.sh

# Delete old files

rm ADER_LARE_1D.o
rm ADER_LARE_1D.dat

clear

# Compile Files and Run the simulation

g++ ADER_LARE_1D.cpp -o ADER_LARE_1D.o

./ADER_LARE_1D.o

# Gnuplot

gnuplot PLOT

