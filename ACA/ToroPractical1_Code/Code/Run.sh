#!/bin/bash

# Give permission: chmod +x Run.sh

# Delete old files

rm EulerFV_1d.o
rm EulerFV_1dOutput.dat

clear

# Compile Files and Run the simulation

g++ EulerFV_1d.cpp -o EulerFV_1d.o

./EulerFV_1d.o

# Gnuplot

gnuplot PLOT

