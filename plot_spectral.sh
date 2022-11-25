#!/bin/bash 

ORBITAL=$1 
Z=$2 #can be a number or "avg"

gnuplot -persist <<EOF
set xrange [-0.5:0.5]
p "spectral/spectral_${ORBITAL}_z${Z}.dat" w lp
EOF
