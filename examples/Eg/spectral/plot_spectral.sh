#!/bin/bash 

gnuplot --persist <<EOF 
p "spectral_Eg_zavg.dat" w l lw 1
EOF
