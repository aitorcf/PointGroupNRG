#!/bin/bash 

ORBITAL=$1

gnuplot -persist << EOF
set logscale x
diffiles=system("ls mag_diff_$(ORBITAL)_z*.dat") 
p for [f in diffiles] f w l, "mag_avg_$(ORBITAL).dat" w l lw 5
EOF
