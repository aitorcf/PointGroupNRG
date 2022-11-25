#!/bin/bash 

SYSTEM="A1g"
COLUMN=2

gnuplot -persist  << EOF
set logscale x 
set yrange [-0.05:0.3]
p "thermodata/th_A.dat" u 1:2 w lp, "thermodata/th_B.dat" u 1:2 w lp
EOF








