#!/bin/bash 

gnuplot --persist <<EOF 

set logscale x 

p "th_zavg_Eg.dat" u 1:6 w l

EOF
