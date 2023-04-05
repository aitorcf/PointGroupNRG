#!/bin/bash 

gnuplot --persist <<EOF 

set logscale x 

p "th_zavg_A1g.dat" u 1:6 w l

EOF
