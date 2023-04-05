#!/bin/bash 

gnuplot --persist <<EOF 

set logscale x 

p "th_zavg_JV.dat" w l

EOF
