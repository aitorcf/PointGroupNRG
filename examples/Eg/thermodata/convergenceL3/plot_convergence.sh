#!/bin/bash 

gnuplot --persist <<EOF 

set logscale x 

p "M100.dat" w l, "M200.dat" w l, "M300.dat" w l

EOF
