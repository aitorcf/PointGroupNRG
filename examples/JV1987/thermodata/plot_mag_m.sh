#!/bin/bash 

gnuplot --persist <<EOF 

set logscale x 

p "th_zavg_JVm.dat" w l

EOF
