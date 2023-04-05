#!/bin/bash 

gnuplot --persist <<EOF 

set logscale x 

p "th_cleanavg_Eg.dat" u 1:6 w l,\
  "th_impavg_Eg.dat" u 1:6 w l,\
  "th_zavg_Eg.dat" u 1:6 w l

EOF 
