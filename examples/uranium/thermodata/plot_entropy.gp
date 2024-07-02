#!/usr/bin/env gnuplot

set term eps enhanced color
set output "uranium_entropy.eps"

set size 0.8,0.8

set logscale x

set yrange [0:1.8]

set ytics ("log(2)" log(2), "log(2)/2" log(2)/2)

p "thermo_U_diff_z0.0.dat" u 1:3 w l,\
  "thermo_U_diff_z0.0.dat" u 1:(log(2)) w l,\
  "thermo_U_diff_z0.0.dat" u 1:(log(2)/2) w l


