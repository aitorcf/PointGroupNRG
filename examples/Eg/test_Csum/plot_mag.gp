#!/usr/bin/gnuplot --persist 

set logscale x 

p "old.dat" w l,\
  "new.dat" w l,\
  "newnew.dat" w l
