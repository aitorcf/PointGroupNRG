#!/bin/bash 

FILENAME=$1

gnuplot --persist <<EOF 
set xrange [0:0.2] 
p "${FILENAME}" u (-\$1):2 w lp
EOF
