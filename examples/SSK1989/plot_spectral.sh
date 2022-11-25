#!/bin/bash 

FILENAME=$1

gnuplot --persist <<EOF 
p "${FILENAME}" u 1:2 w lp
EOF
