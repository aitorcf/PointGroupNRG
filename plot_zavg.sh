#!/bin/bash 

SYSTEM=$1
COLUMN=$2

if [ $COLUMN = "6" ] # entopy -> plot exponential
then

gnuplot -persist  << EOF
set logscale x 
diffiles = system( "ls thermodata/th_diff_${SYSTEM}_*" ) 
p "thermodata/th_zavg_${SYSTEM}.dat" u 1:(exp(\$6)) w l, for [f in diffiles] f u 1:(exp(\$6)) pt 5 
EOF

else

gnuplot -persist  << EOF
set logscale x 
diffiles = system( "ls thermodata/th_diff_${SYSTEM}_*" ) 
p "thermodata/th_zavg_${SYSTEM}.dat" u 1:$COLUMN w l, for [f in diffiles] f u 1:$COLUMN pt 5 
EOF

fi







