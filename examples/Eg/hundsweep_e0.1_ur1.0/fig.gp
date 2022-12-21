#set term postscript enhanced color 16
set term pdf
set output 'hundsweep.pdf'

set logscale x 

set xlabel 'k_B T/D'
set ylabel 'k_b T {/Symbol c} / (g{/Symbol m}_B)^2'

set size 0.7,1.0

set yrange [0:0.7]

set linestyle 1 lc "forest-green" lw 2.5
set linestyle 2 lc rgb "#DC143C" lw 2.5    # pretty red
set linestyle 3 lc rgb "#4169E1" lw 2.5    # pretty blue

set key samplen 1.0
set key box
set key top left
p 'J0.0.dat' w l ls 1 title 'J/U=0',\
  'J0.1.dat' w l ls 2 title 'J/U=0.1',\
  'J0.3.dat' w l ls 3 title 'J/U=0.3'
unset label 

## Orbital 2
#set ytics 40
#set mytics 2
#set origin 0,0
#set label '(a)' at graph 0.05, graph 0.9 font ",20"
#p 'G0.01/spectral_B3ORES_3_zavg.dat' w l ls 1 title '{/Symbol G} = 0.01D',\
#  'G0.05/spectral_B3ORES_3_zavg.dat' w l ls 2 title '{/Symbol G} = 0.05D',\
#   'G0.1/spectral_B3ORES_3_zavg.dat' w l ls 3 title '{/Symbol G} = 0.1D'
#unset label
#
## Orbital 3
#set ytics 10
#set origin 0.5,0
#set label '(b)' at graph 0.05, graph 0.9 font ",20"
#p 'G0.01/spectral_B3ORES_1_zavg.dat' w l ls 1 title '{/Symbol G} = 0.01D',\
#  'G0.05/spectral_B3ORES_1_zavg.dat' w l ls 2 title '{/Symbol G} = 0.05D',\
#   'G0.1/spectral_B3ORES_1_zavg.dat' w l ls 3 title '{/Symbol G} = 0.1D'
#unset label
