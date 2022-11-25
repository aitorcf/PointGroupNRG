#!/usr/bin/env julia 

using DelimitedFiles 
using Plots

FILENAME=ARGS[1]
x_low,x_high = @.parse(Float64,ARGS[2:3])

data = readdlm( FILENAME )
energies = data[:,1]
spectral = reverse(data[:,2])

lowdiff  = map( x->abs(x-x_low)  , energies )
highdiff = map( x->abs(x-x_high) , energies )

lowidx  = findfirst( x->x==minimum(lowdiff)  , lowdiff  )
highidx = findfirst( x->x==minimum(highdiff) , highdiff )

energies = energies[lowidx:highidx]
spectral = spectral[lowidx:highidx]

maxval = maximum(spectral)
idxmax = findfirst( x->x==maxval , spectral )
halfval = maxval/sqrt(2.0)
halfdiff = map( x->abs(x-halfval) , spectral )
idxhalf = findfirst( x->x==minimum(halfdiff) , halfdiff )
halfwidth = abs( energies[idxmax] - energies[idxhalf] )

emax  = energies[idxmax]
ehalf = energies[idxhalf]

p = plot( energies , spectral )
plot!( [emax,  emax]  , [1,0] , arrow=true )
plot!( [ehalf, ehalf] , [1,0] , arrow=true )
display( p )
readline()

println( "halfwidth: $halfwidth ( max=$(energies[idxmax]) , half=$(energies[idxhalf]) )" )
