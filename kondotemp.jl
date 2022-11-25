#!/usr/bin/env julia 

using DelimitedFiles 
using Plots
using Printf

# input determines thermo file
label = ARGS[1] 
filename = "thermodata/th_zavg_$label.dat" 

# χ(T) is used for the calculation of Tk
tm = readdlm(filename)[:,1:2]
t = tm[:,1]
logt = log10.(t)
tau = logt[2]-logt[1]
m = tm[:,2]

# we want the inflexion point, so we need the second derivative
m2 = [(m[i+2]-2*m[i]+m[i-2])/(4*tau^2) for i=3:(length(m)-2)]

# zeros of the second derivative
sign(x) = x>0 ? 1 : -1
ii0 = [i for i=1:(length(m2)-1) if sign(m2[i])==-sign(m2[i+1])]
zeros = [(logt[i+2]+logt[i+3])/2 for i in ii0]

# isolate relevant zeros 
idx_maxmag = findfirst( x->x==maximum(m) , m )
idx_minmag = findfirst( x->x==minimum(m) , m )
zeros = [ zero for zero in zeros 
               if (zero>logt[idx_minmag] && zero<logt[idx_maxmag]) ]

# plot
p = plot( logt[3:(end-2)] , m[3:(end-2)] , 
          xticks=-12:-2 , 
          yticks=0.0:0.05:0.25 ,
          ylims=(-0.0,0.25) )
for zero in zeros 
    plot!( p , 
           [zero, zero] , [0.0,0.02] , 
           arrow=true , 
           color=:black )
end
display(p)
readline()

println( @sprintf( "Kondo temperature = %.2E D", exp10(zeros[1]) ) )
