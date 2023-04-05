#!/usr/bin/env julia 

using DelimitedFiles 
using Plots
using Printf

# input determines thermo file
label = ARGS[1] 
filename = "thermodata/th_$label.dat" 

# χ(T) is used for the calculation of Tk
tm = readdlm(filename)[:,1:2]
t = tm[:,1]
logt = log10.(t)
tau = logt[2]-logt[1]
m = tm[:,2]

# we demand kTχ=0.07
diff = @.(m->abs(m-0.07))(m)
zero = logt[findfirst( x->x==minimum(diff) , diff )]

# plot
p = plot( logt[3:(end-2)] , m[3:(end-2)] , 
          xticks=-12:-2 , 
          yticks=0.0:0.05:0.25 ,
          ylims=(-0.0,0.25) )
plot!( p , 
       [zero, zero] , [0.0,0.02] , 
       arrow=true , 
       color=:black )
display(p)
readline()

toprint =  @sprintf( "%.2E", exp10(zero) ) 
println("Kondo temperature: $(toprint)D ")










