#!/usr/bin/env julia 

using DelimitedFiles

filename = ARGS[1]

data = readdlm(filename)
energies = data[:,1]
spectral = data[:,2]

integral = sum( (spectral[i]+spectral[i+1])*(energies[i+1]-energies[i])/2.0 
                for i=1:(length(energies)-1) )
@show integral
