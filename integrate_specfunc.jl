#!/usr/bin/env julia 

using DelimitedFiles 
include( "modules/spectral.jl" )

filename = ARGS[1] 

spec = readdlm( filename )
integral = integ( spec[:,1] , spec[:,2] ) 

@show integral


