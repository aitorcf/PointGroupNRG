#!/usr/bin/env julia

using DelimitedFiles
using Glob
include("modules/shell.jl")


# input processing
orbital = ARGS[1]
Nz      = parse(Int64,ARGS[2])
calculation = "IMP"

# script and z values
scriptname = "./nrgdynamics_$(orbital).jl"
Z = get_Z(Nz) #.+ (0.5/Nz)

# delete previous files
files = glob( "spectral/spectral_$(orbital)*.dat" )
rm.(files)

println( "& ====================== &" )
println( "& Z-AVERAGED CALCULATION &" )
println( "& ====================== &" )
@show orbital 
@show scriptname
@show Nz 
@show Z
println()

# doing the computation
@time for z in Z  
    cmd = `./$scriptname $calculation $z` 
    run(cmd)
    cmd = `mv spectral/spectral.dat spectral/spectral_$(orbital)_z$(z).dat`
    run(cmd)
end

#%% averaging
println( "Averaging over values of z..." )
data = Dict()
for z in Z 
    filename = "spectral/spectral_$(orbital)_z$(z).dat"
    data[z] = readdlm( filename )
end
omegas = data[Z[1]][:,1]
data_zavg = [0.0 for _=1:length(omegas)]
for z in Z
    data_zavg .+= data[z][:,2]./Nz
end
zavgfile = "spectral/spectral_$(orbital)_zavg.dat"
cmd = `touch $zavgfile`
writedlm( zavgfile , [omegas data_zavg] )

