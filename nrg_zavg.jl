#!/usr/bin/env julia

using DelimitedFiles
using Glob
include("modules/shell.jl")


# input processing
label = ARGS[1]
Nz      = parse(Int64,ARGS[2])
calculations = ("CLEAN","IMP")
spectral = false
if length(ARGS)>2 
    if ARGS[3]=="spectral" 
        spectral=true
        calculations = ("IMP",)
    else
        calculations = ARGS[3]=="BOTH" ? ("CLEAN","IMP") : (ARGS[3],) # CLEAN, IMP OR BOTH
    end
end


# script and z values
scriptname = "./nrg_$(label).jl"
Z = get_Z(Nz) #.+ (0.5/Nz)

# delete previous files
for C in calculations 
    c = lowercase( C )
    files = glob( "thermodata/thermo_$(c)_$(label)_z*.dat" )
    append!( files , glob("thermodata/th_diff_$(label)_z*.dat") )
    spectral && append!( files , glob("spectral/spectral_$(label)_z*.dat") )
    rm.(files)
end

println( "& ====================== &" )
println( "& Z-AVERAGED CALCULATION &" )
println( "& ====================== &" )
@show label 
@show scriptname
@show Nz 
@show Z
@show calculations
println()

# doing the computation
@time begin
for calculation in calculations, z in Z  
    if spectral 
        cmd = `./$scriptname $calculation $z spectral` 
        run(cmd)
        cmd = `cp spectral/spectral.dat spectral/spectral_$(label)_z$(z).dat`
        run(cmd)
    else 
        cmd = `./$scriptname $calculation $z` 
        run(cmd)
    end
end
end

#%% thermo averaging
if (calculations!==("CLEAN",) && !spectral)

    println( "Averaging thermodynamic functions over values of z..." )
    th_tot  = Dict()
    thclean_tot = Dict()
    thimp_tot = Dict()

    for z in Z
        th_z   = readdlm( "thermodata/th_diff_$(label)_z$z.dat" )
        thclean_z = readdlm( "thermodata/thermo_clean_$(label)_z$z.dat" ) 
        thimp_z = readdlm( "thermodata/thermo_imp_$(label)_z$z.dat" ) 
        t = th_z[:,1] 
        th_tot[z] = Dict( round(t[i],sigdigits=2)=>th_z[i,:] for i in 1:length(t) )
        thclean_tot[z] = Dict( round(t[i],sigdigits=2)=>thclean_z[i,:] for i in 1:length(t) ) 
        thimp_tot[z] = Dict( round(t[i],sigdigits=2)=>thimp_z[i,:] for i in 1:length(t) ) 
    end

    th_zavg = Dict()
    thclean_zavg = Dict()
    thimp_zavg = Dict()
    T = sort(collect(keys(th_tot[Z[end]])))
    T = T[Nz:(end-Nz)]
    for t in T
        th_zavg[t]  = sum( th_tot[z][t] for z in Z )/Nz
        thclean_zavg[t] = sum( thclean_tot[z][t] for z in Z )/Nz 
        thimp_zavg[t] = sum( thimp_tot[z][t] for z in Z )/Nz 
    end

    println()
    
    th_zavg_vec = [th_zavg[t] for t in T]
    thclean_zavg_vec = [thclean_zavg[t] for t in T]
    thimp_zavg_vec = [thimp_zavg[t] for t in T]
    writedlm( "thermodata/th_zavg_$(label).dat" , th_zavg_vec )
    writedlm( "thermodata/th_cleanavg_$(label).dat" , thclean_zavg_vec )
    writedlm( "thermodata/th_impavg_$(label).dat" , thimp_zavg_vec )
end
#%% spectral averaging
if spectral 

    println( "Averaging spectral function over values of z..." )
    data = Dict()
    for z in Z 
        filename = "spectral/spectral_$(label)_z$(z).dat"
        data[z] = readdlm( filename )
    end
    omegas = data[Z[1]][:,1]
    data_zavg = zeros(Float64,length(omegas))
    for z in Z
        data_zavg .+= data[z][:,2]./Nz
    end

    zavgfile = "spectral/spectral_$(label)_zavg.dat"
    writedlm( zavgfile , [omegas data_zavg] )

end
