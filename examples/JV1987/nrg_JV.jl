#!/usr/bin/env julia

#   =================   #
#%% NECESSARY MODULES %%#
#   =================   #

using DelimitedFiles
using SparseArrays

include( "../../modules/symbols.jl" )
include( "../../modules/numericals.jl" )
include( "../../modules/compoundoperators.jl" )
include( "../../modules/symmetry.jl" )
include( "../../modules/shell.jl" )
include( "../../modules/thermo.jl" )
include( "../../modules/discretization.jl" )
include( "../../modules/spectral.jl" ) 
include( "../../modules/reddiag.jl" )
include( "../../modules/diagonalization.jl" )
include( "../../modules/automatization.jl" )
    

#   ==========   #
#%% USER INPUT %%#
#   ==========   #   

label = "JV"

# stage: multiplets, spectrum, nrg
stage = "nrg"

# parallel computation
distributed  = false
# number of workers (CPUs)
distworkers = 6
# parallel method: distfor or async
method = "distfor"
# discretization ("standard" or "co2005")
discretization = "lanczos"

# clean system or with impurity added
calculation = "IMP"

# twisting parameter
z = 0.0

# numerical parameters
L = 10.0
betabar = 1.0

# cutoff
cutoff_type = "multiplet" 
cutoff_magnitude = 500
minmult = 0 
mine = 0.0


# one-body coupling parameters
eps_g   = -0.1
eps_u   = -0.1
u  = 0.25
up = 0.0
j  = 0.0
u_g0_11 = u_g0_22 = 0.5*( u + up ) + j
u_g0_12 = u_g0_21 = 0.5*( u - up )
u_u0    = u - j 
u_u1    = up - j
gam_g = 0.007
gam_u = 0.007

kR = pi

iterations = 60

max_spin2 = 8

spectral = false
etafac = 0.4


# directory where the orbital cg info is stored
cg_o_dir = "/home/aitor/Bulegoa/ClebschGordan/Ci/cg/"
# directory with information about antisymmetric combinations
asym_dir = "/home/aitor/Bulegoa/AntiSymmetricPart/Ci/";
# atomic configuration: 1 Eg irrep
atom_config  = Dict( "Ag" => 1,
                     "Au" => 1 )
shell_config = Dict( "Ag" => 1,
                     "Au" => 1 )
# identity irrep 
identityrep = "Ag"
# atomic hamiltonian parameters
epsilon_symparams = Dict( 
    "Ag" => ComplexF64[eps_g] ,
    "Au" => ComplexF64[eps_u]
)
u_symparams = Dict( 
    ("Ag",0) => ComplexF64[u_g0_11 u_g0_12 
                           u_g0_12 u_g0_22 ],
    ("Au",0) => ComplexF64[u_u0][:,:],
    ("Au",1) => ComplexF64[u_u1][:,:]
)
hop_symparams = Dict( 
    "Ag" => sqrt(2*gam_g/pi)*ComplexF64[1][:,:] ,
    "Au" => sqrt(2*gam_u/pi)*ComplexF64[1][:,:]
)

@eval g(x) = $kR*x<0.1 ? 1.0 : sin($kR*x)/($kR*x)
channel_etas = Dict( 
    "Ag" => Function[x->0.5*(1+g(1+x))],
    "Au" => Function[x->0.5*(1-g(1+x))]
)
#channel_etas = Dict( 
#    "Ag" => Function[x->0.5],
#    "Au" => Function[x->0.5]
#)

#   ==============   #
#%% EXTERNAL INPUT %%#
#   ==============   #
if length(ARGS)>0
    # clean/imp
    calculation = "IMP"
    if ARGS[1]=="CLEAN" 
        calculation = "CLEAN"
    end
    # twisting parameter
    if length(ARGS)>1
        z = parse(Float64,ARGS[2])
    end
    # spectral
    spectral = false 
    if length(ARGS)>2 
        spectral = ARGS[3]=="spectral" ? true : false 
    end
end


#   ===========   #
#%% DISTRIBUTED %%#
#   ===========   #
if distributed 

    using Distributed 

    # kill current processes
    for i in workers()
        t = rmprocs(i, waitfor=0)
        wait(t)
    end

    # add requested workers
    if distworkers ≥ nprocs()
        addprocs(distworkers)
    else 
        println( "more workers than processors!" )
    end

    println( "DISTRIBUTED CALCULATION WITH $(nworkers()) WORKERS" )

    @everywhere begin 
        using ProgressMeter
        using PartialWaveFunctions
        using StaticArrays
        include( "../../modules/symmetry.jl" )
        include( "../../modules/diagonalization.jl" )
        include( "../../modules/spectral.jl" ) 
        include( "../../modules/reddiag.jl" )
    end

else 

    println( "SERIAL CALCULATION" )

end
println()


if stage=="multiplets"

    multiplets_2part( 
                cg_o_dir ,
                asym_dir ,
                atom_config ,
                identityrep )

elseif stage=="spectrum"

    atomic_spectrum( 
                label,
                z,
                cg_o_dir,
                asym_dir,
                atom_config,
                identityrep,
                epsilon_symparams,
                u_symparams)

elseif stage=="nrg"

    nrg_full( 
                label,
                calculation,
                L,
                z,
                distributed,
                iterations,
                cutoff_type,
                cutoff_magnitude,
                max_spin2,
                cg_o_dir,
                asym_dir,
                atom_config,
                shell_config,
                identityrep,
                epsilon_symparams ,
                u_symparams,
                hop_symparams;
                discretization=discretization,
                channel_etas=channel_etas,
                spectral=spectral,
                etafac=etafac)

end
