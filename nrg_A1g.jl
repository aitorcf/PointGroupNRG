#!/usr/bin/env julia


#   =================   #
#%% NECESSARY MODULES %%#
#   =================   #

using DelimitedFiles
using SparseArrays

include( "modules/symbols.jl" )
include( "modules/numericals.jl" )
include( "modules/compoundoperators.jl" )
include( "modules/symmetry.jl" )
include( "modules/shell.jl" )
include( "modules/thermo.jl" )
include( "modules/discretization.jl" )
include( "modules/spectral.jl" ) 
include( "modules/reddiag.jl" )
include( "modules/diagonalization.jl" )
include( "modules/automatization.jl" )
include( "modules/templates.jl" )
    

#   ==========   #
#%% USER INPUT %%#
#   ==========   #   

label = "A1g"

# parallel computation
distributed  = false
# number of workers (CPUs)
distworkers = 6
# parallel method: distfor or async
method = "distfor"
# discretization ("standard" or "co2005")
discretization = "standard"

# clean system or with impurity added
calculation = "IMP"

# twisting parameter
z = 0.0

# numerical parameters
L = 10.0
betabar = 1.0

# cutoff
cutoff_type = "multiplet" 
cutoff_magnitude = 30
minmult = 0 
mine = 0.0


# one-body coupling parameters
eps = -0.1
u   = -2*eps
gam = 0.01

gam=sqrt(2*gam/pi)

iterations = 30

max_spin2 = 8

spectral = false
etafac = 1.0


# directory where the orbital cg info is stored
cg_o_dir = "/home/aitor/Bulegoa/ClebschGordan/Oh/cg_symbolic/"
# directory with information about antisymmetric combinations
asym_dir = "/home/aitor/Bulegoa/AntiSymmetricPart/Oh/";
# orbital irreps present in the atom
atom_orbital_irreps = String["A1g"]
# atomic configuration: 1 Eg irrep
atom_config  = Dict( "A1g" => 1 )
shell_config = Dict( "A1g" => 1 )
# identity irrep 
identityrep = "A1g"
# atomic hamiltonian parameters
epsilon_symparams = Dict( 
    ("A1g",1) => ComplexF64(eps)
)
u_symparams = Dict( 
    ("A1g",0) => ComplexF64[u][:,:]
)
hop_symparams = Dict( 
    ("A1g") => gam*ComplexF64[1][:,:]
)

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
        include( "modules/symmetry.jl" )
        include( "modules/diagonalization.jl" )
        include( "modules/spectral.jl" ) 
        include( "modules/reddiag.jl" )
    end

else 

    println( "SERIAL CALCULATION" )

end
println()


#multiplets_2part( 
#            cg_o_dir ,
#            asym_dir ,
#            atom_orbital_irreps ,
#            atom_config ,
#            identityrep )

nrg_full_thermo( 
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
            spectral=spectral,
            etafac=etafac,
            betabar=betabar)
