#!/usr/bin/env julia


#   =================   #
#%% NECESSARY MODULES %%#
#   =================   #

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

label = "Eg"

# stage: multiplets, spectrum, nrg
stage = "nrg"

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
cutoff_magnitude = 1000
minmult = 0 
mine = 0.0


# one-body coupling parameters
fac = 1
eps   = -0.1*fac
u_11  = 1*abs(eps)
u_h   = 0.3*u_11
u_12  = u_11-2*u_h
gam   = u_11/pi*0.1
# coulomb coupling parameters
#u_eg  = u_11
#u_a1g = u_12 + u_h/2.0
#u_a2g = u_12 - u_h/2.0
u_eg  = u_11 - u_h/4.0
u_a1g = u_11 + u_h/4.0
u_a2g = u_11 - 3.0*u_h/4.0

iterations = 10

max_spin2 = 8

spectral = false
etafac = 0.4


# directory where the orbital cg info is stored
cg_o_dir = "/home/aitor/Bulegoa/ClebschGordan/Oh/cg_symbolic/"
# directory with information about antisymmetric combinations
asym_dir = "/home/aitor/Bulegoa/AntiSymmetricPart/Oh/";
# atomic configuration: 1 Eg irrep
atom_config  = Dict( "Eg" => 1 )
shell_config = Dict( "Eg" => 1 )
# identity irrep 
identityrep = "A1g"
# atomic hamiltonian parameters
epsilon_symparams = Dict( 
    "Eg" => ComplexF64[eps]
)
u_symparams = Dict( 
    ("A1g",0) => ComplexF64[u_a1g][:,:],
    ("A2g",1) => ComplexF64[u_a2g][:,:],
    ("Eg", 0) => ComplexF64[u_eg][:,:]
)
hop_symparams = Dict( 
    "Eg" => sqrt(2*gam/pi)*ComplexF64[1][:,:]
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
        include( "../../modules/symmetry.jl" )
        include( "../../modules/diagonalization.jl" )
        include( "../../modules/spectral.jl" ) 
        include( "../../modules/reddiag.jl" )
    end

else 

    println( "SERIAL CALCULATION" )

end
println()

# improve efficiency of dictionary lookup
max_population_atom  = 4
max_population_shell = 4
Ndim = max_population_atom + iterations*max_population_shell
Idim = 3
Sdim = max_spin2
maxISdim = maximum((Idim,Sdim))
ISdim = Idim*Sdim
Gdim = Ndim*Idim*Sdim
dimtup = ( ISdim , Sdim , 1 ) 
@eval function Base.hash( x::NTuple{3,NTuple{3,Int64}} )
    UInt64(sum( sum.((x[1].*$dimtup,x[2].*$dimtup,x[3].*$dimtup)).*($(Gdim^2),$Gdim,1) ))
end

imp_spectrum = Dict( 
        (2,"A2g",1.0) => [0.0],
        (1,"Eg",0.5)  => [1.0] )

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

    nrg_full( label,
              calculation,
              L,
              iterations,
              cutoff_type,
              cutoff_magnitude,
              cg_o_dir,
              asym_dir,
              atom_config,
              shell_config,
              identityrep,
              epsilon_symparams ,
              u_symparams,
              hop_symparams;
              spectral=spectral,
              etafac=etafac,
              z=z,
              precompute_iaj=true)
              #Nz=2)
 #               imp_spectrum=imp_spectrum)

end
