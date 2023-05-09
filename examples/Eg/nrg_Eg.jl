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
include( "../../modules/zavg.jl" )
    

#   ==========   #
#%% USER INPUT %%#
#   ==========   #   

label = "Eg"

# run: multiplets, spectrum, nrg
run = "spectralzavg"
Z = generate_Z( 8 )

# parallel computation
distributed_nrg  = false
distributed_zavg = true
# number of workers (CPUs)
distworkers = length(Z)
# discretization ("standard" or "co2005")
discretization = "lanczos"

# clean system or with impurity added
calculation = "IMP"

# twisting parameter
z = 0.0

# numerical parameters
L = 3.0
betabar = 1.0

# cutoff
cutoff_type = "multiplet" 
cutoff_magnitude = 100
minmult = 0 
mine = 0.0


# one-body coupling parameters
fac = 10
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

iterations = 40

max_spin2 = 8

spectral = true
orbitalresolved = false
etafac = 1.0


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
if (distributed_zavg || distributed_nrg)

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

if run=="multiplets"

    multiplets_2part( 
                cg_o_dir ,
                asym_dir ,
                atom_config ,
                identityrep )

elseif run=="spectrum"

    atomic_spectrum( 
                label,
                z,
                cg_o_dir,
                asym_dir,
                atom_config,
                identityrep,
                epsilon_symparams,
                u_symparams)

elseif run=="nrg"

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
              orbitalresolved=orbitalresolved,
              distributed=distributed_nrg,
              z=z,
              precompute_iaj=true,
              Nz=1)
              #Nz=2)
 #               imp_spectrum=imp_spectrum)

elseif run=="thermozavg"

    if distributed_zavg
        @everywhere include( "../../modules/zavg.jl" )
        @everywhere include( "../../modules/automatization.jl" )

        @everywhere distributed_nrg = $distributed_nrg
        @everywhere label = $label
        @everywhere calculation = $calculation
        @everywhere L = $L
        @everywhere iterations = $iterations
        @everywhere cutoff_type = $cutoff_type
        @everywhere cutoff_magnitude = $cutoff_magnitude
        @everywhere cg_o_dir = $cg_o_dir
        @everywhere asym_dir = $asym_dir
        @everywhere atom_config = $atom_config
        @everywhere shell_config = $shell_config
        @everywhere identityrep = $identityrep
        @everywhere epsilon_symparams = $epsilon_symparams
        @everywhere u_symparams = $u_symparams
        @everywhere hop_symparams = $hop_symparams
        
        for calculation in ["CLEAN","IMP"]
            @sync @distributed for z in Z
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
                          distributed=distributed_nrg,
                          etafac=etafac,
                          z=z,
                          Nz=length(Z))
            end
        end
    else 

        for calculation in ["CLEAN","IMP"],
            z in Z

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
                      distributed=distributed_nrg,
                      etafac=etafac,
                      z=z,
                      Nz=length(Z))

        end

    end

    zavg_thermo( label , Z )

elseif run=="spectralzavg"

    println( "-------------------------------" )
    println( "SPECTRAL Z-AVERAGED CALCULATION" )
    println( "-------------------------------" )
    println()

    if calculation=="CLEAN"
        println( "WARNING: calculation=CLEAN => changing to
                  calculation=IMP" )
        calculation="IMP"
    end


    if distributed_zavg
        @everywhere include( "../../modules/automatization.jl" )
        @everywhere include( "../../modules/spectral.jl" )

        @everywhere distributed_nrg = $distributed_nrg
        @everywhere label = $label
        @everywhere calculation = "IMP"
        @everywhere L = $L
        @everywhere iterations = $iterations
        @everywhere cutoff_type = $cutoff_type
        @everywhere cutoff_magnitude = $cutoff_magnitude
        @everywhere cg_o_dir = $cg_o_dir
        @everywhere asym_dir = $asym_dir
        @everywhere atom_config = $atom_config
        @everywhere shell_config = $shell_config
        @everywhere identityrep = $identityrep
        @everywhere epsilon_symparams = $epsilon_symparams
        @everywhere u_symparams = $u_symparams
        @everywhere hop_symparams = $hop_symparams
        @everywhere etafac = $etafac
        @everywhere Z = $Z
        @everywhere orbitalresolved = $orbitalresolved
        
        @sync @distributed for z in Z

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
                      spectral=true,
                      etafac=etafac,
                      z=z)
        end

    else 

        Nz = length(Z)
        for (i,z) in enumerate(Z)

            println( "********************************" )
            println( "Calculation $(i)/$(Nz): z = $(z)" )
            println( "********************************" )
            println()

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
                      epsilon_symparams,
                      u_symparams,
                      hop_symparams;
                      distributed=distributed_nrg,
                      spectral=true,
                      etafac=etafac,
                      orbitalresolved=orbitalresolved,
                      z=z,
                      precompute_iaj=true,
                      Nz=1)

        end

    end

    include( "../../modules/zavg.jl" )
    zavg_spectral( label , Z )
end
