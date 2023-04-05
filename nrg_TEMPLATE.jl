#!/usr/bin/env julia

#   =================   #
#%% NECESSARY MODULES %%#
#   =================   #

# |modules|: absolute path to module directory
moduledir = "/path/to/modules/"

include( "$(moduledir)/symbols.jl" )
include( "$(moduledir)/numericals.jl" )
include( "$(moduledir)/compoundoperators.jl" )
include( "$(moduledir)/symmetry.jl" )
include( "$(moduledir)/shell.jl" )
include( "$(moduledir)/thermo.jl" )
include( "$(moduledir)/discretization.jl" )
include( "$(moduledir)/spectral.jl" ) 
include( "$(moduledir)/reddiag.jl" )
include( "$(moduledir)/diagonalization.jl" )
include( "$(moduledir)/automatization.jl" )
    

#   ==========   #
#%% USER INPUT %%#
#   ==========   #   

# |label|: system name (user-given)
label = "name_of_the_system"

# [stage] (optional): multiplets, spectrum or nrg
# - multiplets: compute only multiplets. useful
#   for obtaining the two-particle multiplets needed
#   in order to write the coulomb matrix elements.
# - spectrum: compute impurity spectrum.
# - nrg (default): full nrg calculation.
stage = "multiplets"
#stage = "spectrum"
#stage = "nrg"

# [discretization] (optional):
# - "standard" (default): krishna-murthy 1980 
# - "co2005": campo 2005
# - "lanczos": matrix lanczos algorithm 
#   (necessary for energy-dependent hybridization)
discretization = "standard"


# |calculation|:
# - "CLEAN": no impurity
# - "IMP": with impurity
# this variable can be changed by script input.
# nrg_zavg.jl changes it automatically, so there is
# no need to set it manually in that case.
calculation = "IMP"

# [z] (default=0): twisting parameter in the z-averaging procedure
# (also possible to set by input)
z = 0.0

# |L|: discretization parameter Lambda
L = 10.0
# [betabar] (default=1.0): factor for thermodynamic calculations (krishna-murthy 1980)
betabar = 1.0

# cutoff
# |cutoff_type|:
# - "multiplet": keep a fixed number of multiplets 
#   + approximately degenerate multiplets
# - "energy": keep multiplets up to a given energy
#   + approximately degenerate ones
#
cutoff_type = "multiplet" 
# |cutoff_magnitude|: number of multiplets to be kept or maximum
#            energy, depending on "type"
cutoff_magnitude = 300
# [minmult] (default=0): minimum number of multiplets to be kept.
minmult = 0 

# |iterations| in the NRG calculations.
iterations = 40

# [max_spin2] (default=10): maximum spin to be found in the multiplets in the NRG iterations.
# it is used in order to compute the clebsch-gordan coefficients. 10 is a safe
# number to choose, then it can be adjusted to a lower value in order to save
# memory.
max_spin2 = 10

# [spectral] (default=false): whether to compute the spectral function.
spectral = false
# [etafac] (default=1.0): broadening factor in the spectral function calculation
etafac = 0.4

# |cg_o_dir|: directory where the orbital cg info is stored
cg_o_dir = "/path/to/Clebsch-Gordan/coefficients/"
# |asym_dir|: directory with information about multiplet states.
asym_dir = "/path/to/antisymmetric/combinations/";
# |atom_config| and |shell_config|:dictionary with structure
#       one-particle irrep => number of multiplets
# that specifies the configuration of the impurity and shell subspaces.
# - example for an Eg system:
#       atom_config  = Dict( "Eg" => 1 )
#       shell_config = Dict( "Eg" => 1 )
atom_config  = Dict()
shell_config = Dict()
# |identityrep|: name of the identity irrep for the chosen symmetry group
identityrep = "A1g"

# symmetry-adapted parameters that define the hamiltonian.
# |epsilon_symparams|: occupation energies epsilon
epsilon_symparams = Dict( 
    "Eg" => ComplexF64[eps]
)
# |u_symparams|: coulomb parameters U
u_symparams = Dict( 
    ("A1g",0) => ComplexF64[u_a1g][:,:],
    ("A2g",1) => ComplexF64[u_a2g][:,:],
    ("Eg", 0) => ComplexF64[u_eg][:,:]
)
# |hop_symparams|: hybridization parameters V
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

if stage=="multiplets"

    multiplets_2part( 
                cg_o_dir ,
                asym_dir ,
                atom_orbital_irreps ,
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
              max_spin2=max_spin2,
              z=z, 
              spectral=spectral,
              etafac=etafac,
              Nz=2)

end
