module PointGroupNRG

# user interface
export compute_multiplets
export multiplets_2particles
export impurity_spectrum
export nrg_full
export generate_Z 
export generate_Zlaps
export zavg_thermo 
export zavg_spectral

# multiplet calculation submodule
module MultipletCalculator
export compute_multiplets
include( "multiplets.jl" )
include( "multipletsnonsimple.jl" )
end
using .MultipletCalculator

# NRG calculation submodule
module NRGCalculator
export multiplets_2particles
export multiplets_2particles_doublegroups
export impurity_spectrum
export impurity_spectrum_doublegroups
export nrg_full
export nrg_full_doublegroups
export nrg_full_doublegroups_nonsimple
export generate_Z 
export generate_Zlaps
export zavg_thermo 
export zavg_spectral
export discretization_default
export tridiagonalization_default
include( "symbols.jl" )
include( "numericals.jl" )
include( "compoundoperators.jl" )
include( "symmetry.jl" )
include( "shell.jl" )
include( "band.jl" )
include( "spectral.jl" )
include( "thermo.jl" )
include( "reddiag.jl" )
include( "automatization.jl" )
include( "zavg.jl" )
include( "molecule.jl" )
include( "doublegroups.jl" )
include( "clebschgordansums.jl" )
include( "doublegroupsnonsimple.jl" )
end
using .NRGCalculator

# DMFT submodule (experimental)
module DMFT
using ..NRGCalculator
include( "dmft.jl" )
end
using .DMFT


end # module PointGroupNRG
