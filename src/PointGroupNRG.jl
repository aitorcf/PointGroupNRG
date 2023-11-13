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
end

# NRG calculation submodule
module NRGCalculator
export multiplets_2particles
export impurity_spectrum
export nrg_full
export generate_Z 
export generate_Zlaps
export zavg_thermo 
export zavg_spectral
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
end

using .MultipletCalculator
using .NRGCalculator

end # module PointGroupNRG
