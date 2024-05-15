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

    #include( "multiplets.jl" )
    #include( "multipletsnonsimple.jl" )
    include( "multiplets/classes.jl" )
    include( "multiplets/subroutines.jl" )
    include( "multiplets/clebschgordan.jl" )
    include( "multiplets/doublegroup.jl" )
    include( "multiplets/pointspin.jl" )
    include( "multiplets/totalangularmomentum.jl" )
    include( "multiplets/interface.jl" )

end
using .MultipletCalculator

# NRG calculation submodule
module NRGCalculator

    export multiplets_2particles
    #export multiplets_2particles_doublegroups
    export impurity_spectrum
    #export impurity_spectrum_doublegroups
    export nrg_full
    export nrg_full_doublegroup
    export nrg_full_doublegroup_nonsimple
    export nrg_full_totalangularmomentum
    export nrg_full_pointspin
    export nrg_full_allsymmetries
    export generate_Z 
    export generate_Zlaps
    export zavg_thermo 
    export zavg_spectral
    export discretization_default
    export tridiagonalization_default

    include( "nrg/symbols.jl" )
    include( "nrg/numericals.jl" )
    include( "nrg/compoundoperators.jl" )
    include( "nrg/symmetry.jl" )
    include( "nrg/shell.jl" )
    include( "nrg/band.jl" )
    include( "nrg/clebschgordansums.jl" )
    include( "nrg/spectral.jl" )
    include( "nrg/thermo.jl" )
    include( "nrg/reddiag.jl" )
    include( "nrg/automatization.jl" )
    include( "nrg/zavg.jl" )
    include( "nrg/molecule.jl" )
#    include( "nrg/doublegroups_old.jl" )
#    include( "nrg/doublegroupsnonsimple_old.jl" ) # old
    include( "nrg/doublegroup.jl" )
    include( "nrg/pointspin.jl" )
    include( "nrg/totalangularmomentum.jl" )
    include( "nrg/allsymmetries.jl")

end
using .NRGCalculator

# DMFT submodule (experimental)
module DMFT
using ..NRGCalculator
include( "dmft/dmft.jl" )
end
using .DMFT

end # module PointGroupNRG
