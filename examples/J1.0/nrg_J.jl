#!/usr/bin/env julia

# load modules
package_dir = "../.."
import Pkg; Pkg.activate(package_dir)
using PointGroupNRG.NRGCalculator

J = 1.0
label = "J$J"

# input for run=="multiplets"
multiplets_dir = "multiplets"
atom_config = Dict{Float64,Int64}( J => 1 )


# additional input for run=="spectrum"
epsilon_symparams = Dict{Float64,Vector{ComplexF64}}(
    J => [-0.1]
)
u_symparams = Dict{Tuple{String,Float64},Matrix{ComplexF64}}(
    ("A",1.0) => [0.5;;],
)

# additional input for run=="thermo" (and run=="spectral")
cutoff_type = "multiplet"
cutoff_magnitude = 600
L = 10.0
iterations = 42
shell_config = atom_config
hop_symparams = Dict{Float64,Matrix{ComplexF64}}(
    J => [0.07;;]
)

# choose what to calculate
run = "thermo"

if run=="multiplets"

    multiplets_2particles( 
        cg_o_dir ,
        multiplets_dir ,
        atom_config ,
        identityrep
    )

elseif run=="spectrum"

    impurity_spectrum( 
        cg_o_dir ,
        multiplets_dir ,
        atom_config ,
        identityrep ,
        epsilon_symparams ,
        u_symparams
    )

elseif run=="imp"

    nrg_full_totalangularmomentum( 
        label ,
        "IMP" ,
        L ,
        iterations ,
        cutoff_type ,
        cutoff_magnitude ,
        multiplets_dir ,
        atom_config ,
        shell_config ,
        epsilon_symparams ,
        u_symparams ,
        hop_symparams
    )

elseif run=="thermo"

    for calculation in ["CLEAN","IMP"]

        nrg_full_totalangularmomentum( 
            label ,
            calculation ,
            L ,
            iterations ,
            cutoff_type ,
            cutoff_magnitude ,
            multiplets_dir ,
            atom_config ,
            shell_config ,
            epsilon_symparams ,
            u_symparams ,
            hop_symparams
        )

    end

elseif run=="spectral"

    calculation = "IMP"

    nrg_full( 
        label ,
        calculation ,
        L ,
        iterations ,
        cutoff_type ,
        cutoff_magnitude ,
        cg_o_dir ,
        multiplets_dir ,
        atom_config ,
        shell_config ,
        identityrep ,
        epsilon_symparams ,
        u_symparams ,
        hop_symparams ;
        spectral=true
    )

end
