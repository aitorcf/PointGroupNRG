#!/usr/bin/env julia

# load modules
package_dir = "../.."
import Pkg; Pkg.activate(package_dir)
using PointGroupNRG.NRGCalculator

# input for run=="multiplets"
cg_o_dir = "../clebschgordan"
multiplets_dir = "multiplets"
atom_config = Dict{String,Int64}( "Eg" => 1 )
identityrep = "A1g"

# additional input for run=="spectrum"
epsilon_symparams = Dict{String,Vector{ComplexF64}}(
    "Eg" => [-0.0]
)
u_symparams = Dict{Tuple{String,Float64},Matrix{ComplexF64}}(
    ("Eg",0.0) => [0.2;;],
    ("A1g",0.0) => [0.5;;],
    ("A2g",1.0) => [0.05;;]
)

# additional input for run=="thermo" (and run=="spectral")
label = "Eg"
cutoff_type = "multiplet"
cutoff_magnitude = 1000
L = 10.0
iterations = 42
shell_config = atom_config
hop_symparams = Dict{String,Matrix{ComplexF64}}(
    "Eg" => [0.1;;]
)

# choose what to calculate
run = "spectral"

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
    
    nrg_full( 
        label ,
        "IMP" ,
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
        hop_symparams
    )

elseif run=="thermo"

    for calculation in ["CLEAN","IMP"]

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
