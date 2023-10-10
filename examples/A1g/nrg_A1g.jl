#!/usr/bin/env julia

# load modules
package_dir = "/home/aitor/Bulegoa/PointGroupNRG"
#cd(package_dir)
import Pkg; Pkg.activate(package_dir)
using PointGroupNRG

# input necessary for run=="multiplets"
cg_o_dir = "/home/aitor/Bulegoa/PointGroupNRG/examples/clebschgordan"
multiplets_dir = "/home/aitor/Bulegoa/PointGroupNRG/examples/A1g/multiplets"
atom_config = Dict{String,Int64}( "A1g" => 1 )
identityrep = "A1g"

# additional input for run=="spectrum"
epsilon_symparams = Dict{String,Vector{ComplexF64}}(
    "A1g" => [-0.1]
)
u_symparams = Dict{Tuple{String,Float64},Matrix{ComplexF64}}(
    ("A1g",0.0) => [0.4;;]
)

# additional input for run=="thermo"
label = "A1g"
L = 3.0
iterations = 42
cutoff_type = "multiplet"
cutoff_magnitude = 100
shell_config = atom_config 
hop_symparams = Dict{String,Matrix{ComplexF64}}(
    "A1g" => [0.1;;]
)

# additional input for run=="spectral"
K_factor = 3.0

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
        spectral=true ,
        K_factor=K_factor
    )

end
