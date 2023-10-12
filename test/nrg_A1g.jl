#!/usr/bin/env julia

# load package
package_dir = "../"
import Pkg; Pkg.activate(package_dir)
using PointGroupNRG.NRGCalculator

# input necessary for run=="multiplets"
cg_o_dir = "clebschgordan"
multiplets_dir = "multiplets"
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
iterations = 5
cutoff_type = "multiplet"
cutoff_magnitude = 10
shell_config = atom_config 
hop_symparams = Dict{String,Matrix{ComplexF64}}(
    "A1g" => [0.1;;]
)

# additional input for run=="spectral"
K_factor = 3.0

multiplets_2particles( 
    cg_o_dir ,
    multiplets_dir ,
    atom_config ,
    identityrep
)

impurity_spectrum( 
    cg_o_dir ,
    multiplets_dir ,
    atom_config ,
    identityrep ,
    epsilon_symparams ,
    u_symparams
)

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
    hop_symparams ;
    spectral=true ,
    K_factor=K_factor
)
