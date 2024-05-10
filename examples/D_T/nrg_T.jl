#!/usr/bin/env julia

# load modules
package_dir = "../.."
import Pkg; Pkg.activate(package_dir)
using PointGroupNRG.NRGCalculator

using Profile, PProf, ProfileVega

# input for run=="multiplets"
cg_o_dir = "/home/aitor/Bulegoa/ClebschGordan/D_T/cg_symbolic"
multiplets_dir = "multiplets"
atom_config = Dict{String,Int64}( "T" => 1 )
identityrep = "A"

# additional input for run=="spectrum"
#epsilon_symparams = Dict{String,Vector{ComplexF64}}(
#    "T" => [-0.1]
#)
#u_symparams = Dict{Tuple{String,Float64},Matrix{ComplexF64}}(
#    ("T",0.0) => [1.0;;],
#)
epsilon_symparams = Dict{String,Vector{ComplexF64}}(
    "T" => [-0.1]
)
u_symparams = Dict{Tuple{String,Float64},Matrix{ComplexF64}}(
    ("T",0.0) => [3.0;;]
)


# additional input for run=="thermo" (and run=="spectral")
label = "T"
cutoff_type = "multiplet"
cutoff_magnitude = 1000
L = 10.0
iterations = 42
shell_config = atom_config
hop_symparams = Dict{String,Matrix{ComplexF64}}(
    "T" => [0.05;;]
)

# choose what to calculate
run = "profile"

if run=="multiplets"

    multiplets_2particles_doublegroups( 
        cg_o_dir ,
        multiplets_dir ,
        atom_config ,
        identityrep
    )

elseif run=="spectrum"

    impurity_spectrum_doublegroups( 
        cg_o_dir ,
        multiplets_dir ,
        atom_config ,
        identityrep ,
        epsilon_symparams ,
        u_symparams
    )

elseif run=="imp"
    
    nrg_full_doublegroups_nonsimple( 
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

elseif run=="profile"

    nrg_full_doublegroups_nonsimple( 
        label ,
        "IMP" ,
        L ,
        iterations ,
        cutoff_type ,
        20 ,
        cg_o_dir ,
        multiplets_dir ,
        atom_config ,
        shell_config ,
        identityrep ,
        epsilon_symparams ,
        u_symparams ,
        hop_symparams
    )
    Profile.clear()
    @profview nrg_full_doublegroups_nonsimple( 
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

        nrg_full_doublegroups_nonsimple( 
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
