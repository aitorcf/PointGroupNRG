#!/usr/bin/env julia

# Equivalent to Eg in https://doi.org/10.1016/j.cpc.2023.109032

# load modules
package_dir = "../.."
import Pkg; Pkg.activate(package_dir)
using PointGroupNRG.NRGCalculator

# label for the system
label = "E"

# symmetry
symmetry = "PS"
cg_o_dir = "../clebschgordan/O"
identityrep = "A1"

# path to multiplet states
multiplets_dir = "multiplets"

# system configuration
impurity_config = Dict{String,Int64}( "E" => 1 )
shell_config = impurity_config

# standard anderson hamiltonian
ϵ = -0.1
u = 0.3
j = 0.2
epsilon_symparams = Dict{String,Vector{ComplexF64}}(
    "E" => [-0.1]
)
u_symparams = Dict{Tuple{String,Float64},Matrix{ComplexF64}}(
    ("E",0.0) =>  [u+j;;],
    ("A1",0.0) => [u-j;;],
    ("A2",1.0) => [u-3j;;]
)

# ionic anderson model
G_ground  = (2,"A2",1.0)
G_excited = (1,"E",0.5) 
G_hop     = (1,"E",0.5)
spectrum = Dict{Tuple{Int64,String,Float64},Vector{Float64}}(
    G_ground=>[0.0] ,
    G_excited =>[0.1]
)
lehmann_iaj = Dict{NTuple{3,Tuple{Int64,String,Float64}},Array{ComplexF64,4}}(
    (G_ground,G_hop,G_excited)=>ones(ComplexF64,1,1,1,1)
)

# hybridization
Γ = 0.01 # for thermo
hop_symparams = Dict{String,Matrix{ComplexF64}}(
    "E" => [sqrt(2Γ/π);;]
)

# numerical parameters
cutoff_type = "multiplet"
cutoff_magnitude = 500
L = 5.0
iterations = 100

# choose what to calculate
run = "impurityprojections"

if run=="multiplets"

    nrg_full_allsymmetries(
        symmetry,
        label,
        L,
        iterations,
        cutoff_type,
        cutoff_magnitude,
        multiplets_dir,
        shell_config,
        hop_symparams;
        cg_o_dir=cg_o_dir,
        identityrep=identityrep,
        impurity_config=impurity_config,
        until="2-particle multiplets"
    )

elseif run=="impurity spectrum"

    nrg_full_allsymmetries( 
        symmetry,
        label,
        L,
        iterations,
        cutoff_type,
        cutoff_magnitude,
        multiplets_dir,
        shell_config,
        hop_symparams;
        cg_o_dir=cg_o_dir,
        identityrep=identityrep,
        impurity_config=impurity_config,
        until="impurity spectrum",
        epsilon_symparams=epsilon_symparams,
        u_symparams=u_symparams
    )

elseif run=="impurity-shell spectrum"

    nrg_full_allsymmetries( 
        symmetry,
        label,
        L,
        iterations,
        cutoff_type,
        cutoff_magnitude,
        multiplets_dir,
        shell_config,
        hop_symparams;
        cg_o_dir=cg_o_dir,
        identityrep=identityrep,
        impurity_config=impurity_config,
        until="impurity-shell spectrum",
        epsilon_symparams=epsilon_symparams,
        u_symparams=u_symparams
    )

elseif run=="thermo"

    for calculation in ["CLEAN","IMP"]

        nrg_full_allsymmetries( 
            symmetry,
            label,
            L,
            iterations,
            cutoff_type,
            cutoff_magnitude,
            multiplets_dir,
            shell_config,
            hop_symparams;
            calculation=calculation,
            cg_o_dir=cg_o_dir,
            identityrep=identityrep,
            impurity_config=impurity_config,
            epsilon_symparams=epsilon_symparams,
            u_symparams=u_symparams
        )

    end

elseif run=="thermoionic"

    for calculation in ["CLEAN","IMP"]

        nrg_full_allsymmetries( 
            symmetry,
            label,
            L,
            iterations,
            cutoff_type,
            cutoff_magnitude,
            multiplets_dir,
            shell_config,
            hop_symparams;
            calculation=calculation,
            cg_o_dir=cg_o_dir,
            identityrep=identityrep,
            spectrum=spectrum,
            lehmann_iaj=lehmann_iaj
        )

    end

elseif run=="thermoionic"

    for calculation in ["CLEAN","IMP"]

        nrg_full_allsymmetries( 
            "PS",
            label,
            calculation,
            L,
            iterations,
            cutoff_type,
            cutoff_magnitude,
            multiplets_dir,
            shell_config,
            hop_symparams;
            cg_o_dir=cg_o_dir,
            identityrep=identityrep,
            spectrum=spectrum,
            lehmann_iaj=lehmann_iaj
        )

    end

elseif run=="impurityprojections"

    nrg_full_allsymmetries( 
        symmetry,
        label,
        "IMP",
        L,
        iterations,
        cutoff_type,
        cutoff_magnitude,
        multiplets_dir,
        shell_config,
        hop_symparams;
        cg_o_dir=cg_o_dir,
        identityrep=identityrep,
        impurity_config=impurity_config,
        epsilon_symparams=epsilon_symparams,
        u_symparams=u_symparams,
        compute_impurity_projections=true
    )

end
