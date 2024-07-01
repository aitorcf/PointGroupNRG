#!/usr/bin/env julia

# system with J=1 total angular momentum

# load modules
package_dir = "../.."
import Pkg; Pkg.activate(package_dir)
using PointGroupNRG.NRGCalculator

# label for the system
label = "J1"

# symmetry and total angular momentum
symmetry = "J"
J = 1.0

# path to multiplet states
multiplets_dir = "multiplets"

# system configuration
impurity_config = Dict{Float64,Int64}( J => 1 )
shell_config = impurity_config

# hamiltonian parameters
ϵ = -0.1
u = 0.2
Γ = 0.002
Γ = 0.005
onsite = Dict{Float64,Vector{ComplexF64}}(
    J => [ϵ]
)
interaction = Dict{Tuple{String,Float64},Matrix{ComplexF64}}(
    ("B",1.0) => [u;;]
)
tunneling = Dict{Float64,Matrix{ComplexF64}}(
    J => [sqrt(2Γ/π);;]
)

# numerical parameters
cutoff = 100
L = 2.0
iterations = 100

# choose what to calculate
run = "spectral"

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
        impurity_config,
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
        impurity_config,
        epsilon_symparams=epsilon_symparams,
        u_symparams=u_symparams,
        until="impurity spectrum"
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
        impurity_config,
        epsilon_symparams=epsilon_symparams,
        u_symparams=u_symparams,
        until="impurity-shell spectrum"
    )

elseif run=="thermo"

    for calculation in ["CLEAN","IMP"]
        nrgfull( 
            symmetry,
            label,
            L,
            iterations,
            cutoff,
            shell_config,
            tunneling;
            calculation=calculation,
            impurity_config,
            onsite=onsite,
            interaction=interaction,
            max_SJ2=16
        )
    end

elseif run=="spectral"

    nrgfull( 
        symmetry,
        label,
        L,
        iterations,
        cutoff,
        shell_config,
        tunneling;
        impurity_config=impurity_config,
        max_SJ2=16,
        onsite=onsite,
        interaction=interaction,
        spectral=true
    )

end
