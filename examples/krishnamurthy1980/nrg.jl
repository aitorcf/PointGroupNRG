#!/usr/bin/env julia

# Phys. Rev. B 21, 1003

# load modules
package_dir = "../.."
import Pkg; Pkg.activate(package_dir)
using PointGroupNRG.NRGCalculator

# label for the system
label = "KM"

# symmetry and spin value
symmetry = "spin"
S = 0.5

# path to multiplet states
multiplets_dir = "multiplets"

# system configuration
impurity_config = Dict{Float64,Int64}( S => 1 )
shell_config = impurity_config

# hamiltonian parameters
u = 1e-3
ϵ = -0.5*abs(u)
Γ = u/(π*1.013)
Γ = u/(π*12.66)
epsilon_symparams = Dict{Float64,Vector{ComplexF64}}(
    S => [ϵ]
)
u_symparams = Dict{Tuple{String,Float64},Matrix{ComplexF64}}(
    ("A",0.0) => [u;;]
)
hop_symparams = Dict{Float64,Matrix{ComplexF64}}(
    S => [sqrt(2Γ/π);;]
)

# numerical parameters
cutoff_type = "multiplet"
cutoff_magnitude = 500
L = 2.0
iterations = 100

# choose what to calculate
run = "thermo"

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

elseif run=="spectrum"

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
        impurity_config=impurity_config,
        until="impurity spectrum",
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
            impurity_config=impurity_config,
            epsilon_symparams=epsilon_symparams,
            u_symparams=u_symparams,
        )

    end

elseif run=="spectral"

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
        impurity_config=impurity_config,
        epsilon_symparams=epsilon_symparams,
        u_symparams=u_symparams,
        spectral=true
    )

end
