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

# system configuration
impurity_config = Dict{Float64,Int64}( S => 1 )
shell_config = impurity_config

# hamiltonian parameters
u = 1e-3
ϵ = -0.5*abs(u)
Γ = u/(π*1.013)
Γ = u/(π*12.66)
onsite = Dict{Float64,Vector{ComplexF64}}(
    S => [ϵ]
)
interaction = Dict{Tuple{String,Float64},Matrix{ComplexF64}}(
    ("A",0.0) => [u;;]
)
tunneling = Dict{Float64,Matrix{ComplexF64}}(
    S => [sqrt(2Γ/π);;]
)

# numerical parameters
cutoff = 500
L = 2.5
iterations = 100

# choose what to calculate among:
# - "2-particle multiplets"
# - "impurity spectrum"
# - "impurity-shell spectrum"
# - "thermodynamics"
# - "spectral"
run = "spectral"

if run=="2-particle multiplets"

    nrgfull( 
        symmetry,
        label,
        L,
        iterations,
        cutoff,
        shell_config,
        tunneling;
        impurity_config=impurity_config,
        until="2-particle multiplets"
    )

elseif run=="impurity spectrum"

    nrgfull( 
        symmetry,
        label,
        L,
        iterations,
        cutoff,
        shell_config,
        tunneling;
        impurity_config=impurity_config,
        onsite=onsite,
        interaction=interaction,
        until="impurity spectrum"
    )

elseif run=="impurity-shell spectrum"

    nrgfull( 
        symmetry,
        label,
        L,
        iterations,
        cutoff,
        shell_config,
        tunneling;
        impurity_config=impurity_config,
        onsite=onsite,
        interaction=interaction,
        until="impurity-shell spectrum"
    )

elseif run=="thermodynamics"

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
            impurity_config=impurity_config,
            onsite=onsite,
            interaction=interaction
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
        onsite=onsite,
        interaction=interaction,
        spectral=true
    )

end
