#!/usr/bin/env julia

# Infinite-U (where U refers to the Coulomb) interaction
# for an impurity with electrons having total angular
# momentum J. The value of J can be changed (here and
# in multiplets.jl), but beware of increasing it past
# the threshold dimensionality (2J+1=6, equivalent to 
# a standard 3-channel Kondo model, is too much).

# load modules
package_dir = "../.."
import Pkg; Pkg.activate(package_dir)
using PointGroupNRG.NRGCalculator

# label for the system
label = "IU" # infinite U

# symmetry and total angular momentum
symmetry = "J"
J = 1.0

# path to multiplet states
multiplets_dir = "multiplets"

# system configuration
shell_config = Dict{Float64,Int64}( J => 1 )

# ionic model
Δ = 0.1  # spectral
Δ = 1e-3 # thermodynamics
G_ground = (1,"A",J)
G_tunnel = G_ground
G_excited = (0,"A",0.0)
spectrum = Dict{Tuple{Int64,String,Float64},Vector{Float64}}(
    G_ground => [0.0],
    G_excited => [Δ]
)
lehmann_iaj = Dict(
    (G_ground,G_tunnel,G_excited) => ones(ComplexF64,1,1,1,1)
)

# hybridization
Γ = 0.03Δ # thermodynamics
Γ = 0.1Δ  # spectral
tunneling = Dict(
    J => ComplexF64[sqrt(2Γ/π);;]
)

# numerical parameters
cutoff = 500
L = 2.5 # spectral
L = 5.0 # thermodynamics
iterations = 60

# choose what to calculate among:
# - "impurity-shell spectrum"
# - "thermodynamics"
# - "spectral"
run = "thermodynamics"

if run=="impurity-shell spectrum"

    nrgfull( 
        symmetry,
        label,
        L,
        iterations,
        cutoff,
        shell_config,
        tunneling;
        spectrum=spectrum,
        lehmann_iaj=lehmann_iaj,
        until="2-particle multiplets"
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
            max_SJ2=16,
            calculation=calculation,
            spectrum=spectrum,
            lehmann_iaj=lehmann_iaj
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
        max_SJ2=16,
        spectral=true,
        spectrum=spectrum,
        lehmann_iaj=lehmann_iaj
    )

end
