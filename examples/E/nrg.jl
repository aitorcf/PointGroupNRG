#!/usr/bin/env julia

# Equivalent to Eg in https://doi.org/10.1016/j.cpc.2023.109032

# load modules
package_dir = "../.."
import Pkg; Pkg.activate(package_dir)
using PointGroupNRG.NRGCalculator

# label for the system
label = "E"

# symmetry
symmetry = "pointspin"
clebschgordan_path = "../clebschgordan/O"
identityrep = "A1"

# system configuration
impurity_config = Dict{String,Int64}( "E" => 1 )
shell_config = impurity_config

# standard anderson hamiltonian
ϵ = -0.1
u = 0.3
j = 0.2
onsite = Dict{String,Vector{ComplexF64}}(
    "E" => [-0.1]
)
interaction = Dict{Tuple{String,Float64},Matrix{ComplexF64}}(
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
Γ = 0.01
tunneling = Dict{String,Matrix{ComplexF64}}(
    "E" => [sqrt(2Γ/π);;]
)

# numerical parameters
cutoff = 500
L = 5.0
iterations = 100

# choose what to calculate among:
# - "2-particle multiplets"
# - "impurity spectrum"
# - "impurity-shell spectrum"
# - "thermodynamics"
# - "thermoionic"
run = "thermoionic"

if run=="2-particle multiplets"

    nrgfull(
        symmetry,
        label,
        L,
        iterations,
        cutoff,
        shell_config,
        tunneling;
        clebschgordan_path=clebschgordan_path,
        identityrep=identityrep,
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
        clebschgordan_path=clebschgordan_path,
        identityrep=identityrep,
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
        clebschgordan_path=clebschgordan_path,
        identityrep=identityrep,
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
            clebschgordan_path=clebschgordan_path,
            identityrep=identityrep,
            impurity_config=impurity_config,
            onsite=onsite,
            interaction=interaction,
            compute_impurity_projections=true
        )

    end

elseif run=="thermoionic"

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
            clebschgordan_path=clebschgordan_path,
            identityrep=identityrep,
            spectrum=spectrum,
            lehmann_iaj=lehmann_iaj,
            compute_impurity_projections=true
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
        spectral=true,
        clebschgordan_path=clebschgordan_path,
        identityrep=identityrep,
        impurity_config=impurity_config,
        onsite=onsite,
        interaction=interaction,
        compute_impurity_projections=true
    )

end
