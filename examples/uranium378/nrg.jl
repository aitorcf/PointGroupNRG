#!/usr/bin/env julia

# Two multiplet uranium model introduced in Phys. Rev. Lett. 59, 1240

# load modules
package_dir = "../.."
import Pkg; Pkg.activate(package_dir)
using PointGroupNRG.NRGCalculator
using Distributed

# label for the system
label = "U"

# symmetry
symmetry = "D"
identityrep = "A1"
clebschgordan_path = "../clebschgordan/O"

# multiplet states
multiplets_dir = "multiplets"

# system configuration
shell_config = Dict{String,Int64}( "F32"=>1 )

# ionic anderson model
Δ = 0.1
G_hop     = (1,"F32",0.0)
G_ground  = (2,"E",0.0)
G_excited = (1,"E12",0.0)
spectrum = Dict{Tuple{Int64,String,Float64},Vector{Float64}}(
    G_ground=>[0.0],
    G_excited =>[Δ]
)
lehmann_iaj = Dict{NTuple{3,Tuple{Int64,String,Float64}},Array{ComplexF64,4}}(
    (G_ground,G_hop,G_excited)=>ones(ComplexF64,1,1,1,1)
)

# hybridization
Γ = 0.01
tunneling = Dict{String,Matrix{ComplexF64}}(
    "F32" => [sqrt(2Γ/π);;]
)

# channel DOS
f(x) = 0.5#+0.01x
channels_dos = Dict{String,Vector{Function}}(
    "F32" => [f]
)

# numerical parameters
cutoff = 200
L = 10.0
iterations = 50

# choose what to calculate among:
# - "impurity-shell spectrum"
# - "thermodynamics"
# - "thermozavg"
#       *compare results with the same small cutoff
#       with and without averaging
run = "thermozavg"

if run=="impurity-shell spectrum"

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
        spectrum=spectrum,
        lehmann_iaj=lehmann_iaj,
        until="impurity-shell spectrum",
        channels_dos=channels_dos
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
            spectrum=spectrum,
            lehmann_iaj=lehmann_iaj,
            compute_impurity_projections=true
        )

    end

elseif run=="thermozavg"

    # number of z values (= number of parallel processes)
    Nz = 4

    # generate z values
    Z = generate_Z(Nz)

    # generate one (worker) process per z value
    addprocs(Nz)

    # define variables and load package on all processes
    @everywhere begin
        using PointGroupNRG.NRGCalculator
        symmetry=$symmetry
        label=$label
        L=$L
        iterations=$iterations
        cutoff=$cutoff
        shell_config=$shell_config
        tunneling=$tunneling
        clebschgordan_path=$clebschgordan_path
        identityrep=$identityrep
        spectrum=$spectrum
        lehmann_iaj=$lehmann_iaj
    end

    # parallel loop with @distributed 
    # @sync prevents from asyncronously continuing with the script
    @sync @distributed for z in Z
        for calculation in ["CLEAN","IMP"]

            nrgfull(
                symmetry,
                label,
                L,
                iterations,
                cutoff,
                shell_config,
                tunneling;
                z=z,
                calculation=calculation,
                clebschgordan_path=clebschgordan_path,
                identityrep=identityrep,
                spectrum=spectrum,
                lehmann_iaj=lehmann_iaj,
                compute_impurity_projections=true
            )

        end
    end

    # average over results for various z
    zavg_thermo(label,Z)

    # remove created processes. necessary when loading the script
    # from a julia session with include("nrg.jl") to avoid process
    # overpopulation.
    rmprocs(workers())
end
