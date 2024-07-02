#/usr/bin/env julia

# load modules
package_dir = "../.."
import Pkg; Pkg.activate(package_dir)
using PointGroupNRG.NRGCalculator
using Profile, ProfileVega

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
# Cox, original paper
G_ground  = (2,"E",0.0)
G_excited = (1,"E12",0.0)
# Koga 1999
G_ground  = (2,"T2",0.0)
G_excited = (1,"E52",0.0)
spectrum = Dict{Tuple{Int64,String,Float64},Vector{Float64}}(
    G_ground=>[0.0],
    G_excited =>[Δ]
)
lehmann_iaj = Dict{NTuple{3,Tuple{Int64,String,Float64}},Array{ComplexF64,4}}(
    (G_ground,G_hop,G_excited)=>ones(ComplexF64,1,1,1,1)
)

# hybridization
Γ = 0.01 # with L=15, Δ=0.1, and 200 iterations
Γ = 0.02
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
cutoff_type = "multiplet"
cutoff = 200
L = 10.0
iterations = 50

# choose what to calculate
run = "thermo"

if run=="spectrum"

    nrg_full_allsymmetries( 
        symmetry,
        label,
        L,
        iterations,
        cutoff,
        multiplets_dir,
        shell_config,
        tunneling;
        clebschgordan_path=clebschgordan_path,
        identityrep=identityrep,
        spectrum=spectrum,
        lehmann_iaj=lehmann_iaj,
        until="impurity-shell spectrum",
        channels_dos=channels_dos
    )

elseif run=="impurityprojections"

    nrg_full_allsymmetries( 
        symmetry,
        label,
        L,
        iterations,
        cutoff,
        multiplets_dir,
        shell_config,
        tunneling;
        spectrum=spectrum,
        lehmann_iaj=lehmann_iaj,
        compute_impurity_projections=true,
        identityrep=identityrep,
        clebschgordan_path=clebschgordan_path,
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
            clebschgordan_path=clebschgordan_path,
            identityrep=identityrep,
            spectrum=spectrum,
            lehmann_iaj=lehmann_iaj,
            compute_impurity_projections=true,
            print_spectrum_levels=10
        )

    end

elseif run=="convergence"

    for cutoff in collect(100:200:900)
        for calculation in ["CLEAN","IMP"]

            nrg_full_allsymmetries(
                symmetry,
                label,
                L,
                iterations,
                cutoff_type,
                cutoff,
                multiplets_dir,
                shell_config,
                tunneling;
                calculation=calculation,
                clebschgordan_path=clebschgordan_path,
                identityrep=identityrep,
                spectrum=spectrum,
                lehmann_iaj=lehmann_iaj,
                compute_impurity_projections=true,
                print_spectrum_levels=10
            )

        end
        cp(
            "thermodata/thermo_U_diff_z0.0.dat",
            "convergence/m$(cutoff)";
            force=true
        )
    end

elseif run=="gammasweep"

    for Γ in 1.0:1.0:5.0

        tunneling = Dict{String,Matrix{ComplexF64}}(
            "F32" => [sqrt(2Γ/π);;]
        )

        for calculation in ["CLEAN","IMP"]

            nrg_full_allsymmetries(
                symmetry,
                label,
                L,
                iterations,
                cutoff,
                multiplets_dir,
                shell_config,
                tunneling;
                calculation=calculation,
                clebschgordan_path=clebschgordan_path,
                identityrep=identityrep,
                spectrum=spectrum,
                lehmann_iaj=lehmann_iaj,
                compute_impurity_projections=true,
                print_spectrum_levels=10
            )

        end
        cp(
            "thermodata/thermo_U_diff_z0.0.dat",
            "gammasweep/g$(Γ)_thermo";
            force=true
        )
        cp(
            "impurityprojections/imp_proj_U_z0.0.dat",
            "gammasweep/g$(Γ)_improj";
            force=true
        )
    end
end
