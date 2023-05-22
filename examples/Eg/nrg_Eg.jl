#!/usr/bin/env julia

# load modules
moduledir = "/home/aitor/Bulegoa/PointGroupNRG/modules"
include( "$(moduledir)/modules.jl" )

# name given to the system
label = "Eg"

# numerical parameters
L = 3.0
cutoff_type = "multiplet"
cutoff_magnitude = 100
iterations = 10

# z averaging 
Nz = 4
Z = generate_Z(Nz)

# clebsch-gordan and multiplet directories
cg_o_dir = "/home/aitor/Bulegoa/PointGroupNRG/examples/clebschgordan_reduced"
multiplets_dir = "multiplets"
identityrep = "A1g"

# configuration of impurity and shells
imp_config = Dict{String,Int64}( "Eg" => 1 )
shell_config = imp_config

# parameters of the hamiltonian
epsilon_symparams = Dict{String,Vector{ComplexF64}}(
    "Eg" => [-0.1]
)
u_symparams = Dict{ Tuple{String,Float64} , Matrix{ComplexF64} }(
    ("Eg",0.0)  => [0.6;;],
    ("A1g",0.0) => [0.5;;],
    ("A2g",1.0) => [0.0;;]
)
hop_symparams = Dict{ String , Matrix{ComplexF64} }(
    "Eg" => [0.5;;]
)

# choose what to calculate
run = "thermozavg"

if run=="multiplets"

    multiplets_2part( 
        cg_o_dir ,
        multiplets_dir ,
        imp_config ,
        identityrep
    )

elseif run=="spectrum"

    atomic_spectrum( 
        cg_o_dir ,
        multiplets_dir ,
        imp_config ,
        identityrep ,
        epsilon_symparams ,
        u_symparams
    )

elseif run=="thermo"

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
            imp_config ,
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
        imp_config ,
        shell_config ,
        identityrep ,
        epsilon_symparams ,
        u_symparams ,
        hop_symparams ;
        spectral=true
    )

elseif run=="thermozavg"


    for calculation in ["CLEAN","IMP"],
        z in generate_Z(Nz)

        nrg_full( 
            label ,
            calculation ,
            L ,
            iterations ,
            cutoff_type ,
            cutoff_magnitude ,
            cg_o_dir ,
            multiplets_dir ,
            imp_config ,
            shell_config ,
            identityrep ,
            epsilon_symparams ,
            u_symparams ,
            hop_symparams
        )

    end

    zavg_thermo( label , Z )

end
