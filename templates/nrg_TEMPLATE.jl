#!/usr/bin/env julia

# load modules
moduledir = "/path/to/modules/"
include( "$(moduledir)/modules.jl" )

# choose what to calculate
run = "run_choice"

if run=="multiplets"

    multiplets_2particles( 
        cg_o_dir ,
        multiplets_dir ,
        atom_config ,
        identityrep
    )

elseif run=="spectrum"

    impurity_spectrum( 
        cg_o_dir ,
        multiplets_dir ,
        atom_config ,
        identityrep ,
        epsilon_symparams ,
        u_symparams
    )

elseif run=="nrg"

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
        hop_symparams
    )


end
