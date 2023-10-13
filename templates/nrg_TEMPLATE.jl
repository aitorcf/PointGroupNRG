#!/usr/bin/env julia

# activate environment where PointGroupNRG is installed
# (if necessary) and load module
environment_pointgroupnrg = "<path_to_environment>"
import Pkg; Pkg.activate(environment_pointgroupnrg)
using PointGroupNRG.NRGCalculator

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
            atom_config ,
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
        atom_config ,
        shell_config ,
        identityrep ,
        epsilon_symparams ,
        u_symparams ,
        hop_symparams ;
        spectral=true
    )


end
