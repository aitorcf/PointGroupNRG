#!/usr/bin/env julia 

# load multiplet module
import Pkg; Pkg.activate("../..")
using PointGroupNRG.MultipletCalculator

# path to clebsch-gordan coefficients for the chosen orbital
cg_path = "clebschgordan"

# orbital irrep for which to compute the multiplets.
orbital = "T1g"

# directory where multiplet folder will be stored
multiplets_path =  "multiplets"

# Compute multiplet states.
compute_multiplets( orbital ,
                    cg_path ,
                    multiplets_path ;
                    verbose=true ,
                    doublegroup=true )
