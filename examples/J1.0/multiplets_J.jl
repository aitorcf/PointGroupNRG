#!/usr/bin/env julia 

# load multiplet module
import Pkg; Pkg.activate("../..")
using PointGroupNRG.MultipletCalculator

# orbital irrep for which to compute the multiplets.
J = 1.0

# directory where multiplet folder will be stored
multiplets_path =  "multiplets"

# Compute multiplet states.
compute_multiplets( 
    "J",
    multiplets_path;
    irrep=J,
    verbose=true
)
