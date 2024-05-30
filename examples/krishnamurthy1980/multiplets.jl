#!/usr/bin/env julia 

# load multiplet module
import Pkg; Pkg.activate("../..")
using PointGroupNRG.MultipletCalculator

# spin SU(2)
symmetry = "spin"

# spin
S = 0.5

# directory where multiplet folder will be stored
multiplets_path =  "multiplets"

# Compute multiplet states.
compute_multiplets( 
    symmetry,
    multiplets_path;
    irrep=S,
    verbose=true
)
