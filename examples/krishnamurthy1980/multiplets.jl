#!/usr/bin/env julia 

# load multiplet module
import Pkg; Pkg.activate("../..")
using PointGroupNRG.MultipletCalculator

# spin SU(2)
symmetry = "spin"

# spin
spin = 0.5

# Compute multiplet states.
compute_multiplets( 
    symmetry;
    irrep=spin
)
