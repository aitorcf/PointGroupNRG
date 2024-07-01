#!/usr/bin/env julia 

# load multiplet module
import Pkg; Pkg.activate("../..")
using PointGroupNRG.MultipletCalculator

symmetry = "totalangularmomentum"
J = 1.0

# Compute multiplet states.
compute_multiplets( 
    symmetry;
    irrep=J
)
