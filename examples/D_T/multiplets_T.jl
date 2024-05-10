#!/usr/bin/env julia 

# load multiplet module
import Pkg; Pkg.activate("../..")
using PointGroupNRG.MultipletCalculator

# path to clebsch-gordan coefficients for the chosen orbital
clebschgordan_path = "/home/aitor/Bulegoa/ClebschGordan/D_T/cg_symbolic"

# orbital irrep for which to compute the multiplets.
orbital = "T"

# directory where multiplet folder will be stored
multiplets_path =  "multiplets"

# Compute multiplet states.
compute_multiplets( 
    "D",
    multiplets_path;
    orbital=orbital,
    clebschgordan_path=clebschgordan_path,
    verbose=true
)
