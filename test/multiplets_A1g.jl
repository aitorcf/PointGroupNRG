#!/usr/bin/env julia 

# load package
using PointGroupNRG.MultipletCalculator

# path to clebsch-gordan coefficients for the chosen orblets
cg_path = "clebschgordan"

# orbital irrep for which to compute the multiplets.
orbital = "A1g"

# directory where multiplet folder will be stored
multiplets_path = "multiplets"

# Compute multiplet states.
compute_multiplets( orbital , cg_path , multiplets_path ; verbose=true )
