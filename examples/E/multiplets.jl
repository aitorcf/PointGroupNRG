#!/usr/bin/env julia 

# load multiplet module
import Pkg; Pkg.activate("../..")
using PointGroupNRG.MultipletCalculator

# path to clebsch-gordan coefficients for the chosen orbital
clebschgordan_path = "../clebschgordan/O"

# spin-orbital symmetry of the system: O ⊗ SU(2)
symmetry = "PS"

# orbital irrep for which to compute the multiplets.
orbital = "E"

# directory where multiplet folder will be stored
multiplets_path =  "multiplets"

# Compute multiplet states.
compute_multiplets(
    symmetry,
    multiplets_path;
    irrep=orbital,
    clebschgordan_path=clebschgordan_path,
    verbose=false
)
