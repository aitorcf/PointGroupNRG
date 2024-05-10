#!/usr/bin/env julia 

# load multiplet module
import Pkg; Pkg.activate("../..")
using PointGroupNRG.MultipletCalculator

# path to clebsch-gordan coefficients for the chosen orbital
#cg_path = "../clebschgordan"
clebschgordan_path = "../clebschgordan_nonsimple/"

# spin-orbital symmetry of the system: P ⊗ SU(2)
symmetry = "PS"

# orbital irrep for which to compute the multiplets.
#orbital = "Eg"
orbital = "Eg"

# directory where multiplet folder will be stored
multiplets_path =  "multiplets"

# Compute multiplet states.
compute_multiplets(
    symmetry,
    multiplets_path;
    irrep=orbital,
    clebschgordan_path=clebschgordan_path,
    verbose=true
)
#compute_multiplets( orbital , cg_path , multiplets_path ; verbose=true )
