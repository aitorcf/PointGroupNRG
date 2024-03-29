#!/usr/bin/env julia 

# activate environment where PointGroupNRG is installed
# (if necessary) and load module
environment_pointgroupnrg = "<path_to_environment>"
import Pkg; Pkg.activate(environment_pointgroupnrg)
using PointGroupNRG.MultipletCalculator

# path to clebsch-gordan coefficients for the chosen orbital
cg_path = "<path/to/clebsch/gordan/coefficients>"

# orbital irrep for which to compute the multiplets.
orbital = "<orbital_irrep>"

# directory where multiplet folder will be stored
multiplets_path =  "<path/to/multiplet/folder>"

# Compute multiplet states.
compute_multiplets( orbital , cg_path , multiplets_path ; verbose=true )
