#!/usr/bin/env julia 

# load package
package_dir = "."
import Pkg; Pkg.activate(package_dir)
using PointGroupNRG

# path to clebsch-gordan coefficients for the chosen orblets
cg_path = "precompilescripts/clebschgordan"

# orbital irrep for which to compute the multiplets.
orbital = "A1g"

# directory where multiplet folder will be stored
multiplets_path =  "precompilescripts/multiplets"

# Compute multiplet states.
compute_multiplets( orbital , cg_path , multiplets_path ; verbose=true )
