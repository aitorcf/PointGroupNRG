#!/usr/bin/env julia 

# load multiplet module
moduledir = "/path/to/modules/directory"
include( "$(moduledir)/multiplets.jl" )

# orbital irrep for which to compute the multiplets.
orbital = "orbital_irrep"

# directory where multiplet folder will be stored
multiplets_path =  "path/to/multiplet/folder"

# path to clebsch-gordan coefficients for the chosen orbital
cg_path = "path/to/clebsch/gordan/coefficients"

# Compute multiplet states.
compute_asymstates_allN( orbital , cg_path , multiplets_path ; verbose=true )
