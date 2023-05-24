#!/usr/bin/env julia 

# load multiplet module
moduledir = "/home/aitor/Bulegoa/PointGroupNRG/modules"
include( "$(moduledir)/multiplets.jl" )

# path to clebsch-gordan coefficients for the chosen orbital
cg_path = "/home/aitor/Bulegoa/PointGroupNRG/examples/clebschgordan"

# orbital irrep for which to compute the multiplets.
orbital = "Eg"

# directory where multiplet folder will be stored
multiplets_path =  "./multiplets"

# Compute multiplet states.
compute_asymstates_allN( orbital , cg_path , multiplets_path ; verbose=true )
