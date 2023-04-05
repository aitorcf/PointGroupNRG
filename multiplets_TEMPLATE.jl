#!/usr/bin/env julia 

# ------------------------------------------------------------------------------
#
# USER INPUT 
#

# Orbital irrep for which to compute the multiplets.
orbital = "Eg"

# Paths and directories.

# main modules directory
moduledir = "/path/to/modules"
# path to clebsch-gordan coefficients for the chosen orbital
cg_path = "/home/acalvo/Bulegoa/ClebschGordan/Oh/cg_symbolic/"
# main directory where the folder for this orbital's irreps are to be stored.
multiplets_dir = "/home/acalvo/Bulegoa/AntiSymmetricPart/Oh/"

# ------------------------------------------------------------------------------


# Create multiplet directory if it does not exist.
multiplets_path = "/home/acalvo/Bulegoa/AntiSymmetricPart/Oh/$(orbital)_julia/"
if !isdir(multiplets_path) 
    mkdir( multiplets_path )
end

# Source multiplet calculator module.
include( "$(moduledir)/multiplets.jl" )

# Compute multiplet states.
compute_asymstates_allN( orbital , cg_path , multiplets_path ; verbose=true )
