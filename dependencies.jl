#!/usr/bin/env julia

using Pkg

# formatted printing 
Pkg.add("Prinft")
# file management
Pkg.add("Glob") 

# I/O
Pkg.add("DelimitedFiles")

# vector-matrix operations
Pkg.add("LinearAlgebra")

# basis generation
Pkg.add("Combinatorics")

# spin Clebsch-Gordan coefficients
Pkg.add("PartialWaveFunctions")

# interpolations for the spectral function
Pkg.add("Interpolations")

# parallelization
Pkg.add("Distributed")

