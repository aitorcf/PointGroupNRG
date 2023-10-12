#!/usr/bin/env julia

import Pkg 
Pkg.activate(".")
using PackageCompiler

include("multiplets_A1g.jl")

PackageCompiler.create_sysimage(
    ["PointGroupNRG"]; 
    sysimage_path="PointGroupNRGCalculator.so",
    precompile_execution_file="precompile/nrg_A1g.jl"
)
rm("precompile/multiplets";recursive=true)
rm("thermodata";recursive=true)
rm("spectral";recursive=true)
