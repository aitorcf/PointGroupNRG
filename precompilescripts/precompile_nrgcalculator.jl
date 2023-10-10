#!/usr/bin/env julia

import Pkg 
Pkg.activate(".")
using PackageCompiler

PackageCompiler.create_sysimage(
    ["PointGroupNRG"]; 
    sysimage_path="PointGroupNRGCalculator.so",
    precompile_execution_file="precompilescripts/nrg_A1g.jl"
)
