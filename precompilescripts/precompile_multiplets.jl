#!/usr/bin/env julia

import Pkg 
Pkg.activate(".")
using PackageCompiler

PackageCompiler.create_sysimage(
    ["PointGroupNRG"]; 
    sysimage_path="PointGroupNRGMultiplets.so",
    precompile_execution_file="precompilescripts/multiplets_A1g.jl"
)
