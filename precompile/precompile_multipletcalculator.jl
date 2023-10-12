#!/usr/bin/env julia

import Pkg 
Pkg.activate(".")
using PackageCompiler

create_sysimage(
    ["PointGroupNRG"]; 
    sysimage_path="PointGroupNRGMultiplets.so",
    precompile_execution_file="precompile/multiplets_A1g.jl"
)

rm("precompile/multiplets";recursive=true)
