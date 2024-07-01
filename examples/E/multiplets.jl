#!/usr/bin/env julia 

# load multiplet module
import Pkg; Pkg.activate("../..")
using PointGroupNRG.MultipletCalculator

clebschgordan_path = "../clebschgordan/O"
symmetry = "pointspin"
orbital = "E"

# Compute multiplet states.
compute_multiplets(
    symmetry;
    irrep=orbital,
    clebschgordan_path=clebschgordan_path
)
