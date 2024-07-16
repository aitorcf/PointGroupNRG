# Overview
PointGroupNRG is a Julia package designed to perform
Numerical Renormalization Group calculations for magnetic
impurity Anderson Hamiltonians with several charge, orbital
and spin symmetries, including finite point and double group
symmetries. A comprehensive user's manual can be found in
`manual/`. The presentation of the code and a comprehensive
exposition of the implemented procedure can be found in
[this paper](https://arxiv.org/abs/2307.03658) and [NEW
ARXIV].

# Compatibility
- Tested and working for Julia 1.10.0
- Issues have been found in Julia 1.8.x

# Changes
- 1/7/2024
    - Double group symmetries and total angular momentum conservation
    available.
    - Modified documentation:
        - `doc/Manual.md` replaced by `manual/manual.pdf`.
        - `doc/Tutorial.md` removed.
    - Interface changed, see `manual/manual.pdf`.
    - New examples in `examples/`.
    - Precompilation scripts in `precompile/` removed.
    - Working version with the functionality described in [this paper](https://arxiv.org/abs/2307.03658)
    moved to branch `v1.0`.

# Installation
PointGroupNRG is not in the General Registry. To use it from
your own Julia environment, download `PointGroupNRG` and
add it to your project:

    pkg> activate <path_to_project> 
    pkg> add <path_to_PointGroupNRG>

# Examples
Examples for several systems are provided in the `examples/`
directory. To run them, go to the corresponding directory,
_e.g._ `examples/krishnamurthy1980`, and run the scripts:

    $ cd <...>/PointGroupNRG&/examples/krishnamurthy1980
    $ ./multiplets.jl

This computes the multiplet states to be used in the
NRG calculation. Afterwards, edit the `nrg.jl` script to set
the variable `run` to the desired value, perform any other
modifications, and run the script:

    $ julia
    julia> include("nrg.jl")

Loading the script within a Julia session allows to skip
package loading time for subsequent runs. Check the output
in the generated directories (see `manual/manual.pdf`).

# Precompilation
To avoid the package loading latency, it is possible to
create a sysimage using the `PackageCompiler` package. 
See the [PackageCompiler.jl documentation](https://julialang.github.io/PackageCompiler.jl/stable/) and, in particular, the [sysimage section](https://julialang.github.io/PackageCompiler.jl/stable/sysimages.html).
