# Overview
PointGroupNRG is a Julia package designed to perform
Numerical Renormalization Group calculations for magnetic
impurity Anderson Hamiltonians with several charge, orbital
and spin symmetries, including finite point and double group
symmetries. A comprehensive user's manual can be found in
`manual/`. The presentation of the code and a comprehensive
exposition of the implemented procedure can be found in
[this paper](https://arxiv.org/abs/2307.03658) and [this paper](https://doi.org/10.48550/arXiv.2409.12050).

# Compatibility
- Tested and working for Julia 1.10.0
- Issues have been found in Julia 1.8.x

# Changes

Specific stable versions vX.Y.Z of the code are contained in
branches. The `master` branch contains the latest changes.

#### v1.0.0 (1/7/2024)
- Working version with the functionality described in [this paper](https://arxiv.org/abs/2307.03658)

#### master (1/7/2024)
- Double group symmetries and total angular momentum conservation
available.
- Modified documentation:
    - `doc/Manual.md` replaced by `manual/manual.pdf`.
    - `doc/Tutorial.md` removed.
- Interface changed, see `manual/manual.pdf`.
- New examples in `examples/`.
- Precompilation scripts in `precompile/` removed.
- Incompatible with v1.0.0.

#### v2.0.0 (4/8/2024):
- Stable version with functionality introduced in
[1/7/2024](#master-(1/7/2024)).

# Installation
PointGroupNRG is not in the General Registry. To use it from
your own Julia environment, download or clone `PointGroupNRG` and
add it to your project:

    pkg> activate <path_to_project> 
    pkg> add <path_to_PointGroupNRG>

The `pkg>` mode is entered by typing `]` in the `julia>`
shell. If a new version is pulled from GitHub, for the
changes to apply it is necessary to update the added
PointGroupNRG package in the local project using

    pkg> up

Instead of adding the package from a local project, the
downloaded/cloned repository can be activated and
instantiated:

    pkg> activate <path_to_PointGroupNRG>
    pkg> instantiate # only the first time

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
