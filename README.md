# Overview
PointGroupNRG is a Julia pacakage designed to perform
Numerical Renormalization Group calculations for 
magnetic impurity Anderson Hamiltonians with finite point
or double group symmetries. A comprehensive user's manual and a tutorial
can be found in the `doc/` directory. The presentation of
the code and a comprehensive exposition of the implemented
procedure can be found in [this paper](https://arxiv.org/abs/2307.03658)
and [NEW ARXIV].

# Changes
- MERGE DATE
    - Double group symmetries and total angular momentum conservation
    available.
    - Modified documentation:
        - `doc/Manual.md` replaced by `manual/manual.pdf`.
        - `doc/Tutorial.md` removed.
    - Interface changed, see `manual/manual.pdf`.
    - New examples in `examples/`.
    - Precompilation scripts in `precompile/` removed.

# Installation
PointGroupNRG is not in the General Registry, so it cannot
be installed using the standard  

    julia> using Pkg; Pkg.add("PointGroupNRG")

procedure. Instead, the code has to be downloaded/cloned. To
install it, (i) create an environment where `PointGroupNRG`
and all its depenencies are going to be installed (optional,
but recommended), (ii) add the `PointGroupNRG`
package, and (iii) test it (also optional, but recommended).
As an example, one could do the following:

    bash> mkdir NRGCalculations
    bash> cd NRGCalculations
    bash> julia

    julia> ]activate .
    julia> ]add <path_to_PointGroupNRG>
    julia> ]test PointGroupNRG

(Note that `]` is the prefix to enter `Pkg` mode.) To use
the package, (i) load the environment where it is installed
and (ii) include the `using` sentence:

    julia> activate NRGCalculations
    julia> using PointGroupNRG

# Precompilation
To avoid the package loading latency, it is possible to
create a sysimage using the `PackageCompiler` package. 
See the [PackageCompiler.jl documentation](https://julialang.github.io/PackageCompiler.jl/stable/) and, in particular, the [sysimage section](https://julialang.github.io/PackageCompiler.jl/stable/sysimages.html).
