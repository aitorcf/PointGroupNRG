PointGroupNRG is a Julia pacakage designed to perform
Numerical Renormalization Group calculations for 
magnetic impurity Anderson Hamiltonians with finite orbital
symmetries. A comprehensive user's manual and a tutorial
can be found in the `doc/` directory. The presentation of
the code and a comprehensive exposition of the implemented
procedure can be found in [this paper](https://arxiv.org/abs/2307.03658)

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
    julia> ]add "<path_to_PointGroupNRG>"
    julia> ]test PointGroupNRG

(Note that `]` is the prefix to enter `Pkg` mode.) To use
the package, (i) load the environment where it is installed
and (ii) include the `using` sentence:

    julia> activate NRGCalculations
    julia> using PointGroupNRG

# Precompilation
To avoid the latency caused by the loading of the package
and the compilation of functions, two scripts are provided
in order to precompile the necessary functions into sysimage 
files. To use them, run the following commands from the
shell in the package directory (`<PGNRG>`):

    bash> pwd
        <PGNRG>
    > julia precompilescripts/precompile_multiplets.jl
    > julia precompilescripts/precompile_nrgcalculator.jl

These commands will generate the sysimage files
`<PGNRG>/PointGroupNRGMultiplets.so` and
`<PGNRG>/PointGroupNRGCalculator.so`, respectively.
They can then be used to run the multiplet calculations
with 

    > julia -J <PGNRG>/PointGroupNRGMultiplets.so <script>

and NRG calculations with

    > julia -J <PGNRG>/PointGroupNRGCalculator.so <script>
 
with no overhead from package loading and function
compilation. This precompilation process needs to be
repeated if the code is updated or any change is made to it.
For more information, see the [PackageCompiler.jl
documentation](https://julialang.github.io/PackageCompiler.jl/stable/).
