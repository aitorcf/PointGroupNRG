# Overview 

The necessary dependencies can be found in the file
`dependencies.jl`. In order to install them, run

    julia dependencies.jl

A calculation typically consists of four steps.

1. Compute the multiplet states for a given symmetry and
   orbital.
2. Obtain the two-electron atomic multiplets, which will
   help when selecting the Coulomb parameters.
3. Compute atomic spectra with varying parameters until the
   desired ground state is found.
4. Perform NRG calculations with the chosen parameters.

Step 1 is best carried out in a separate script, since it
will not be done more than once. It is recommended that
steps 2, 3, and 4 are integrated into a single script.

## 1. Multiplet calculation

To start with, create a working directory for the
calculations. Copy the file `multiplet_TEMPLATE.jl` to that
directory. The file is a template for running multiplet
calculations for a single orbital irrep, so it has to be run
once for every orbital irrep needed. The main functioncalled at the end of the script, has the form 

     function compute_asymstates_allN( orbital::String ,
                                       cg_path::String , 
                                       multiplets_path::String ;
                                       verbose::Bool=false )

The variables are:

`orbital::String` : Name of the orbital irrep (_e.g._
        "A1g").

`cg_path::String` : Path to the directory where the
Clebsch-Gordan coefficients are stored.

`multiplets_path::String` : Path to the directory where the
multiplet states will be stored.

`verbose::Bool` : `True` for extra verbosity, `False` for
default.

## 2. Two-particle multiplet calculation 

Copy the template file `nrg_TEMPLATE.jl` to the working
directory. The file contains contains the functions
`multiplets_2part`, `atomic_spectrum` and `nrg_full`
corresponding to steps 2, 3 and 4, respectively. In order to
compute the multiplets, we set `run="multiplets"` so that
the code runs the function 

    
    multiplets_2part( cg_o_dir::String , 
                      multiplet_dir::String ,
                      atom_config::Dict{String,Int64} ,
                      identityrep::String ;
                      max_spin2::Int64 )

that computes the multiplets. The input parameters are:

`cg_o_dir::String` : Path to orbital Clebsch-Gordan
coefficients.

`multiplet_dir::String` : Path to multiplet information
(computed in step 1).

`atom_config::Dict{String,Int64}` : Configuration of the
atom, which is given as a dictionary of pairs
    
                Orbital irrep => Number of multiplets

specifying how many one-electron multiplets are there per
orbital irrep. For example, the single-impurity Anderson
model first studied by Wilson has the configuration

        atom_config = Dict{String,Int64}( "A1g" => 1 )

if we use a cubic symmetry group where "A1g" is the
identity.

`identityrep::String` : Identity irrep in the chosen orbital
                        point-group.

`max_spin2::Int64` (optional): Upper bound for $2S$, where
                               $S$ is the total spin of the
                               impurity. If the $2S$ ends up
                               being larger than the chosen
                               value, the calculation will
                               fail. The default value is
                               usually good in this case. 

## 3. Atomic spectrum calculation

We compute the atomic spectrum with the function

    function atomic_spectrum( 
                cg_o_dir::String ,
                asym_dir::String ,
                atom_config::Dict{String,Int64} ,
                identityrep::String ,
                epsilon_symparams::Dict{ String , Vector{ComplexF64} } ,
                u_symparams::Dict{ Tuple{String,Int64} , Matrix{ComplexF64} } ;
                max_spin2::Int64=10 ,
                calculation::String="IMP" )

The new input parameters that we have to add on top of those
needed for step 2 are:

`epsilon_symparams::Dict{String,Vector{ComplexF64}}` :
Occupation energies $\epsilon _r(I)$, where $I$ is the
orbital irrep and $r$ is the outer multiplicity label,
in the format 

        I => [multiplet energies]

where the index of the vector of energies corresponding to
$I$ is the outer multiplicity $r$.

`u_symparams::Dict{Tuple{String,Int64},Matrix{ComplexF64}}`
: Coulomb parameters $U_{r _1 r _2}(I,2S)$, where $I$ is the
orbital irrep, $2S$ is twice the total spin, and $r _1$ and
$r_2$ are outer multiplicity labels. The format is 

        I,2S => Matrix elements

where each matrix element $U _{r_1,r_2}(I,2S)$ with indices $r _1,
r _2$ is given by 

$$
U _{r _1,r _2}(I,2S) = < s_1 | H | s_2 >,
$$

where $|s _1>$ and $|s _2>$ are states belonging to
multiplets with quantum numbers $I$ and $S$ and with outer
multiplicities $r _1$ and $r _2$, respectively. The
definition of this parameter requires knowledge of the
two-particle multiplets, which can be obtained by following
step 2.


## 4. Full NRG calculation

The last step is to perform a full NRG calculation. To start
with, try the main function with only the mandatory input
parameter, which looks like

    function nrg_full( 
                label::String ,
                calculation::String ,
                L::Float64 ,
                iterations::Int64 ,
                cutoff_type::String ,
                cutoff_magnitude::R ,
                cg_o_dir::String ,
                asym_dir::String ,
                atom_config::Dict{String,Int64} ,
                shell_config::Dict{String,Int64} ,
                identityrep::String ,
                epsilon_symparams::Dict{ String , Vector{ComplexF64} } ,
                u_symparams::Dict{ Tuple{String,Int64} , Matrix{ComplexF64} } ,
                hop_symparams::Dict{ String , Matrix{ComplexF64} } )

In this form, the function performs a calculation of
thermodynamic properties which are stored into the
`thermodata` directory, which is created automatically
inside the working directory. The new input parameters are

`label::String` : Name for the system/calculation, which
will be included in the output files.

`calculation::String` : It can be "IMP" for a calculation
with the defined impurity structure, or "CLEAN" for a
calculation without the impurity. If "IMP" is run after
"CLEAN", the thermodynamic contribution of the impurity is
computed by subtracting the "CLEAN" results to the "IMP"
results and storing them as "thermodata/*diff*.dat".

`L::Float64` : Discretization parameter $\Lambda$.

`iterations::Int64` : Number of NRG iterations to perform.

`cutoff_type::String` : It allows to choose between a
`"multiplet"` or an `"energy"` cutoff.

`cutoff_magnitude::R where {R<:Real}` : For
`cutoff="multiplet"`, it specifies the maximum number of
multiplets to be kept (additional multiplets will be kept if
their energy is in the vicinity of the highest-energy
multiplets kept by following this specification). If
`cutoff="energy"`, the maximum energy is set instead (not
suited for $z$-averaged calculations).

`hop_symparams::Dict{String,Matrix{ComplexF64}}` :
Atom-shell hopping amplitudes (hybridization amplitudes)
$V_{r _1 r _2}(\Gamma)$ between the one-particle atomic
multiplet with outer multiplicity $r_1$ and the shell
multiplet with outer multiplicity $r _2$, both belonging to
irrep $\Gamma$. It is given in the format 

        Orbital irrep => [matrix elements;;]





