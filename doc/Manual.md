# Table of contents 
1. [Hamiltonian and Symmetry](#hamiltonian)
2. [Main functions and parameters](#functions)
3. [Orbital Clebsch-Gordan coefficients](#clebsch-gordan)
4. [Thermodynamic calculations](#thermodynamics)
5. [Spectral function calculations](#spectral)
6. [Averaging over the twisting parameter $z$](#z-averaging)

    6.1. [Parallelization](#parallelization)

7. [Template scripts and workflow](#templates)
8. [Examples and tutorial](#examples)


# Hamiltonian and Symmetry <a name="hamiltonian"></a>
The PointGroupNRG code is designed to solve Anderson impurity
models with symmetry 

```math
U(1)_{\text{C}} \otimes P_{\text{O}} \otimes SU(2)_{\text{S}},
```

where $\`U(1)_{\text{C}}\`$ corresponds to particle conservation, $P_{\text O}$
is the a finite orbital point-group, and $SU(2)_{\text S}$
represents spin isotropy.

The general Anderson Hamiltonian is

$$
\begin{equation}
H = H_{o} + H_{C} + H_{h} + H_{c},
\end{equation}
$$

where 

$$
\begin{equation}
H_{o} = \sum_{\alpha} \epsilon_\alpha 
f^\dagger_\alpha f_\alpha
\end{equation}
$$

contains the occupation energies,

$$
\begin{equation}
H_{C} = \sum_{\alpha_{1} \alpha_{2} \alpha_{3} \alpha_{4}} 
U_{\alpha_{1} \alpha_{2} \alpha_{3} \alpha_{4}}
f^\dagger_{\alpha_{1}} f^\dagger_{\alpha_{2}}
f_{\alpha_{3}} f_{\alpha_{4}}
\end{equation}
$$

is the Coulomb interaction,

$$
\begin{equation}
H_{h} = \int d\epsilon\sum_{\alpha \beta} V_{\alpha\beta}(\epsilon) f^\dagger_\alpha
c_{\beta \epsilon}
\end{equation}
$$

is the hybridization term, and 

$$
\begin{equation}
H_{c} = \int 
\epsilon
c^\dagger_{\beta \epsilon}
c_{\beta \epsilon}
d\epsilon  
\end{equation}
$$

is the energy of the conduction electrons. In order to ensure that the Hamiltonian has the 
desired symmetry, the parameters $\epsilon_\alpha$,
$U_{\alpha_{1}\alpha_{2}\alpha_{3}\alpha_{4}}$ and
$V_{\alpha\beta}(\epsilon)$ have to fullfil certain
symmetry restrictions. (An in-depth
review of the model and the notation used can be found in
[REF_TO_PAPER]).The code provides a tool for
constructing a Hamiltonian with those restrictions for a
given crystal point-group and solving it using the NRG
method.  


# Main functions and parameters <a name="functions"></a>

The code provides four functions to the user, one for each
step in a typical NRG calculation:

1. `compute_multiplets`<a name="step_1"></a>: Compute the many-particle multiplet
   states composed of electrons in states belonging to the
   same orbital irrep.
2. `multiplets_2particles`<a name="step_2"></a>: Compute the two-particle
   impurity multiplet states for the defined system. This
   information helps when defining the parameters of the
   impurity Anderson Hamiltonian.
3. `impurity_spectrum`<a name="step_3"></a>: Compute the impurity spectrum.
4. `nrg_full`<a name="step_4"></a>: Perform a complete NRG
   calculation and write the thermodynamic or spectral
   functions to files.

These functions are described in detail in the following.

## `compute_multiplets`

This function is defined in the `modules/multiplets.jl`
module. Given the appropriately formatted orbital
Clebsch-Gordan coefficients as input (see the corresponding
[section](#clebsch-gordan)), it computes the
multiplet states for a given orbital irrep and stores them
in a directory created in the desired path. This directory
is then read by the other functions. This function has to be
run for each orbital irrep present in the system: for
example, in an impurity with orbital irreps $E_g$ and
$T_{2g}$, the function must be run once for $E_g$ and once
for $T_{2g}$. 

The function is defined as follows:

    compute_multiplets( orbital::String ,
                        cg_o_dir::String ,
                        multiplets_path::String ;
                        verbose::Bool=false )

The mandatory input parameters are:

* `orbital::String`: Name of the orbital irrep for which the
multiplets are going to be computed. It has to coincide with
the name of the orbital irrep in the Clebsch-Gordan
coefficient files (see [here](#clebsch-gordan)).
* `cg_o_dir::String`<a name="cg_o_dir"></a>: Path to the
orbital Clebsch-Gordan coefficient files.
* `multiplets_path` <a name="multiplets_path"></a>: Path
where the information about the multiplet states is going to
be stored. A folder with the appropriate name will be
created in that path. If the path does not exist, it will be
created.

The optional parameter is:

* `verbose::Bool=false`: Extra verbose output if `true`,
    normal verbosity if `false`. Mostly for testing
    purposes.

## `multiplets_2particles` 

This function is defined in the `modules/automatization.jl`
module. Given the Clebsch-Gordan coefficients and the
one-irrep multiplet states computed in the previous step, it
computes all the impurity multiplet states for a given
impurity configuration. 

The function is defined as follows:

    multiplets_2particles( cg_o_dir::String ,
                           multiplets_path::String ,
                           impurity_config::Dict{String,Int64} ,
                           identityrep::String ;
                           max_spin2::Int64=10 )

The mandatory input parameters are:

* `cg_o_dir::String`: see [here](#cg_o_dir).
* `multiplets_path::String`: see [here](#multiplets_path)
* `impurity_config::Dict{String,Int64}`<a name="impurity_config"></a>: 
For each orbital irrep, how many orbitals levels there are
belonging to that same irrep. For example, in a system with
two impurities, each with an $s$ orbital that we represent
with the irrep $A_{1g}$ of the cubic group $\mathcal O_{h}$,
we would have:

        impurity_config = Dict{String,Int64}( "A1g" => 2 )

* `identityrep::String`<a name="identityrep"></a>: 
The name of the identity orbital irrep. For example, for the
orbital group $\mathcal O_{h}$ it is $A_{1g}$, so
`identityrep="A1g"`.

Optional parameter:

* `max_spin2::Int64=10`<a name="max_spin2"></a>: 
Maximum value of twice the total spin, $2S$, for which the
spin Clebsch-Gordan coefficients are going to be computed.
If the value is too low, the program will likely run into
`NaN` problems as a consequence of segmentation faults,
since it will try to access out-of-bounds array
elements. If the value is too large, storing all the
coefficients will require a lot of memory. This is
especially relevant for the `nrg_full` function, which
also uses this variable and is very
performance-sensitive.

## `impurity_spectrum`

This function is defined in the `modules/automatization.jl`
module. It contructs and solves the impurity Hamiltonian,
giving its spectrum as output. The main goal of this
function is to provide a way to investigate the impurity
without having to run a complete NRG calculation.

The function is defined as follows:

    impurity_spectrum( cg_o_dir::String ,
                       multiplets_dir::String ,
                       impurity_config::Dict{String,Int64} ,
                       identityrep::String ,
                       epsilon_symparams::Dict{ String , Vector{ComplexF64} } ,
                       u_symparams::Dict{ Tuple{String,Float64} , Matrix{ComplexF64} } ;
                       max_spin2::Int64=10 )

The mandatory input parameters are:
* `cg_o_dir::String`: Described [here](#cg_o_dir).
* `multiplets_dir::String`: Described [here](#multiplets_dir).
* `impurity_config::Dict{String,Int64}`: See
[here](#impurity_config).
* `identityrep::String`: Described [here](#identityrep).
* `epsilon_symparams::Dict{String,Vector{Float64}}` 
<a = name="epsilon_symparams"></a>:
Symmetry-adapted occupation energy parameters
$\epsilon_{\alpha}=\epsilon_{r_\alpha}(\Gamma_\alpha)$ of
the impurity Anderson Hamiltonian. The keys are the names of
the one-electron orbital irreps $I_{\alpha}$, which
completely determine the irrep
$\Gamma_{\alpha}=(N_\alpha=1,I_\alpha,S_\alpha=1/2)$. The
values are vectors of occupations energies for each
one-electron multiplet belonging to irrep $\Gamma_\alpha$
and with outer multiplicity $r_\alpha$: the $r_\alpha$-th
element of each of these vectors is the occupation energy of
the corresponding multiplet. Each dictionary entry can be
expressed as follows:

$$
\Gamma \implies \epsilon_{r}(\Gamma)
$$
    
                                        I::String => eps_r::Vector{Float64}

* `u_symparams::Dict{Tuple{String,Float64},Matrix{ComplexF64}}`
<a name="u_symparams"></a>:
Symmetry-adapted Coulomb parameters
$U_{r_{\alpha}r_{\alpha}'}(\Gamma_\alpha)$of the impurity Anderson
Hamiltonian (see PAPER). The keys are tuples containing the
irreps $I_\alpha$ and the spin irreps $S_\alpha$ (total
spin) of the two-particle impurity multiplets. The values
are matrices with indices $r_\alpha$ and $r_\alpha'$. The
dictionary entries can be expressed accordingly as follows:

$$
(I,S) \implies U_{r r'}(\Gamma)
$$

                    (I,S)::Tuple{String,Float64} => U_rrp::Matrix{ComplexF64}

The optional argument is:
* `max_spin2::Int64`: Described [here](#max_spin2).


## `nrg_full`

This function is defined in the `modules/automatization.jl`
module. It is the main function provided by the code. Its
aim is to perform a complete NRG calculation: (i) build
the impurity Hamiltonian and solve it, (ii) couple the
Hamiltonian to the conduction channels, and (iii) perform
the NRG iterations while computing all the information
needed for giving the thermodynamic curves or the spectral
function.

The function is defined as follows:

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
        u_symparams::Dict{ Tuple{String,Float64} , Matrix{ComplexF64} } ,
        hop_symparams::Dict{ String , Matrix{ComplexF64} } ;
        distributed::Bool=false,
        z::Float64=0.0 ,
        max_spin2::Int64=10 ,
        channel_etas::Dict{ String , Vector{Function} }=Dict{ String , Vector{Function} }() ,
        discretization::String="standard" ,
        mine::Float64=0.0 ,
        betabar::Float64=1.0 ,
        spectral::Bool=false ,
        spectral_broadening::Float64=1.0 ,
        K_factor=2.0 ,
        orbitalresolved::Bool=false,
        compute_impmults::Bool=false 
    ) where {R<:Real}

The mandatory input parameters are:

* `label::String`<a name="label"></a>: 
Name given to the system. It is used as part of the name of
the output files containing data about thermodynamic or
spectral functions.
* `calculation::String`: It defines whether to perform a
calculation with the impurity (`calculation="IMP"`) or
without the impurity (`calculation="CLEAN"`). The latter
option is used in the thermodynamic calculations in order to
subtract the results to those obtained with the impurity and
thereby isolate the impurity contribution to the
thermodynamics of the system.
* `L::Float64`: Discretization parameter $\Lambda$. It is
usually chosen between $\Lambda=2$ and $\Lambda=3$, although
values as large as $\Lambda=10$ can be used for
thermodynamic calculations with $z$-averaging (see
[corresponding section](#z-averaging)).
* `iterations::Int64`: Number of NRG iterations to perform.
* `cutoff_type::String`: Choose whether to keep a given
number of multiplets per iteration (`cutoff_type="multiplet"`)
or to keep all the multiplets with eigenenergy below an
energy cutoff (`cutoff_type="energy"`).
* `cutoff_magnitude::R where {R<:Real}`: For
`cutoff_type="multiplet"`, the number of multiplets to keep
at each iteration. For `cutoff_type="energy"`, the cutoff
energy. The code automatically keeps multiplets with
eigenenergy within a small range above the cutoff so as to
avoid breaking accidental degeneracies.
* `cg_o_dir::String`: Described [here](#cg_o_dir).
* `multiplets_dir::String`: Described [here](#multiplets_dir).
* `impurity_config::Dict{String,Int64}`: Described [here](#impurity_config)
* `shell_config::Dict{String,Int64}`: Same as `impurity_config`, but
for the conduction shell sites (every conduction shell site has the
same configuration). For example, in we have a one-orbital
impurity connected to a single conduction channel,
we would have `shell_config=impurity_config`. For an
impurity with one $s$ orbital connected to two $s$ channels,
we would have 

    impurity_config = Dict{String,Int64}( "A1g" => 1 )
    shell_config    = Dict{String,Int64}( "A1g" => 2 )

* `identityrep::String`: Described [here](#identityrep).
* `epsilon_symparams::Dict{String,Vector{ComplexF64}}`: See 
[here](#epsilon_symparams).
* `u_symparams::Dict{Tuple{String,Float64},Matrix{ComplexF64}}`: 
Described [here](#u_symparams).
* `hop_symparams::Dict{String,Matrix{ComplexF64}}`:
Hybridization amplitudes (hopping parameters)
$V_{\alpha\beta}=V_{r_\alpha r_\beta}(\Gamma_\alpha)$, where
$\Gamma_\alpha=\Gamma_\beta$ is a one-particle irrep
completely determined by the orbital irrep $I_\alpha$. The
dictionary entries can be expressed as follows:

$$
I \implies V_{r r'}(\Gamma)
$$

            I::String => V_rrp::Matrix{ComplexF64}

The optional input parameters are:

* `distributed::Bool=false`: Whether to perform the
diagonalization at each NRG step in parallel
(`distributed=true`) or not (`distributed=false`).
Parallelization is implemented using the `@distributed` 
macro provided by the `Distributed` package. As of
now, serial calculations have been found to be faster 
for every test performed and therefore it has not been
extensively used or tested. It is adivised to leave it to
its default value `distributed=false`.
* `z::Float64=0.0`<a name="z"></a>: 
Value of the twisting parameter $z$ used for averaging over
various discretization grids in order to improve smoothness
and resolution in thermodynamic and spectral calculations,
and also to reduce overbroadening effects in the latter.
See [corresponding section](#z-averaging).
* `max_spin2::Int64=10`: Described [here](#max_spin2).
* `channel_etas::Dict{String,Vector{Function}}=Dict{String,Vector{Function}}()`:
Parameters $\gamma_{r_\alpha}(\Gamma_\alpha;\epsilon)$ 
that define energy-dependent hybridization functions (see parameter 
$\gamma(\epsilon)$ used by Campo and Oliveira in Phys. Rev.
B 72, 104432). The dictionary has a similar structure as
`hop_symparams`, the differences being that (i) the vectors
contain functions instead of numbers and (ii) we have a
single index $r_\alpha$ instead of double indices $r_\alpha
r_\alpha'$ because the energy-dependent hybridization is
only implemented for systems where each orbital is coupled
to one channel at most. The default value results in
constant functions $\gamma(\epsilon)=1/2$ for every channel.
If a value other than the default is used for this
parameter, choose `discretization="lanczos"` (see
[`discretization`](#discretization))
* `discretization::String="standard"`: Discretization
method. If `discretization="standard"`, the original scheme
proposed by Wilson (Phys. Rev. B 21, 1003) is used and the
rescaled channel hopping amplitudes are computed using the
analytic formula. If `discretization="lanczos"`, a
matrix Lanczos algorith is used to compute the hoppings for
the first iterations and the rest are computed using the
asymptotic formula (see Phys. Rev. B 41, 9403). 
* `mine::Float64=0.0`: If `cutoff_type="multiplet"`, `mine`
sets a minimum energy of the multiplet with the largest
eigenenergy, in such a way that the multiplet cutoff is
augmented as necessary to meet this requirement.
* `betabar::Float64=1.0`<a name="betabar"></a>: Parameter $\bar\beta$ used in
thermodynamic calculations (see Phys. Rev. B 21, 1003). It
defines the temperature for each NRG iteration.
* `spectral::Bool=false`: Whether to compute spectral
functions (`spectral=true`) or not (`spectral=false`).
* `spectral_broadening::Float64=1.0`<a
name="broadening"></a>: 
Broadening factor applied to the Gaussian broadening of
spectral functions. The default value is recommended unless
the $z$-averaging is used (see [corresponding
section](#z-averaging)), in which case the
broadening should be around $1/N_{z}$, where $N_{z}$ is the
number of $z$ values. This reduction allows to reduce
overbroadening errors.
* `K_factor::Float64=2.0`<a name="K"></a>: 
Parameter $K$ used in spectral function calculations to
define the energy $\omega$ for which to compute the spectral
function at each NRG iteration, which is defined as
$\omega=K\omega_{N}$, where $\omega_{N}$ is the energy scale
associated to the $N$-th iteration.
* `orbitalresolved::Bool=false`<a name="orbitalresolved"></a>: 
Whether to compute orbital-resolved spectral functions
(`orbitalresolved=true`) or not
(`orbitalresolved=false`). In orbital-resolved
calculations, the spectral function is computed
separately for excitations belonging to each
one-electron impurity multiplet.
* `compute_impmults::Bool=false`: Whether to compute the
thermodynamic weights $W_{m}(T)$ for each impurity multiplet
$m$. This quantity is defined as

$$
W_{m}(T) = \text{Tr}\{ D_{m} P_{m} e^{-H/k_{B} T}\},
$$

where the $P_{m}$ is a projector onto the impurity multiplet
space, defined as 

$$
P_{m} = \sum_{\psi\in m} \sum_{c}
|\psi;c\rangle \langle\psi;c|,
$$

where the sums run over projections onto states where
the impurity is in a state $|\psi\rangle$ belonging to the
impurity multiplet $m$ and the conduction electrons are in a
state $|c\rangle$.

# Orbital Clebsch-Gordan coefficients <a name="clebsch-gordan"></a>

The directory of Clebsch-Gordan coefficients, specified by
the variable [`cg_o_dir`](#cg_o_dir), must contain files
named `NxM_IxJ.tex`. These files contain information about 
the decomposition of $I\otimes J$, the product of orbital
irreps `I` and `J` with associated indices `N` and `M`,
respectively (the names of the irreps are defined by the
user). It is sufficient to provide files for the
irrep products that will be used in the calculation. A
single file has to be provided for each orbital irrep pair:
if there is a file for $I\otimes J$, there must be no file
for $J\otimes I$.

The content of each file is organized in blocks separated by
blank lines. Each block contains the coefficients of the
states belonging to an irrep $K$ in the product $I\otimes J$,
_i.e._ $I\otimes J = ... \oplus K \oplus \dots$, 
in terms of states belonging to irreps $I$ and $J$ (these
are the  Clebsch-Gordan coefficients). The first line of the
block contains the index of $K$ and its name, and the
following lines contain the Clebsch-Gordan coefficients. If
$i$, $j$ and $k$ are the partner labels of states belonging
to $I$, $J$ and $K$, respectively, then the Clebsch-Gordan
coefficients are formatted as 

    ( i j | k ) = c

where `c` is the value of the coefficient.

As an example, consider a file `2x2_EgxEg.txt` containing
the Clebsch-Gordan coefficients for the decomposition of the
irrep product $E_{g}\otimes E_{g}$ involving the irrep $E_{g}$
of the cubic group $O_{h}$. Then the contents of the file
can be

    0 A1g 
    ( 2 1 | 1 ) = 1/2*sqrt(2)
    ( 1 2 | 1 ) = 1/2*sqrt(2)

    2 A2g
    ( 2 1 | 1 ) = -1/2*sqrt(2)
    ( 1 2 | 1 ) = 1/2*sqrt(2)

    4 Eg
    ( 2 2 | 1 ) = 1
    ( 1 1 | 2 ) = 1

This is an example taken from the `examples/clebschgordan`
directory included in the code. The irreps are indexed as
`A1g->0`, `A2g->1` and `Eg->2`. The first block, for
instance, tells us that the only state belonging to $A_{1g}$
resulting from the decomposition of products of states in
$E_{g}$ is 

$$
|A_{1g},1\rangle = \frac{1}{\sqrt 2}\left(
        |E_g,2\rangle \otimes |E_g,1\rangle
        +
        |E_g,1\rangle \otimes |E_g,2\rangle
        \right).
$$

Notice that the code can parse expressions such as `sqrt(2)`
provided that they are valid as expressions in the Julia
language. More example files can be found in the
`examples/clebschgordan` directory, which contains all the
Clebsch-Gordan coefficients necessary to perform
calculations for (i) a simple system with orbitals
belonging to the identity representation (notice that
$A_{1g}$ is the identity representation of $\mathcal
O_h$, but it can be used in this case for any other identity
representation such as the irrep $s$ of the rotation group
$O(3)$), and (ii) a system with two $E_{g}$ orbitals, since
any combination of states belonging to $E_{g}$ can only
yield states belonging to $A_{1g}$, $A_{2g}$ and $E_{g}$.


# Thermodynamic calculations <a name="thermodynamics"></a>

Thermodynamic calculations are performed in every NRG
calculation due to the low computational cost involved. 
The only parameter that targets thermodynamic calculations
in particular is [`betabar`](#betabar), which defines the
temperature scale associated to each NRG step. In order to
obtain the final thermodynamic curves, results from even and
odd iterations are stored separately and then linearly 
interpolated and averaged in order to reduce oscillations
(see Rev. Mod. Phys. 80, 395).

The data obtained at each calculation is stored in a file
called `thermodata/thermo_$label_$datatype_z$z.dat`, where
the interpolated values `label` and `z` are those given as
input to `nrg_full` (see [`label`](#label) and [`z`](#z)),
and `datatype` is `clean` for `calculation="CLEAN"`,
`imp` for `calculation="IMP"` or `diff`, which
contains the impurity contribution as the subtraction
of the functions in `imp` and `clean`. In $z$-averaged
calculations, another file with `z=avg` is created,
which contains the average over files with the
specified values of `z`. Every file contains a `#`-commented
header indicating the thermodynamic quantity that
corresponds to each column. The container directory
`thermodata` is created automatically in the working
directory if it does not already exist. 

The computed thermodynamic quantities are the following,
in the order where they appear as columns in the data files.
* Temperature $T$.
* Magnetic susceptibility as $k_{B} T \chi / (g\mu_{B})$, where
$k_{B}$ is the Boltzmann constant, $T$ is the temperature,
$\chi$ is the actual magnetic susceptibility, $g$ is the
Landé factor, and $\mu_{B}$ is the Bohr magneton.
* Entropy as $S/k_{B}$, where $S$ is the actual entropy.
* Heat capacity as $C/k_{B}$, where $C$ is the actual heat
capacity.
* Free energy as $F/k_{B} T$, where $F$ is the actual free
energy.
* Average number of particles $\langle \hat N \rangle$.
* Energy as $\langle H \rangle$.
* Partition function $Z$.


# Spectral function calculations <a name="spectral"></a>

The spectral function is computed when `spectral=true`. The
calculations are performed as in J. Phys. Soc. Jpn. 58, pp.
3666-3678, and a Gaussian broadening kernel is used to
obtain smooth spectra. The parameters that affect the 
calculation of spectral functions are
[`spectral_broadening`](#broadening), [`K_factor`](#K), and 
[`orbitalresolved`](#orbitalresolved).

In orbital-resolved calculations, results for one-electron
excitations are grouped by multiplets: excitations
$f_\alpha^{(\dagger)}$ belonging to the multiplet
$m_\alpha$, 

$$
\langle F | f_\alpha^{(\dagger)} | I \rangle,
    \;\; \alpha\in m_\alpha,
$$
where $I$ and $F$ are the initial and final states, respectively,
are in the same group. Spectral funcions for each of
these groups are stored in separate files for calculations
with [individual $z$](#orbitalresolvedfile) and
[$z$-averaged](#orbitalresolvedfilezavg). Notice that the
one-electron multiplet $m_\alpha=(\Gamma_\alpha,r_\alpha)$,
where $\Gamma_\alpha=(N_\alpha=1,I_\alpha,S_\alpha=1/2)$, is
completely specified by the orbital irrep $I_\alpha$ and
the outer multiplicity $r_\alpha$, hence the name
"orbital-resolved".

The main result of the calculations are stored in files with
the following names (similar convention as in [here](#thermo)):
* `spectral/spectral_$label_z$z.dat` contains the main
result of each individual calculation, which results from
averaging over even and odd iterations.
* `spectral/spectral_$label_z$z_even.dat` and 
`spectral/spectral_$label_z$z_odd.dat` contain data
obtained from even and odd iterations, respectively. 
* `spectral/spectral_$label_z$z_splined.dat` contains the
spline interpolation of the data in
`spectral/spectral_$label_z$z.dat`.
* `spectral/spectral_$label_zavg.dat` contains $z$-averaged
values as in [thermodynamic calculations](#thermo).
* `spectral/spectral_$label_z$z_o$o.dat`<a
name="orbitalresolvedfile"></a> contains the
orbital-resolved spectral function for excitations belonging
the multiplet with index `o`. The index `o` is assigned by the
code, and the multiplet to which it corresponds is indicated
in the header line of the file. 
* `spectral/spectral_$label_zavg_o$o.dat`<a
name="orbitalresolvedfilezavg"></a> contains the
$z$-average of data from [individual $z$](#orbitalresolvedfile).

# Averaging over the twisting parameter $z$ <a name="z-averaging"></a>
In some cases, the discretization of the conduction band
leads to errors that are computationally very expensive to
eliminate by varying the parameters in a single calculation.
The problem usually revolves around the need to choose a
larger value of the discretization parameter $\Lambda$ in
order to improve the calculation in one way or another: this
forces one to impose a larger multiplet cutoff to obtain
converged calculations, which in turn increases the
computational cost of the calculation exponentially. 
The $z$-averaging procedures provides a way around this
problem by introducing a way to take several calculations
with different values of $z$ into account, which increases
the computational cost linearly with the number of
calculations.

In the calculation of thermodynamic functions, the use of a
large value of $\Lambda$ introduces spurious oscillations in
the thermodynamic functions. These can be eliminated by
averaging over several values of $z$, which makes it
possible to obtain smooth curves with a lower cutoff
requirement. See Phys. Rev. B 49, 11986.

In the calculation of spectral functions, the problem with
picking a large value of $\Lambda$ is two fold: on the one
hand, the resolution of the spectral function is very poor
for large energy ranges in the linear scale; on the other
hand, a larger broadening has to be chosen in
order to include contributions in broader energy windows.
The $z$-averaging improves the resolution because by
providing a different energy grid for each calculation. On
top of that, it works with broadening parameters of the
order of $\eta\approx 1/N_{z}$, where $N_{z}$ is the number of
values of $z$, because the energy windows from different $z$
"patch" together to cover the whole spectrum.


## Parallelization<a name="parallelization"></a>

The $z$-averaging procedure is where parallelization can be
exploited most optimally. Since calculations for each $z$ are
completely independent from each other, the linear increase
in computational time introduced by the method can be
greatly reduced, if not almost eliminated, by running
calculations for different $z$ in parallel.

Parallelization is implement using the tools in the
`Distributed` Julia package, which is documented in the
[Julia docs](https://docs.julialang.org/en/v1/manual/distributed-computing/).
This package has to be loaded including the `using
Distributed` line and adding the processors with the
`addprocs` function. All the modules and variables
have to be loaded in the various processes using the
`@everywhere` macro. Parallel processes can then spanned
using a `@distributed` for loop, with the added `@sync` macro 
to ensure that the code waits until every process
is completed.

For thermodynamic calculations, $z$-averaging can be simply
achieved by first running calculations for various $z$ without the
impurity (`calculation="CLEAN"`) and then calculations
with the impurity (`calculation="IMP"`). The code has the
following structure:

    Z = generate_Z( Nz )

    addprocs(number_of_processes)

    @everywhere begin
        # load all modules
    end

    for calculation in ["CLEAN","IMP"]

        @everywhere begin
            # load all variables
        end

        @sync @distributed for z in Z
            nrg_full(...)
        end

    end
    zavg_thermo( label , Z )

    rmprocs(number_of_processes)

The variable `Z` is a vector containing `Nz` values of $z$
and is constructed in this example by the function
`generate_Z`, included in the module `modules/zavg.jl`.
The `@distributed` macro instructs the program to run the
various `nrg_full` calculations in parallel, and the `@sync`
macro ensures that the program waits until every calculation
is completed before proceding. Finally, the `zavg_thermo`
function, also part of the `modules/zavg.jl` module,
averages over thermodynamic calculations for every $z$
value in `Z` for the system labeled as
[`label`](#label).

For spectral function calculations, the structure is very
similar, the only difference being that the variable
`calculation="IMP"` does not change and that the averaging
function is `zavg_spectral`:

    @everywhere calculation = "IMP"

    @sync @distributed for z in Z
        nrg_full(...)
    end

    zavg_zavg( 
        label , 
        Z ,
        orbitalresolved_number=number_of_orbitals 
    )

The `zavg_spectral` function takes the additional optional
argument `orbitalresolved_number::Float64=0`. Its default
value `0` is used when the calculation is not
orbital-resolved. For orbital-resolved calculations, the
number of orbitals `number_of_orbitals` has to be provided.

# Template scripts and workflow<a name="templates"></a>

Two template scripts are provided in the `templates`
directory. These are designed to be used as a starting point
for any calculation. A typical workflow would start by
creating a directory for a given system, copying the
template scripts and using them to perform the steps
described [above](#functions).

The `multiplets_TEMPLATE.jl` performs the multiplet
calculation ([step 1](#step_1)). This has be run once for
every orbital irrep.

The `nrg_TEMPLATE.jl` performs the two-particle multiplet
calculation ([step 2](#step_2)), the impurity spectrum
calculation ([step 3](#step_3)), and the NRG calculation 
([step 4](#step_4)). These are grouped together in the same
script because (i) they share many common parameters and
(ii) a typical workflow can involve back-and-forth running
of the various steps included.

# Examples and tutorial<a name="examples"></a> 

A directory named `examples` is provided with the code. The
aim of these examples is to provide guidance and ideas as to
how to structure calculations using the code. It
contains two subdirectories, each containg scripts and data
files for calculations performed for a different system:

* `A1g`: Calculations for a system with one orbital
transforming as the $A_{1g}$ irrep of the cubic group
$O_{h}$ and one channel of the same symmetry coupled to that
orbital. The irrep $A_{1g}$ is the identity irrep, so this
system is equivalent to the single-impurity Anderson model
first studied by Wilson and co-workers in Phys. Rev. B 21, 1003.

* `Eg`: Calculations for a system with two orbitals
transforming as the $E_{g}$ irrep of the cubic group $O_{h}$
and two channels of the same symmetry, each coupled to its
corresponding orbital. 

In both examples we find modified versions of the
[template scripts](#templates). The Clebsch-Gordan
coefficient directory is the same for both: the one shipped
with the code. 

A tutorial on how to run calculations following the same
scheme proposed in the examples is explained in
`doc/Tutorial.md`, which can be used as additional guidance
for performing other calculations.
