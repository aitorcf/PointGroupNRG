# Table of contents 
[Overview](#overview)
1. [Multiplet calculation](#multiplets)
2. [Two-particle multiplet calculation](#2mults)
3. [Impurity spectrum calculation](#spectrum)
4. [Full NRG calculation](#nrg)
    4.1. [Thermodynamic calculation](#thermo)
    4.2. [Spectral function calculation](#spectral)

# Overview<a name="overview"></a>

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

This document is a guide for carrying out each of the four
steps. Detailed information about the functions and their
arguments can be found in `doc/Manual.md`.

# 1. Multiplet calculation<a name="multiplets"></a>

To start with, create a working directory for the
calculations. Copy the file `multiplets_TEMPLATE.jl` to that
directory. The file is a template for running multiplet
calculations for a single orbital irrep, so it has to be run
once for every orbital irrep needed. The main function,
called at the end of the script, is `compute_multiplets`.
Run the script for each one-electron orbital irrep. The
results are stored in the specified directory.

# 2. Two-particle multiplet calculation <a name="2mults"></a>

Copy the template file `nrg_TEMPLATE.jl` to the working
directory. The file contains contains the functions
`multiplets_2particles`, `impurity_spectrum` and `nrg_full`
corresponding to steps 2, 3 and 4, respectively. In order to
compute the multiplets, we set `run="multiplets"` so that
the code runs the function `multiplets_2particles` and 
we define every input parameter required by the function.
We then run the script and check its output, which prints
the multiplets in the format

    (N,I,S,r)

where `N` is the number of particles, `I` is the orbital
irrep, `S` is the total spin, and `r` is the outer
multiplicity number.

# 3. Impurity spectrum calculation<a name="spectrum"></a>

We compute the atomic spectrum by setting `run="spectrum"`
in the `nrg_TEMPLATE.jl` file and defining the input
parameters that have not yet been defined:
`epsilon_symparams` and `u_symparams`. Then we run the
script and we read the spectrum in the output. The spectrum
is formatted in two columns: the first column specifies the
multiplet (in the format described in the [previous
section](#2mults) and the second column specifies the energy:

# 4. Full NRG calculation<a name="nrg"></a>

The last step is to perform a full NRG calculation, for
which we set `run="nrg"`, which can be used to compute
thermodynamic or spectral functions.

## 4.1. Thermodynamic calculation<a name="thermo"></a>

In order to perform a thermodynamic calculation, we have to
define only the mandatory input parameters of the function
`nrg_full` that will be used in the calculation. We have to
run two calculations: one for the conduction band only, for
which we set `calculation="CLEAN"`, and one for the system
with the impurity, for which we set `calculation="IMP"`. An
easy way to achieve this is to loop over them instead of
running the script twice (this also saves compilation time):

    for calculation in ["CLEAN","IMP"]
        nrg_full( ... , label , ... )
    end

It is important to run the calculations in that specific
order: first `calculation="CLEAN"`, then
`calculation="IMP"`. Otherwise, the impurity contribution
will not be calculated.

Once the calculations are completed, the output files are stored
in the `thermodata` directory. The impurity contribution
to the thermodynamic functions is found in the file with the
`diff` label in its name.

## 4.2. Spectral function calculation<a name="spectral"></a>

In order to perform a spectral function calculation, we have
to proceed as for the thermodynamic calculation, but with
two main differences.
1. The only calculation we have to perform in the one for
   the system with the impurity, so we set
   `calculation="IMP"`.
2. We have to set `spectral=true` and, optionally, specify
   some of the optional parameters: `spectral_broadening`,
   `K_factor`, `orbitalresolved`. If they are not defined,
   the calculation will be performed with the default
   values.

