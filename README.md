# Overview 

This package is designed to work with magnetic impurity
models with symmetry 

$$
\mathbb G = U(1)_c\otimes P_o \otimes SU(2)_s,
$$

where $U(1)_c$ is the charge symmetry (particle
conservation), $P_o$ is the finite (orbital) point group
symmetry for the orbital, $SU(2)_s$ is the full rotational
symmetry in spin space. We can define arbitrary models for
any given simply reducible point group $P_o$.

A calculation consists of two main steps 

1. Calculation of multiplet symmetry-adapted states. This is
   perform once for each system.
   
2. NRG calculation.

# Calculation of multiplet states 





# Definition of a model

A model is defined by specifying how many orbital 
multiplets we have in the impurity and in the shells.
For example, a two-orbital $E_g$ impurity where each orbital
is connected to its corresponding channel is defined as
follows:

    # Two-orbital Eg model
    atom_config  = Dict( "Eg" => 1 ) 
    shell_config = Dict( "Eg" => 1 )

By this we mean that we have an atomic (impurity) part with
one group of $\text{dim}(E_g)=2$ orbitals belonging to the irrep
$E_g$, and the same for the shell. We can also define a
two-impurity, one-channel system by including two $s$ orbitals
coupled to a single $s$ channel:

    # Two-impurity, one-channel model 
    atom_config  = Dict( "A1g" => 2 ) 
    shell_config = Dict( "A1g" => 1 )

We have used the irrep $A_{1g}$ as an equivalent to $s$
because both are identity representations.

# Anderson Hamiltonian 

The Anderson Hamiltonian for the system is defined by the
occupation energies $\epsilon_{r}(\Gamma)$, the Coulomb
parameters $U_{rr'}(\Gamma)$, and the hybridization
parameters $V_{r}(\Gamma)$. For instance, we can define an 
$E_g$ system as follows:

    # Occupation energies 
    epsilon_symparams = Dict( 
        ("Eg",1) => ComplexF64(eps)
    )
    # Coulomb parameters
    u_symparams = Dict( 
        ("A1g",0) => ComplexF64[u_a1g][:,:],
        ("A2g",1) => ComplexF64[u_a2g][:,:],
        ("Eg", 0) => ComplexF64[u_eg][:,:]
    )
    # Hybridization parameters
    hop_symparams = Dict( 
        ("Eg") => sqrt(2*gam/pi)*ComplexF64[1][:,:]
    )


# Calculation parameters
