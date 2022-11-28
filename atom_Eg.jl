#!/usr/bin/env julia

#   ======================   #
#%% CALCULATION PARAMETERS %%#
#   ======================   #   

orbital = "Eg"
norbitals = 2

# construction and diagonalization in a single step
onestep = true
# use precomputed spin-CG coefficients
spinarray = true
# parallel computation
distributed = parallel = true
# parallel method: distfor or async
method = "distfor"
# discretization ("standard" or "co2005")
discretization = "co2005"

# clean system or with impurity added
calculation = "IMP"

# twisting parameter
z = 0.0

# numerical parameters
L = 10.0
betabar = 1.0

# one-body coupling parameters
eps  = -0.1
# relative coulomb parameters
u_11_r = 0.9
u_h_r  = 0.2

## standard parameters 
#U   = 0.2
#U12 = 0.8*U
#J   = 0.5*( U - U12 )
#u_eg =  U - J
#u_a1 =  U + J
#u_a2 =  U12 - J 

# coulomb coupling parameters
u_11 = u_11_r*abs(2*eps)
u_12 = u_11 
u_h  = u_h_r*( 2*u_12 )
u_eg = u_11
u_a1 = u_12 + u_h/2.0
u_a2 = u_12 - u_h/2.0

iterations = 20

max_spin2 = 8

distworkers = 6

# ---------- #
# PARAMETERS #
# ---------- #
println("ATOMIC PARAMETERS") 
println("standard parameters")
@show eps
@show u_11 
@show u_12 
@show u_h
println("multiplet parameters")
@show u_eg 
@show u_a1 
@show u_a2 
println()


#   =======================   #
#%% MODULES AND DISTRIBUTED %%#
#   =======================   #

using DelimitedFiles
using ProfileVega
using Profile

include( "modules/symbols.jl" )
include( "modules/numericals.jl" )
include( "modules/compoundoperators.jl" )
include( "modules/shell.jl" )
include( "modules/multithread.jl" )
include( "modules/thermo.jl" )
include( "modules/discretization.jl" )

if parallel 

    using Distributed 

    # kill current processes
    for i in workers()
        t = rmprocs(i, waitfor=0)
        wait(t)
    end

    # add requested workers
    if distworkers ≥ nprocs()
        addprocs(distworkers)
    else 
        println( "more workers than processors!" )
    end

    println( "DISTRIBUTED CALCULATION WITH $(nworkers()) WORKERS" )

    @everywhere begin 
        using ProgressMeter
        using PartialWaveFunctions
        using StaticArrays
        include( "modules/symmetry.jl" )
        include( "modules/distributed.jl" )
    end

else 

    println( "SERIAL CALCULATION" )

    using ProgressMeter
    using PartialWaveFunctions
    using StaticArrays

    include( "modules/symmetry.jl" )
    include( "modules/distributed.jl" )

end
println()


#   ======================   #
#%% ORBITAL CLEBSCH-GORDAN %%#
#   ======================   #

CG_PATH = "/home/acalvo/Bulegoa/ClebschGordan/Oh/cg_symbolic/"
ASYM_PATH = "/home/acalvo/Bulegoa/AntiSymmetricPart/Oh/";

ORBITAL_IRREPS = ( "Eg" , )
oirreps = cg_shortcircuit( CG_PATH , ORBITAL_IRREPS... )
oirreps2indices = Dict( o=>i for (i,o) in enumerate(oirreps) )
oirreps2dimensions = Dict( "Eg" => 2 ,
                           "A1g"=> 1 ,
                           "A2g"=> 1 )
oindex2dimensions = collect( oirreps2dimensions[I] for I in oirreps )
println( "ORBITAL IRREPS FOR THE PROBLEM" )
println( oirreps )
println()

cg_o_full = get_cg_o_fulldict( oirreps , CG_PATH )
cg_o_fullmatint = get_cg_o_fullmatint( cg_o_full , oirreps )
#println( "CLEBSCH-GORDAN MATRIX" )
#print_dict( cg_o_fullmatint )
#println()


#   ===================   #
#%% SPIN CLEBSCH-GORDAN %%#
#   ===================   #
cg_s_fullmatint = get_cg_s_fullmatint( max_spin2 );


#   ==================================   #
#&& ONE-SHELL SYMSTATES AND MULTIPLETS &&#
#   ==================================   #

#   ----------------------- #
#%% basis and hilbert space #
#   ----------------------- #
hilbert_0 = HilbertSpace( Eg_states(0) )
basis_0 = CanonicalBasis( hilbert_0 )
println( "BASIS FOR ATOMIC SHELL" )
println( basis_0 )


#   ------------------------ #
#%% symstates and multiplets #
#   ------------------------ #
hiztegia = Dict( 
    "e" => "Eg",
    "u" => 0.5,
    "d" =>-0.5
)

symstates_0_nor = oneirrep_symstates( hilbert_0 , 
                                      hiztegia , 
                                      "A1g" , 
                                      "$(ASYM_PATH)$(orbital)_julia/" )
symstates_0 = Dict( (q[1:5]...,1)=>s for (q,s) in symstates_0_nor )
println( "SYMSTATES" )
for q in sort(collect(keys(symstates_0)),by=x->x[1])
    @show q 
    println( symstates_0[q] )
    println()
end
println()

multiplets_0 = get_multiplets( symstates_0 )
println( "ATOMIC MULTIPLETS" )
for m in multiplets_0 
    @show m 
end
println()


#   ---------   #
#%% operators %%#
#   ---------   #
epsilon_symparams = Dict( 
    "Eg" => eps
)
epsilon = epsilon_sym( symstates_0 , epsilon_symparams ; verbose=false )
#println( "EPSILON OPERATOR" )
#print_diagonals( epsilon )

u_symparams = Dict( 
    ("Eg", 0) => [u_eg][:,:],
    ("A1g",0) => [u_a1][:,:],
    ("A2g",1) => [u_a2][:,:]
)
coulomb = u_sym( symstates_0 , u_symparams ; verbose=false )
println( "COULOMB OPERATOR" )
print_diagonals( coulomb )
println()


#   ---------------------------   #
#%% atomic hamiltonian (scaled) %%#
#   ---------------------------   #
H = epsilon + coulomb 
# α = ϵ_0^z (scaling factor)
if discretization=="standard"
    global α = 0.5 * L^z * (1+L^-1)
elseif discretization=="co2005"
    global α = compute_ebar0_z(z,L;discretization=discretization)
end
H.matrix ./= α
#println( "H OPERATOR" )
#print_diagonals( H )
#println()

#   -----   #
#%% irreu %%#
#   -----   #
# imp
if calculation=="IMP"
    irreps_0 = get_irreps( symstates_0 )
    irrEU = symdiag( irreps_0 , symstates_0 , H )
    println( "irrEU for IMP" )
    print_dict( irrEU )
    println()
end

levels = []
global minE = 10
for (G,(E,U)) in irrEU 
    for r in 1:length(E) 
        m = (G...,r)
        if E[r]<minE
            global minE = E[r]
        end
        push!( levels , [m,E[r]] )
    end
end
levels = [[l[1],l[2]-minE] for l in levels]
levels = map( x->[replace("$(x[1])", " "=>"" ),x[2]], levels )
levels = map( x->[replace("$(x[1])", "\""=>"" ),x[2]], levels )
println( "LEVELS" )
for l in sort(levels,by=x->x[2])
    m = l[1] 
    e = l[2] 
    println("$m => $e")
end

writedlm( "atomic_levels_Eg.dat" , levels )

