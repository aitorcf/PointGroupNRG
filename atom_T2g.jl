#!/usr/bin/env julia

#   ======================   #
#%% CALCULATION PARAMETERS %%#
#   ======================   #   

orbital = "T2g"

# construction and diagonalization in a single step
onestep = true
# use precomputed spin-CG coefficients
spinarray = true
# parallel computation
distributed = parallel = false
# parallel method: distfor or async
method = "distfor"
# discretization ("standard" or "co2005")
discretization = "co2005"

# twisting parameter
z = 0.0

calculation="IMP"

# numerical parameters
L = 10.0
betabar = 1.0

iterations = 20

# cutoff
cutoff_type = "multiplet" 
cutoff_magnitude = 100
minmult = 200

# one-particle parameters
eps  = -0.1 
gam  =  0.01
# standard coulomb parameters 
U_11 = 0.4
U_12 = 0.3 
J    = 0.3
if length(ARGS)==3
    U_11 = parse( Float64 , ARGS[1] )
    U_12 = parse( Float64 , ARGS[2] )
    J    = parse( Float64 , ARGS[3] )
elseif length(ARGS)==4 
    # U_11 is given but not read
    U_12 = parse( Float64 , ARGS[2] )
    J    = parse( Float64 , ARGS[3] ) 
    D    = parse( Float64 , ARGS[4] )
    U_11 = U_12 + 2*J + D
end

# relative coulomb parameters 
u_11 = 1.0
u_12 = 0.1
j    = 0.0
U_11 = u_11*2*abs(eps)
U_12 = u_12*U_11
J    = j*U_12
# multiplet coulomb parameters
u_a1g = U_11 + 2*J # on-site, orbital-symmetric
u_eg  = U_11 - J   # on-site, orbital-symmetric
u_t1g = U_12 - J # inter-site, spin-aligned
u_t2g = U_12 + J # inter-site, spin-antialigned

max_spin2 = 10;

distworkers  = 6

println( "SETUP AND PARAMETERS" )
println( "====================" )
@show calculation
@show distributed 
distributed && @show distworkers 
@show method
@show eps   
@show u_a1g
@show u_eg
@show u_t1g
@show u_t2g
@show gam  
@show L
@show betabar
@show cutoff_type
@show cutoff_magnitude
@show iterations
@show max_spin2;



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
include( "modules/symmetry.jl" );

if parallel 

    using Distributed 

    # remove existing workers
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

ORBITAL_IRREPS = ( orbital , )
oirreps = cg_shortcircuit( CG_PATH , ORBITAL_IRREPS... )
sort!( oirreps )

#%%
oirreps2indices = Dict( o=>i for (i,o) in enumerate(oirreps) )
oirreps2dimensions = Dict( "Eg" => 2 ,
                           "A1g"=> 1 ,
                           "A2g"=> 1 , 
                           "T1g"=> 3 , 
                           "T2g"=> 3 )
oindex2dimensions = collect( oirreps2dimensions[I] 
                             for I in oirreps )
println( "ORBITAL IRREPS FOR THE PROBLEM" )
println( oirreps )
println()

println( "Computing orbital CG coefficients..." )
@time begin
    global cg_o_full = get_cg_o_fulldict( oirreps , CG_PATH )
    global cg_o_fullmatint = get_cg_o_fullmatint( cg_o_full , oirreps )
end
println()
println( "CLEBSCH-GORDAN MATRIX" )
print_dict( cg_o_fullmatint )
println()


#   ===================   #
#%% SPIN CLEBSCH-GORDAN %%#
#   ===================   #
println( "Computing spin CG coefficients..." )
@time cg_s_fullmatint = get_cg_s_fullmatint( max_spin2 );
println()


#   ==================================   #
#&& ONE-SHELL SYMSTATES AND MULTIPLETS &&#
#   ==================================   #

#   ----------------------- #
#%% basis and hilbert space #
#   ----------------------- #
hilbert_0 = HilbertSpace( T2g_states(0) )
println( "Computing atomic basis..." )
@time basis_0 = CanonicalBasis( hilbert_0 )
println()
println( "BASIS FOR ATOMIC SHELL" )
println( basis_0 )
println()


#   ------------------------ #
#%% symstates and multiplets #
#   ------------------------ #
hiztegia = Dict( 
    "t" =>  "T2g",
    "u" =>  0.5,
    "d" => -0.5
)

print( "Computing atomic symstates..." )
@time begin
    global symstates_0_nor = oneirrep_symstates( 
                                hilbert_0 , 
                                hiztegia , 
                                "A1g" , 
                                "$(ASYM_PATH)$(orbital)_julia/" )
    global symstates_0 = Dict( (q[1:5]...,1)=>s for (q,s) in symstates_0_nor )
end
println( " done" )
println()
println( "ATOMIC SYMSTATES" )
for (k,s) in symstates_0 
    k[1]!==2 && continue
    @show k 
    println( s )
    println()
end
println()

println( "ATOMIC MULTIPLETS" )
multiplets_0 = get_multiplets( symstates_0 )
nmultiplets = collect(multiplets_0) 
sort!( nmultiplets , by=x->x[1] )
for m in nmultiplets
    @show m 
end
println()


#   ---------   #
#%% operators %%#
#   ---------   #
epsilon_symparams = Dict( 
    orbital => eps
)
print( "Computing ϵ operator..." )
@time epsilon = epsilon_sym( symstates_0 , epsilon_symparams ; verbose=false )
println( "done" )
#println()
#println( "EPSILON OPERATOR" )
#print_diagonals( epsilon )

u_symparams = Dict( 
    ("A1g", 0) => [u_a1g][:,:],
    ("Eg",  0) => [u_eg ][:,:],
    ("T1g", 1) => [u_t1g][:,:],
    ("T2g", 0) => [u_t2g][:,:]
)
print( "Computing Coulomb operator..." )
@time coulomb = u_sym( symstates_0 , u_symparams ; verbose=false )
println( " done")
println()
#println( "COULOMB OPERATOR" )
#print_diagonals( coulomb )
#println()


#   -----------------------------   #
#%% atomic hamiltonian (rescaled) %%#
#   -----------------------------   #
H = epsilon + coulomb 
# α = ϵ_0^z (scaling factor)
if discretization=="standard"
    global α = 0.5 * L^z * (1+L^-1)
elseif discretization=="co2005"
    global α = compute_ebar0_z(z,L;discretization=discretization)
end
H.matrix ./= α
gam = sqrt(2.0*gam/pi) / α # for later
#println( "H OPERATOR" )
#print_diagonals( H )
#println()
#hermitize!( H )
#println( "Diagonalizing atomic Hamiltonian without symmetries..." )
#@time (e,u) = diagonalize( H )
#println()
#eigenbasis = GenericBasis( basis_0 , u )
#println( "EIGENENERGIES" )
#println( e )
#println()


#   ================================================   #
#&& NUMERICAL ADDITION OF INNERMOST CONDUCTION SHELL &&#
#   ================================================   #

#   -----   #
#%% irreu %%#
#   -----   #
# imp
irreps_0 = get_irreps( symstates_0 )
println( "Diagonalizing atomic Hamiltonian with symmetries..." )
@time irrEU_imp = symdiag( irreps_0 , symstates_0 , H )
println()
println( "irrEU for IMP" )
print_dict( irrEU_imp )
println()
# clean
irrEU_clean = get_irrEU_clean( "A1g" )
println( "irrEU for CLEAN" )
print_dict( irrEU_clean )
println()

irrEU = irrEU_imp

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
for l in sort(levels,by=x->x[1][2])
    m = l[1] 
    e = l[2] 
    println("$m => $e")
end

writedlm( "atomic_levels_T2g.dat" , levels )
writedlm( "phasediag_T2g_brokensym/u11_$(U_11)_u12_$(U_12)_j_$(J).dat" , levels )


