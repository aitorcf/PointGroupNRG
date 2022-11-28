#!/usr/bin/env julia

#   ==========   #
#%% PARAMETERS %%#
#   ==========   #

orbital = "A1"
norbitals = 1 

identityrep = "A1"

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
# calculation 
calculation = "IMP"

z = 0.0
if length(ARGS)>0
    # clean/imp
    calculation = "IMP"
    if ARGS[1]=="CLEAN" 
        calculation = "CLEAN"
    end
    # twisting parameter
    if length(ARGS)>1
        z = parse(Float64,ARGS[2])
    end
end

# Free orbital --> Local moment ==> Strong coupling
eps = -0.1
cou =  0.15
gam =  0.01
#gam = 0.0

L = 3.0
betabar = 1.0

cutoff_type = "multiplet" # multiplet/energy
cutoff_magnitude = 50

iterations = 3

spinarray = true  # faster
max_spin2 = 6

distworkers = 6

println( "================" )
println( "SETUP PARAMETERS" )
println( "================" )
@show onestep 
@show distributed 
@show calculation 
@show method 
@show z 
@show eps 
@show cou 
@show gam 
@show L 
@show betabar 
@show cutoff_type 
@show cutoff_magnitude 
@show iterations 
@show spinarray 
@show max_spin2 
@show distworkers


#   =======   #
#%% MODULES %%#
#   =======   #

using Distributed 


using BenchmarkTools 
using DelimitedFiles
using ProfileVega
using Profile

include( "modules/symbols.jl" )
include( "modules/numericals.jl" )
include( "modules/compoundoperators.jl" )
include( "modules/shell.jl" )
include( "modules/thermo.jl" )
include( "modules//spectral.jl" )

if distributed 

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

    @everywhere begin 
        using ProgressMeter
        using PartialWaveFunctions
        using StaticArrays
        include( "modules/symmetry.jl" )
        include( "modules/diagonalization.jl" )
    end

    println( "DISTRIBUTED CALCULATION WITH $(nworkers()) WORKERS" )

else

    using ProgressMeter
    using PartialWaveFunctions
    using StaticArrays
    include( "modules/symmetry.jl" ) 
    include( "modules/diagonalization.jl" )

    println( "SERIAL CALCULATION" )

end
println()


#   ======================   #
#%% ORBITAL CLEBSCH-GORDAN %%#
#   ======================   #

CG_PATH = "/home/aitor/Bulegoa/ClebschGordan/C4v/cg_symbolic/"
ASYM_PATH = "/home/aitor/Bulegoa/AntiSymmetricPart/C4v/";

ORBITAL_IRREPS = ( "A1" ,)

oirreps = cg_shortcircuit( CG_PATH , ORBITAL_IRREPS... )
oirreps2indices = Dict( o=>i for (i,o) in enumerate(oirreps) )
oirreps2dimensions = Dict( "A1" => 1 ) 
oindex2dimensions = collect( oirreps2dimensions[I] for I in oirreps )

cg_o_full = get_cg_o_fulldict( oirreps , CG_PATH )
cg_o_fullmat = get_cg_o_fullmat( cg_o_full )
cg_o_fullmatint = get_cg_o_fullmatint( cg_o_full , oirreps )


#   ===================   #
#%% SPIN CLEBSCH-GORDAN %%#
#   ===================   #
cg_s_fullmatint = get_cg_s_fullmatint( max_spin2 );

#   ===========================   #
#&& SYMBOLIC CALCULATION OF N=0 &&#
#   ===========================   #

#   ------------ #
#%% state tuples # 
#   ------------ #
a_1 = ( 0 , "a" , 1 )
a_1u = ( a_1... , "u" )
a_1d = ( a_1... , "d" )

tuples_0 = [ a_1u , a_1d ]

#   ----- #
#%% basis #
#   ----- #
hilbert_0 = HilbertSpace( tuples_0... )
basis_0 = CanonicalBasis( hilbert_0 )
println( "HILBERT SPACE" )
println( hilbert_0 )
println()
println( "BASIS" )
println( basis_0 )


#   --------- #
#%% symstates #
#   --------- #
hiztegia = Dict( 
    "a" =>  "A1",
    "u" =>  0.5 ,
    "d" => -0.5
)

symstates_0 = oneirrep_symstates( hilbert_0 , hiztegia , "A1" , "$(ASYM_PATH)A1/" )
symstates_0_new = Dict()
for (q,s) in symstates_0 
    push!( symstates_0_new , (q[1:5]...,1)=>s ) 
end
symstates_0 = symstates_0_new
println( "SYMSTATES" )
print_dict( symstates_0 )

multiplets_0 = get_multiplets( symstates_0 )

irrepsubs_0 = get_multiplets( symstates_0 )
println( "MULTIPLETS" )
for m in irrepsubs_0 
    println( m )
end
println()

#   --------- #
#%% operators #
#   --------- #
epsilon_symparams = Dict( 
    "A1" => eps
)
epsilon = epsilon_sym( symstates_0 , epsilon_symparams ; verbose=false )
println( "EPSILON OPERATOR" )
print_diagonals( epsilon )

u_symparams = Dict( 
    ("A1",0) => [cou][:,:]
)
coulomb = u_sym( symstates_0 , u_symparams ; verbose=false );
println( "COULOMB OPERATOR" )
print_diagonals( coulomb )


#   ---------------- #
#%% full hamiltonian #
#   ---------------- #
H = epsilon + coulomb 
# α = ϵ_0^z (scaling factor)
if discretization=="standard"
    global α = 0.5 * L^z * (1+L^-1)
elseif discretization=="co2005"
    global α = compute_ebar0_z(z,L;discretization=discretization)
end
H.matrix ./= α
gam = sqrt(2.0*gam/pi) / α # for later
#gam /= α


#   ================================================   #
#&& NUMERICAL ADDITION OF INNERMOST CONDUCTION SHELL &&#
#   ================================================   #

#   -----   #
#%% irreu %%#
#   -----   #
# imp
if calculation=="IMP"
    irreps_0 = get_irreps( symstates_0 )
    irrEU_imp = symdiag( irreps_0 , symstates_0 , H )
    normalize_irrEU( irrEU_imp )
    println( "irrEU for IMP" )
    print_dict( irrEU_imp )
    println()
# clean
elseif calculation=="CLEAN"
    irrEU_clean = get_irrEU_clean( "A1g" )
    normalize_irrEU( irrEU_clean )
    println( "irrEU for CLEAN" )
    print_dict( irrEU_clean )
    println()
end


#   ---------------   #
#%% combinations u' %%#
#   ---------------   #
combinations_uprima = 
    Dict{ Tuple{Int64,Int64,Int64,Int64} , NTuple{2,Tuple{Int64,Int64,Int64,Int64}} }()
m_vac = (0,oirreps2indices["A1"],0,1)
#print( "COMBINATIONS U' " )
if calculation=="IMP"
    #println( "for IMP" )
    for m_mu in multiplets_0
        mint_mu = convert_to_int( m_mu , oirreps2indices )
        push!( combinations_uprima , mint_mu=>(mint_mu,m_vac) )
    end 
elseif calculation=="CLEAN" 
    #println( "for CLEAN" )
    push!( combinations_uprima , m_vac=>(m_vac,m_vac) )
end

irreps_uprima = Set( k[1:3] for k in keys(combinations_uprima) )
combinations_uprima = 
        Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} }(
            G => NTuple{3,NTuple{4,Int64}}[
                (m_u,m_mu,m_i)
                for (m_u,(m_mu,m_i)) in combinations_uprima 
                if m_u[1:3]==G
                ]
            for G in irreps_uprima
        ) 
println( "COMBINATIONS U' FOR N=0" )
for (G,combs) in combinations_uprima
    println( "$G => $combs" )
end
println()

#   ----------------------   #
#%% pseudo-CG coefficients %%#
#   ----------------------   #
pcg = get_pseudoCG( symstates_0 , 
                    basis_0 , 
                    hiztegia , 
                    oirreps2indices )
println( "PCG" )
print_dict( pcg ) 
println()

# PCG matrix in int format
irrmult_0 = get_irreps( multiplets_0 ; multiplicity=true )
pcgmat = get_pseudoCG_mat( pcg , 
                           irrmult_0 , 
                           oindex2dimensions , 
                           oirreps2indices );

#   -------------------   #
#%% shell cops and qq_a %%#
#   -------------------   #
shell_cops = shell_coperators( basis_0 , hiztegia )
qq_a = collect( convert_to_int(q_a,oirreps2indices) 
                for q_a in keys(shell_cops) )
println( "ONE-PARTICLE HOPPERS" )
@show qq_a;



#   -----------------   #
#%% hopping parameter %%#
#   -----------------   #
hop = Dict( 
    (oirreps2indices["A1"],1) => ComplexF64(gam) 
)
println( "ATOM-SHELL HOPPING PARAMETER" )
@show hop
println()
println()

#   --------------------------   #
#%% block and shell multiplets %%#
#   --------------------------   #
if calculation=="IMP"
    symstates_block = symstates_0 
    multiplets_block = get_multiplets( symstates_block )
    global irrEU = Dict( (convert_to_int(G,oirreps2indices),(E,U))
                   for (G,(E,U)) in irrEU_imp )
    println( "irrEU before adding first shell (to IMP)" )
    for (G,(E,U)) in irrEU 
        @show G, E
    end
    println()
elseif calculation=="CLEAN"
    multiplets_block = Set([(0,"A1",0.0,1)])
    global irrEU = Dict( (convert_to_int(G,oirreps2indices),(E,U))
                   for (G,(E,U)) in irrEU_clean )
    println( "irrEU before adding first shell (to CLEAN)" )
end 
symstates_shell  = symstates_0 
multiplets_shell = multiplets_0

multiplets_block = Set( convert_to_int(m,oirreps2indices) 
                        for m in multiplets_block )
multiplets_shell = Set( convert_to_int(m,oirreps2indices)
                        for m in multiplets_shell )

#   -----------------   #
#%% excitation matrix %%#
#   -----------------   #
#M = get_M0( pcg , qq_a )
#println( "===================" )
#println( "M excitation matrix" )
#println( "===================" )
#print_dict( M )
#println()
#println( "M size: $(length(collect(keys(M[qq_a[1]]))))" )
#AA = [get_A(M,irrEU)]

multiplets_a = collect(filter( x->x[1]==1 , multiplets_block ))
Mred, AA = setup_redmat_AA(
            pcg ,
            multiplets_block ,
            multiplets_a ,
            cg_o_fullmatint ,
            cg_s_fullmatint ,
            irrEU ;
            verbose=false )


#   ------------------------   #
#%% impurity quantum numbers %%#
#   ------------------------   #
println("IMPURITY MULTIPLET SPACE")
omults = ordered_multiplets(multiplets_block)
mult2index = Dict( m=>i for (i,m) in 
                   enumerate(omults))
mm_i = Dict( 
            m=>[(i==mult2index[m] ? 1.0 : 0.0)
                for i in 1:length(multiplets_block)] 
                for m in omults
           )
@show mm_i 
m_imp = mult_thermo( irrEU ,
                     betabar ,
                     oindex2dimensions ,
                     mm_i )
println()

#   ---------------------------------------   #
#%% matrix construction and diagonalization %%#
#   ---------------------------------------   #
if distributed && onestep
    @time begin
    global (irrEU,combinations_uprima) = matdiag_distributed_mat( 
                multiplets_block , 
                multiplets_shell ,
                irrEU , 
                hop , 
                cg_o_fullmatint , 
                cg_s_fullmatint ,
                pcg , 
                pcgmat, 
                qq_a , 
                combinations_uprima , 
                oindex2dimensions ;
                verbose=false );
    end
else
    @time begin
    global (irrEU,combinations_uprima) = matdiag_serial_mat( 
                multiplets_block , 
                multiplets_shell ,
                irrEU , 
                hop , 
                cg_o_fullmatint , 
                cg_s_fullmatint ,
                pcg , 
                pcgmat, 
                qq_a , 
                combinations_uprima , 
                oindex2dimensions ;
                verbose=false );
    end
end
minE = minimum([e for (E,U) in values(irrEU) for e in E])
irrEU = Dict( G=>(E.-minE,U) for (G,(E,U)) in irrEU )
energies = sort([e for (E,U) in values(irrEU) for e in E])
@show energies
@show length(energies)

mm_i = imp_mults( irrEU ,
                  oindex2dimensions ,
                  combinations_uprima ,
                  mm_i )
@show mm_i
m_imp = mult_thermo( irrEU ,
                     betabar ,
                     oindex2dimensions ,
                     mm_i )


#   ~~~~~~~~~~~~~~~~~~~   #
#%% transformation of M %%#
#   ~~~~~~~~~~~~~~~~~~~   #
@show combinations_uprima

#println( "===================" )
#println( "TRANSFORMATION OF M" )
#println( "===================" )
#println()
#q_a = qq_a[1]
#M = get_new_M( M , 
#               qq_a , 
#               irrEU , 
#               oindex2dimensions , 
#               combinations_uprima ,
#               cg_o_fullmatint , 
#               cg_s_fullmatint ;
#               verbose=true )
#println( "M matrix" )
#print_dict( M )
#println()
#println( "M size: $(length(collect(keys(M[qq_a[1]]))))" )
#
#push!( AA , get_A(M,irrEU) )

function update_redmat_AA(
            Mred ,
            irrEU ,
            combinations_uprima ,
            multiplets_a ,
            cg_o_fullmatint ,
            cg_s_fullmatint ,
            AA )

    Mred = get_new_blockredmat( 
                Mred , 
                irrEU ,
                combinations_uprima ,
                multiplets_a ,
                cg_o_fullmatint ,
                cg_s_fullmatint )

    # spectral thermo 
    G0 = [G for (G,(E,U)) in irrEU if E[1]==0][1]
    I0,S0 = G0[2:3]
    D0s = S0+1
    D0o = oindex2dimensions[I0]
    part0 = D0o*D0s

    push!(AA,
         redM2A(Mred,
                multiplets_a,
                cg_o_fullmatint,
                cg_s_fullmatint,
                irrEU,part0;
                verbose=true) 
         )
    return ( Mred , AA )

end

Mred, AA = update_redmat_AA(
            Mred ,
            irrEU ,
            combinations_uprima ,
            multiplets_a ,
            cg_o_fullmatint ,
            cg_s_fullmatint ,
            AA )


#   =============   #
#%% NRG PROCEDURE %%#
#   =============   #
hopchannels = collect(keys( hop ))
nrg = NRG_mat( 
           iterations , 
           cutoff_type , 
           cutoff_magnitude , 
           L , 
           hopchannels , 
           irrEU , 
           multiplets_shell ,
           cg_o_fullmatint , 
           pcg ,
           pcgmat ,
           qq_a ,
           combinations_uprima ,
           betabar ,
           oindex2dimensions ,
           mm_i ;
           spinarray=true ,
           cg_s_fullmatint=cg_s_fullmatint ,
           distributed=distributed ,
           z=z ,
           spectral=true ,
           M=Mred ,
           AA=AA ,
           etafac=1.0 );
println()


#   ===========   #
#%% PERFORMANCE %%#
#   ===========   #

println( "===========" )
println( "PERFORMANCE" )
println( "===========" )
print_performance_onestep(nrg)
println()


##   ==============   #
##%% SAVING TO FILE %%#
##   ==============   #
#
## impurity properties 
#if calculation=="IMP" 
#    write_impurity_info( nrg , omults , mult2index , orbital , z )
#end
#
## thermodata for this given value of z
#write_thermodata_onez( nrg , calculation , orbital , z )
#
## thermo diff
#if calculation=="IMP"
#    write_thermodiff( orbital , z )
#end


#   ~~~~~~~~~~~~~~~~~   #
#%% SPECTRAL FUNCTION %%#
#   ~~~~~~~~~~~~~~~~~   #
writedlm( "spectral/spectral.dat" , nrg.specfunc )
