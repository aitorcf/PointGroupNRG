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
L = 4.0
betabar = 1.0

# cutoff
cutoff_type = "multiplet" 
cutoff_magnitude = 100
minmult = 0

# command-line input
u_11_r = 1.0
u_h_r  = 1.0
gam_r  = 0.8
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
    # relative parameters
    if length(ARGS)>2
        u_11_r = parse(Float64,ARGS[3])
        u_h_r  = parse(Float64,ARGS[4]) 
        gam_r  = parse(Float64,ARGS[5])
    end
end

## standard parameters 
#U   = 0.2
#U12 = 0.8*U
#J   = 0.5*( U - U12 )
#u_eg =  U - J
#u_a1 =  U + J
#u_a2 =  U12 - J 

# one-body coupling parameters
eps  = -0.2
gam = gam_r*abs(eps) # gam = √2Δ/π 
# coulomb coupling parameters
u_11 = u_11_r*abs(2*eps)
u_12 = u_11 
u_h  = u_h_r*( 2*u_12 )
u_eg = u_11
u_a1 = u_12 + u_h/2.0
u_a2 = u_12 - u_h/2.0

iterations = 50

max_spin2 = 8

distworkers = 6

println( "====================" )
println( "SETUP AND PARAMETERS" )
println( "====================" )
@show calculation
@show distributed 
distributed && @show distworkers
@show method
@show discretization
@show z
@show eps   
@show u_eg 
@show u_a1 
@show u_a2 
@show gam  
@show L
@show betabar
@show cutoff_type
@show cutoff_magnitude
@show iterations
@show max_spin2
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
include( "modules/thermo.jl" )
include( "modules/discretization.jl" )
include( "modules/spectral.jl" )

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
        include( "modules/diagonalization.jl" )
    end

else 

    println( "SERIAL CALCULATION" )

    using ProgressMeter
    using PartialWaveFunctions
    using StaticArrays

    include( "modules/symmetry.jl" )
    include( "modules/diagonalization.jl" )

end
println()


#   ======================   #
#%% ORBITAL CLEBSCH-GORDAN %%#
#   ======================   #

CG_PATH = "/home/aitor/Bulegoa/ClebschGordan/Oh/cg_symbolic/"
ASYM_PATH = "/home/aitor/Bulegoa/AntiSymmetricPart/Oh/";

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
gam /= α
#gam = sqrt(2.0*gam/pi) / α # for later
#println( "H OPERATOR" )
#print_diagonals( H )
#println()

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
    println( "irrEU for IMP" )
    print_dict( irrEU_imp )
    println()
# clean
elseif calculation=="CLEAN"
    irrEU_clean = get_irrEU_clean( "A1g" )
    println( "irrEU for CLEAN" )
    print_dict( irrEU_clean )
    println()
end


#   ---------------   #
#%% combinations u' %%#
#   ---------------   #
combinations_uprima = 
    Dict{ Tuple{Int64,Int64,Int64,Int64} , NTuple{2,Tuple{Int64,Int64,Int64,Int64}} }()
m_vac = (0,oirreps2indices["A1g"],0,1)
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
println()

#   -----------------   #
#%% hopping parameter %%#
#   -----------------   #
hop = Dict( 
    (oirreps2indices["Eg"],1) => ComplexF64(gam) 
)
println( "ATOM-SHELL HOPPING PARAMETER" )
@show hop
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
    multiplets_block = Set([(0,"A1g",0.0,1)])
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
multiplets_a = collect(filter( x->x[1]==1 , multiplets_block ))
Mred, AA = setup_redmat_AA(
            pcg ,
            multiplets_block ,
            multiplets_a ,
            cg_o_fullmatint ,
            cg_s_fullmatint ,
            irrEU ;
            verbose=false )

#   -------------------------   #
#%% effective multiplet model %%#
#   -------------------------   #

#(multiplets_block, discarded) = cut_off!( 
#                                    irrEU ; 
#                                    type="multiplet" , 
#                                    cutoff=2 , 
#                                    safeguard=false ,
#                                    verbose=true )

#   ------------------------   #
#%% impurity quantum numbers %%#
#   ------------------------   #
#println("IMPURITY QUANTUM NUMBERS")
#ss_i = Dict( m_b=>(m_b[3]/2.0*(m_b[3]/2.0+1)) 
#             for m_b in multiplets_block )
#nn_i = Dict( m_b=>m_b[1] for m_b in multiplets_block )
#s_imp = impspin( irrEU ,
#                 betabar , 
#                 oindex2dimensions ,
#                 ss_i )
#n_imp = impnum( irrEU ,
#                betabar , 
#                oindex2dimensions ,
#                nn_i )
#@show ss_i
#@show nn_i
#@show s_imp 
#@show n_imp
#println()

function ordered_multiplets( mults ) 
    max_N = maximum(k[1] for k in mults) 
    omults = []
    for N in 0:max_N 
        for mult in mults
            mult[1]==N || continue
            push!(omults,mult)
        end
    end
    return omults
end

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
m_imp = mult_thermo( irrEU ,
                     betabar ,
                     oindex2dimensions ,
                     mm_i )

#   ~~~~~~~~~~~~~~~~~~~   #
#%% transformation of M %%#
#   ~~~~~~~~~~~~~~~~~~~   #
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
nrg = NRG_mat( iterations , 
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
           spinarray=spinarray ,
           cg_s_fullmatint=cg_s_fullmatint ,
           distributed=parallel ,
           minmult=minmult ,
           discretization=discretization ,
           z=z ,
           spectral=true ,
           M=Mred ,
           AA=AA ,
           etafac=0.8 );
println()


#   ===========   #
#%% PERFORMANCE %%#
#   ===========   #

println( "===========" )
println( "PERFORMANCE" )
println( "===========" )
print_performance_onestep(nrg)
println()

#   ==============   #
#%% SAVING TO FILE %%#
#   ==============   #

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

writedlm( "spectral/spectral.dat" , nrg.specfunc ) 
