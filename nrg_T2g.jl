#!/usr/bin/env julia


#   ======================   #
#%% CALCULATION PARAMETERS %%#
#   ======================   #   

orbital = "T2g"
label = "T2g"

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

# command-line input
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

# numerical parameters
L = 10.0
betabar = 1.0

iterations = 30

# cutoff
cutoff_type = "multiplet" 
cutoff_magnitude = 200
minmult = 100
mine = 1.5

# one-particle parameters
eps  = -0.5
gam  =  0.1
# relative coulomb parameters 
u_11 = 10.0
u_12 = 0.9
j    = 0.1
# absolute coulomb parameters
U_11 = u_11*2*abs(eps)
U_12 = u_12*U_11
J    = j*U_12

# multiplet coulomb parameters
u_a1g = U_11 + 2*J # on-site, orbital-symmetric
u_eg  = U_11 - J   # on-site, orbital-eg 
u_t1g = U_12 - J   # inter-site, spin-aligned
u_t2g = U_12 + J   # inter-site, spin-antialigned

max_spin2 = 10;

distworkers  = 6

println( "====================" )
println( "SETUP AND PARAMETERS" )
println( "====================" )
@show calculation
@show z
@show distributed 
distributed && @show distworkers 
@show method
@show eps   
@show u_11 
@show u_12 
@show j    
@show U_11 
@show U_12 
@show J    
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
include( "modules/thermo.jl" )
include( "modules/symmetry.jl" );
include( "modules/reddiag.jl" )
include( "modules/automatization.jl" )

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


# ~~~~~~~~~~~~~~~~~~~~~~~~ #
# EXTERNAL INPUT VARIABLES #
# ~~~~~~~~~~~~~~~~~~~~~~~~ #
# directory where the orbital cg info is stored
cg_o_dir = "/home/aitor/Bulegoa/ClebschGordan/Oh/cg_symbolic/"
# directory with information about antisymmetric combinations
asym_dir = "/home/aitor/Bulegoa/AntiSymmetricPart/Oh/";
# orbital irreps present in the atom
atom_orbital_irreps = ( "T2g" , )
# atomic configuration: 1 Eg multiplet
atom_config = Dict( "T2g" => 1 )
# identity irrep 
identityrep = "A1g"
# hiztegia
hiztegia = Dict{String,Any}( o=>o for (o,_) in atom_config )
merge!( hiztegia , Dict( "u"=>0.5, "d"=>-0.5 ) )


#   ===========   #
#%% ATOMIC PART %%#
#   ===========   #

# symmetry information
(oirreps,
 oirreps2indices,
 oirreps2dimensions,
 oindex2dimensions,
 cg_o_fullmatint) = get_cg_o_info( cg_o_dir , atom_orbital_irreps )
cg_s_fullmatint = get_cg_s_fullmatint( max_spin2 )

# symstates, basis and multiplets
symstates_0,basis_0,multiplets_0,multiplets_a = 
    get_symstates_basis_multiplets( 
            atom_config,
            oirreps2dimensions,
            identityrep,
            asym_dir,
            cg_o_dir ;
            verbose=false )


#   ---------   #
#%% operators %%#
#   ---------   #
epsilon_symparams = Dict( 
    orbital => eps
)
epsilon = epsilon_sym( symstates_0 , epsilon_symparams ; verbose=false )
println( "EPSILON OPERATOR" )
print_diagonals( epsilon )
println()


u_symparams = Dict( 
    ("A1g", 0) => [u_a1g][:,:],
    ("Eg",  0) => [u_eg ][:,:],
    ("T1g", 1) => [u_t1g][:,:],
    ("T2g", 0) => [u_t2g][:,:]
)
coulomb = u_sym( symstates_0 , u_symparams ; verbose=false )
println( "COULOMBS OPERATOR" )
print_diagonals( coulomb )
println()


#   ---------------------------   #
#%% atomic hamiltonian (scaled) %%#
#   ---------------------------   #
H = epsilon + coulomb 
H   = rescale( H ,   L , z , discretization )
gam = rescale( gam , L , z , discretization )

#   -----   #
#%% irreu %%#
#   -----   #
irrEU = get_irrEU_initial(symstates_0,H,calculation;verbose=false)
print_spectrum( irrEU )
x=y

#   -----------------   #
#%% conversion to int %%#
#   -----------------   #
irrEU = irrEU2int( irrEU , oirreps2indices ) 
mints_0 = multiplets2int( multiplets_0 , oirreps2indices )
multiplets_a = multiplets2int( multiplets_a , oirreps2indices )


#   ================================================   #
#&& NUMERICAL ADDITION OF INNERMOST CONDUCTION SHELL &&#
#   ================================================   #

#   ---------------   #
#%% combinations u' %%#
#   ---------------   #
combinations_uprima = get_combinations_uprima_initial(
            identityrep ,
            calculation ,
            mints_0 ,
            oirreps2indices ;
            verbose=false )

#   -----------------   #
#%% hopping parameter %%#
#   -----------------   #
hop = Dict( 
    (oirreps2indices[orbital],1) => ComplexF64(gam) 
)

#   --------------------------   #
#%% block and shell multiplets %%#
#   --------------------------   #
multiplets_block = get_multiplets_block(
                        calculation ,
                        mints_0 ,
                        identityrep ,
                        oirreps2indices )
multiplets_shell = mints_0
println( "BLOCK/SHELL MULTIPLETS" )
print_multiplets_Nordered( multiplets_block )

#   ----------------------   #
#%% pseudo-CG coefficients %%#
#   ----------------------   #
pcgred_atom,pcgred_shell = prepare_pcgred( 
            symstates_0 ,
            basis_0 ,
            hiztegia ,
            oirreps2indices ,
            cg_o_fullmatint ,
            cg_s_fullmatint ,
            multiplets_block ,
            multiplets_shell ;
            verbose=false )

#   -------------------------   #
#%% effective multiplet model %%#
#   -------------------------   #

print( "Cutting off multiplets...\nNEW " )
(multiplets_block, discarded) = cut_off!( 
                                    irrEU ; 
                                    type="multiplet" , 
                                    cutoff=2 , 
                                    safeguard=false ,
                                    verbose=true )
print_spectrum( irrEU )
omults = ordered_multiplets(multiplets_block)
mult2index = Dict( m=>i for (i,m) in 
                   enumerate(omults))

#   ------------------------   #
#%% impurity quantum numbers %%#
#   ------------------------   #
mm_i,m_imp = setup_impmultinfo( 
                multiplets_block ,
                irrEU ,
                betabar ,
                oindex2dimensions )


# ~~~~~~~~~~~ #
#%% cg sums %%#
# ~~~~~~~~~~~ #
Bsum_o_dict,Bsum_s_dict,Csum_o_dict,Csum_s_dict =
    precompute_CGsums(
            oirreps ,
            multiplets_a ,
            multiplets_shell ,
            max_spin2 ,
            oindex2dimensions ,
            cg_o_fullmatint ,
            cg_s_fullmatint ;
            verbose=false )
Bsum_o_array,Bsum_s_array,Csum_o_array,Csum_s_array = 
    CGsums_dict2array( Bsum_o_dict,
                       Bsum_s_dict,
                       Csum_o_dict,
                       Csum_s_dict ) 

println( "@ irrEU @" )
println( "=========" )
print_dict( irrEU )

#   ---------------------------------------   #
#%% matrix construction and diagonalization %%#
#   ---------------------------------------   #
(irrEU,combinations_uprima) = matdiag_redmat( 
                multiplets_block , 
                multiplets_shell ,
                irrEU , 
                hop , 
                cg_o_fullmatint , 
                cg_s_fullmatint ,
                Csum_o_array ,
                Csum_s_array ,
                Bsum_o_array ,
                Bsum_s_array ,
                pcgred_atom ,
                pcgred_shell ,
                collect(multiplets_a) , 
                combinations_uprima , 
                oindex2dimensions ;
                verbose=false ,
                distributed=distributed );
mm_i,m_imp = update_impmultinfo( 
                mm_i ,
                irrEU ,
                betabar ,
                oindex2dimensions ,
                combinations_uprima )

print( "AFTER FIRST DIAGONALIZATION, " )
print_spectrum( irrEU )


#   =============   #
#%% NRG PROCEDURE %%# 
#   =============   #
hopchannels = collect(keys( hop ))
nrg = NRG_pcgred( iterations,
                   cutoff_type,
                   cutoff_magnitude,
                   L,
                   hopchannels,
                   irrEU,
                   multiplets_shell,
                   cg_o_fullmatint,
                   cg_s_fullmatint,
                   Csum_o_array ,
                   Csum_s_array ,
                   Bsum_o_array ,
                   Bsum_s_array ,
                   pcgred_shell,
                   collect(multiplets_a), 
                   combinations_uprima,
                   betabar,
                   oindex2dimensions,
                   mm_i ;
                   distributed=distributed ,
                   mine=mine ,
                   method=method ,
                   z=z ,
                   discretization=discretization ,
                   verbose=false )
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

# impurity properties 
if calculation=="IMP" 
    write_impurity_info( nrg , omults , mult2index , label , z )
end

# thermodata for this given value of z
write_thermodata_onez( nrg , calculation , label , z )

# thermo diff
if calculation=="IMP"
    write_thermodiff( label , z )
end
println()
