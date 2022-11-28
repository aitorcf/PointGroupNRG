#!/usr/bin/env julia

#   ======================   #
#%% CALCULATION PARAMETERS %%#
#   ======================   #   

orbital = "Eg"
label = "EgPCGRED"
norbitals = 2

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

# clean system or with impurity added
calculation = "IMP"

# twisting parameter
z = 0.0

# numerical parameters
L = 10.0
betabar = 1.0

# cutoff
cutoff_type = "multiplet" 
cutoff_magnitude = 100
minmult = 0 
mine = 0.0

# command-line input
u_11_r = 0.7
u_h_r  = 0.5
gam_r  = 0.3
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
eps  = -0.1
gam = gam_r*abs(eps) # gam = √2Δ/π 
# coulomb coupling parameters
u_11 = u_11_r*abs(2*eps)
u_12 = u_11 
u_h  = u_h_r*( 2*u_12 )
u_eg = u_11
u_a1 = u_12 + u_h/2.0
u_a2 = u_12 - u_h/2.0

iterations = 30

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
using SparseArrays

include( "modules/symbols.jl" )
include( "modules/numericals.jl" )
include( "modules/compoundoperators.jl" )
include( "modules/shell.jl" )
include( "modules/thermo.jl" )
include( "modules/discretization.jl" )
include( "modules/spectral.jl" ) 
include( "modules/reddiag.jl" )
include( "modules/diagonalization.jl" )
include( "modules/automatization.jl" )

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
        include( "modules/spectral.jl" ) 
        include( "modules/reddiag.jl" )
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
# orbital irreps present in the atom
atom_orbital_irreps = ( "Eg" , )
# atomic configuration: 1 Eg irrep
atom_config = Dict( "Eg" => 1 )
# identity irrep 
identityrep = "A1g"
# directory with information about antisymmetric combinations
asym_dir = "/home/aitor/Bulegoa/AntiSymmetricPart/Oh/";
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
cg_s_fullmatint = get_cg_s_fullmatint( max_spin2 );

# symstates, basis and multiplets
symstates_0,basis_0,multiplets_0,multiplets_a = 
    get_symstates_basis_multiplets( 
            atom_config,
            oirreps2dimensions,
            identityrep,
            asym_dir,
            cg_o_dir ;
            verbose=true )

#   ---------   #
#%% operators %%#
#   ---------   #
epsilon_symparams = Dict( 
    "Eg" => eps
)
epsilon = epsilon_sym( symstates_0 , epsilon_symparams ; verbose=false )
println( "EPSILON OPERATOR" )
print_diagonals( epsilon )


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
H   = rescale( H ,   L , z , discretization )
gam = rescale( gam , L , z , discretization )
print_diagonals(H)


#   -----   #
#%% irreu %%#
#   -----   #
irrEU = get_irrEU_initial(symstates_0,H,calculation;verbose=true)
print_spectrum( irrEU )

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
            verbose=true )

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
multiplets_block = get_multiplets_block(
                        calculation ,
                        mints_0 ,
                        identityrep ,
                        oirreps2indices )
omults = ordered_multiplets(multiplets_block)
mult2index = Dict( m=>i for (i,m) in 
                   enumerate(omults))
multiplets_shell = mints_0
println( "BLOCK MULTIPLETS" )
for m in multiplets_block 
    println( m )
end
println()
println( "SHELL MULTIPLETS" )
for m in multiplets_shell 
    println( m )
end
println()

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
            verbose=true )

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
mm_i,m_imp = setup_impmultinfo( 
                multiplets_block ,
                irrEU ,
                betabar ,
                oindex2dimensions )

#multiplets_shell = filter( x->(x[1]!==0 && x[1]!==4) , multiplets_shell )
#@show multiplets_shell

# ~~~~~~~~~~~~~~~~~~ #
# PRECOMPUTE CG SUMS #
# ~~~~~~~~~~~~~~~~~~ #
Bsum_o_dict,Bsum_s_dict,Csum_o_dict,Csum_s_dict =
    precompute_CGsums(
            oirreps ,
            multiplets_a ,
            multiplets_shell ,
            max_spin2 ,
            oindex2dimensions ,
            cg_o_fullmatint ,
            cg_s_fullmatint )
Bsum_o_array,Bsum_s_array,Csum_o_array,Csum_s_array = 
    CGsums_dict2array( Bsum_o_dict,
                       Bsum_s_dict,
                       Csum_o_dict,
                       Csum_s_dict ) 


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
                #Csum_o_dict ,
                #Csum_s_dict ,
                #Bsum_o_dict ,
                #Bsum_s_dict ,
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
                   mine=mine ,
                   distributed=distributed ,
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
