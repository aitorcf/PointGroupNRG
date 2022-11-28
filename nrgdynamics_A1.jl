#!/usr/bin/env julia

#   ==========   #
#%% PARAMETERS %%#
#   ==========   #

orbital = "A1"
label = "A1"
norbitals = 1 

identityrep = "A1"

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
# calculation
calculation = "IMP"

z = 0.0
etafac = 1.0

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
cou =  0.2
gam =  0.1
#gam = 0.0

L = 2.0
betabar = 1.0

cutoff_type = "multiplet" # multiplet/energy
cutoff_magnitude = 50

iterations = 100

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
include( "modules/automatization.jl" )
include( "modules/reddiag.jl" )
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


# ~~~~~ #
# SETUP #
# ~~~~~ #
# directory where the orbital cg info is stored
cg_o_dir = "/home/aitor/Bulegoa/ClebschGordan/C4v/cg_symbolic/"
# orbital irreps present in the atom
atom_orbital_irreps = ( "A1" , )
# atomic configuration: 1 Eg irrep
atom_config = Dict( "A1" => 1 )
# identity irrep 
identityrep = "A1"
# directory with information about antisymmetric combinations
asym_dir = "/home/aitor/Bulegoa/AntiSymmetricPart/C4v/";
# hiztegia
hiztegia = Dict{String,Any}( o=>o for (o,_) in atom_config )
merge!( hiztegia , Dict( "u"=>0.5, "d"=>-0.5 ) )


#   ===========   #
#%% ATOMIC PART %%#
#   ===========   #

# symmetry-related information
print( "Obtaining Clebsch-Gordan coefficients... " )
@time begin 
(oirreps,
 oirreps2indices,
 oirreps2dimensions,
 oindex2dimensions,
 cg_o_fullmatint) = get_cg_o_info( cg_o_dir , atom_orbital_irreps  )
cg_s_fullmatint = get_cg_s_fullmatint( max_spin2 );
end; println()

# symstates, basis and multiplets
print( "Constructing atomic states... " )
@time symstates_0,basis_0,multiplets_0,multiplets_a = 
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
    "A1" => eps
)
epsilon = epsilon_sym( symstates_0 , epsilon_symparams ; verbose=false )
println( "EPSILON OPERATOR" )
print_diagonals( epsilon )


u_symparams = Dict( 
    ("A1",0) => [cou][:,:]
)
coulomb = u_sym( symstates_0 , u_symparams ; verbose=false )
println( "COULOMB OPERATOR" )
print_diagonals( coulomb )
println()

#   ---------------------------   #
#%% atomic hamiltonian (scaled) %%#
#   ---------------------------   #
H = epsilon + coulomb 
alpha = discretization=="co2005" ? compute_ebar0_z(z,L;discretization=discretization) : 0.5 * L^z * (1+L^-1)
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
    (oirreps2indices["A1"],1) => ComplexF64(gam) 
)
println( "ATOM-SHELL HOPPING PARAMETER" )
@show hop
hopchannels = collect(keys(hop))
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

#   ------------------------   #
#%% impurity quantum numbers %%#
#   ------------------------   #
mm_i,m_imp = setup_impmultinfo( 
                multiplets_block ,
                irrEU ,
                betabar ,
                oindex2dimensions )

#   -------------------   #
#%% clebsch-gordan sums %%#
#   -------------------   #
print( "Precomputing Clebsch-Gordan sums..." )
@time begin
Bsum_o_dict,Bsum_s_dict,Csum_o_dict,Csum_s_dict =
    precompute_CGsums(
            oirreps ,
            multiplets_a,
            multiplets_block ,
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
(Karray_orbital,Karray_spin) = compute_Ksum_arrays(
            oindex2dimensions,
            cg_o_fullmatint,
            cg_s_fullmatint,
            1 ,
            [1] ,
            max_spin2 ,
            1 )
end

#   ------------------ #
#%% reduced pcg matrix #
#   ------------------ #
pcg_atom = get_pseudoCG( symstates_0 , 
                         basis_0 , 
                         hiztegia , 
                         oirreps2indices )
Mred, AA = setup_redmat_AA(
            pcg_atom,
            multiplets_block ,
            multiplets_a ,
            cg_o_fullmatint ,
            cg_s_fullmatint ,
            irrEU ;
            verbose=false )

#   --------------------------------------- #
#%% diagonalization: atom + innermost shell #
#   --------------------------------------- #
print( "Diagonlizing atom + innermost shell... " )
@time begin
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

# impurity thermodynamics 
mm_i,m_imp = update_impmultinfo( 
                mm_i ,
                irrEU ,
                betabar ,
                oindex2dimensions ,
                combinations_uprima )

Mred, AA = update_redmat_AA_CGsummethod(
            Mred ,
            irrEU ,
            combinations_uprima ,
            collect(multiplets_a) ,
            cg_o_fullmatint ,
            cg_s_fullmatint ,
            Karray_orbital ,
            Karray_spin ,
            AA ;
            verbose=false )
end #timing
println()


#   =============== #
#%% NRG CALCULATION #
#   =============== #

@time begin
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
                   method=method ,
                   z=z ,
                   discretization=discretization ,
                   verbose=false ,
                   spectral=true ,
                   M=Mred ,
                   AA=AA , 
                   etafac=etafac ,
                   Karray_orbital=Karray_orbital ,
                   Karray_spin=Karray_spin ,
                   multiplets_atomhop=collect(multiplets_a) ,
                   alpha=Float64(alpha) )
end


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

writedlm( "spectral/spectral.dat" , nrg.specfunc )

##   ===========================   #
##&& SYMBOLIC CALCULATION OF N=0 &&#
##   ===========================   #
#
##   ------------ #
##%% state tuples # 
##   ------------ #
#a_1 = ( 0 , "a" , 1 )
#a_1u = ( a_1... , "u" )
#a_1d = ( a_1... , "d" )
#
#tuples_0 = [ a_1u , a_1d ]
#
##   ----- #
##%% basis #
##   ----- #
#hilbert_0 = HilbertSpace( tuples_0... )
#basis_0 = CanonicalBasis( hilbert_0 )
#println( "HILBERT SPACE" )
#println( hilbert_0 )
#println()
#println( "BASIS" )
#println( basis_0 )
#
#
##   --------- #
##%% symstates #
##   --------- #
#hiztegia = Dict( 
#    "a" =>  "A1",
#    "u" =>  0.5 ,
#    "d" => -0.5
#)
#
#symstates_0 = oneirrep_symstates( hilbert_0 , hiztegia , "A1" , "$(ASYM_PATH)A1/" )
#symstates_0_new = Dict()
#for (q,s) in symstates_0 
#    push!( symstates_0_new , (q[1:5]...,1)=>s ) 
#end
#symstates_0 = symstates_0_new
#println( "SYMSTATES" )
#print_dict( symstates_0 )
#
#multiplets_0 = get_multiplets( symstates_0 )
#
#irrepsubs_0 = get_multiplets( symstates_0 )
#println( "MULTIPLETS" )
#for m in irrepsubs_0 
#    println( m )
#end
#println()
#
##   --------- #
##%% operators #
##   --------- #
#epsilon_symparams = Dict( 
#    "A1" => eps
#)
#epsilon = epsilon_sym( symstates_0 , epsilon_symparams ; verbose=false )
#println( "EPSILON OPERATOR" )
#print_diagonals( epsilon )
#
#u_symparams = Dict( 
#    ("A1",0) => [cou][:,:]
#)
#coulomb = u_sym( symstates_0 , u_symparams ; verbose=false );
#println( "COULOMB OPERATOR" )
#print_diagonals( coulomb )
#
#
##   ---------------- #
##%% full hamiltonian #
##   ---------------- #
#H = epsilon + coulomb 
## α = ϵ_0^z (scaling factor)
#if discretization=="standard"
#    global α = 0.5 * L^z * (1+L^-1)
#elseif discretization=="co2005"
#    global α = compute_ebar0_z(z,L;discretization=discretization)
#end
#H.matrix ./= α
#gam = sqrt(2.0*gam/pi) / α # for later
##gam /= α
#
#
##   ================================================   #
##&& NUMERICAL ADDITION OF INNERMOST CONDUCTION SHELL &&#
##   ================================================   #
#
#function normalize_irrEU( irrEU )
#    minE = minimum([e for (G,(E,U)) in irrEU for e in E])
#    return Dict( G=>(E.-minE,U) for (G,(E,U)) in irrEU )
#end
#
##   -----   #
##%% irreu %%#
##   -----   #
## imp
#if calculation=="IMP"
#    irreps_0 = get_irreps( symstates_0 )
#    irrEU_imp = symdiag( irreps_0 , symstates_0 , H )
#    irrEU_imp = normalize_irrEU( irrEU_imp )
#    println( "irrEU for IMP" )
#    print_dict( irrEU_imp )
#    println()
## clean
#elseif calculation=="CLEAN"
#    irrEU_clean = get_irrEU_clean( "A1g" )
#    println( "irrEU for CLEAN" )
#    print_dict( irrEU_clean )
#    println()
#end
#
#
##   ---------------   #
##%% combinations u' %%#
##   ---------------   #
#combinations_uprima = 
#    Dict{ Tuple{Int64,Int64,Int64,Int64} , NTuple{2,Tuple{Int64,Int64,Int64,Int64}} }()
#m_vac = (0,oirreps2indices["A1"],0,1)
##print( "COMBINATIONS U' " )
#if calculation=="IMP"
#    #println( "for IMP" )
#    for m_mu in multiplets_0
#        mint_mu = convert_to_int( m_mu , oirreps2indices )
#        push!( combinations_uprima , mint_mu=>(mint_mu,m_vac) )
#    end 
#elseif calculation=="CLEAN" 
#    #println( "for CLEAN" )
#    push!( combinations_uprima , m_vac=>(m_vac,m_vac) )
#end
#
#irreps_uprima = Set( k[1:3] for k in keys(combinations_uprima) )
#combinations_uprima = 
#        Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} }(
#            G => NTuple{3,NTuple{4,Int64}}[
#                (m_u,m_mu,m_i)
#                for (m_u,(m_mu,m_i)) in combinations_uprima 
#                if m_u[1:3]==G
#                ]
#            for G in irreps_uprima
#        ) 
#println( "COMBINATIONS U' FOR N=0" )
#for (G,combs) in combinations_uprima
#    println( "$G => $combs" )
#end
#println()
#
##   ----------------------   #
##%% pseudo-CG coefficients %%#
##   ----------------------   #
#pcg = get_pseudoCG( symstates_0 , 
#                    basis_0 , 
#                    hiztegia , 
#                    oirreps2indices )
#println( "PCG" )
#print_dict( pcg ) 
#println()
#
## PCG matrix in int format
#irrmult_0 = get_irreps( multiplets_0 ; multiplicity=true )
#pcgmat = get_pseudoCG_mat( pcg , 
#                           irrmult_0 , 
#                           oindex2dimensions , 
#                           oirreps2indices );
#
##   -------------------   #
##%% shell cops and qq_a %%#
##   -------------------   #
#shell_cops = shell_coperators( basis_0 , hiztegia )
#qq_a = collect( convert_to_int(q_a,oirreps2indices) 
#                for q_a in keys(shell_cops) )
#println( "ONE-PARTICLE HOPPERS" )
#@show qq_a;
#
#
#
##   -----------------   #
##%% hopping parameter %%#
##   -----------------   #
#hop = Dict( 
#    (oirreps2indices["A1"],1) => ComplexF64(gam) 
#)
#println( "ATOM-SHELL HOPPING PARAMETER" )
#@show hop
#println()
#println()
#
##   --------------------------   #
##%% block and shell multiplets %%#
##   --------------------------   #
#if calculation=="IMP"
#    symstates_block = symstates_0 
#    multiplets_block = get_multiplets( symstates_block )
#    global irrEU = Dict( (convert_to_int(G,oirreps2indices),(E,U))
#                   for (G,(E,U)) in irrEU_imp )
#    println( "irrEU before adding first shell (to IMP)" )
#    for (G,(E,U)) in irrEU 
#        @show G, E
#    end
#    println()
#elseif calculation=="CLEAN"
#    multiplets_block = Set([(0,"A1",0.0,1)])
#    global irrEU = Dict( (convert_to_int(G,oirreps2indices),(E,U))
#                   for (G,(E,U)) in irrEU_clean )
#    println( "irrEU before adding first shell (to CLEAN)" )
#end 
#symstates_shell  = symstates_0 
#multiplets_shell = multiplets_0
#
#multiplets_block = Set( convert_to_int(m,oirreps2indices) 
#                        for m in multiplets_block )
#multiplets_shell = Set( convert_to_int(m,oirreps2indices)
#                        for m in multiplets_shell )
#
##   -----------------   #
##%% excitation matrix %%#
##   -----------------   #
#M = get_M0( pcg , qq_a )
#println( "===================" )
#println( "M excitation matrix" )
#println( "===================" )
#print_dict( M )
#println()
#println( "M size: $(length(collect(keys(M[qq_a[1]]))))" )
#AA = [get_A(M,irrEU)]
#
##   ------------------------   #
##%% impurity quantum numbers %%#
##   ------------------------   #
#println("IMPURITY MULTIPLET SPACE")
#omults = ordered_multiplets(multiplets_block)
#mult2index = Dict( m=>i for (i,m) in 
#                   enumerate(omults))
#mm_i = Dict( 
#            m=>[(i==mult2index[m] ? 1.0 : 0.0)
#                for i in 1:length(multiplets_block)] 
#                for m in omults
#           )
#@show mm_i 
#m_imp = mult_thermo( irrEU ,
#                     betabar ,
#                     oindex2dimensions ,
#                     mm_i )
#println()
#
##   ---------------------------------------   #
##%% matrix construction and diagonalization %%#
##   ---------------------------------------   #
#if distributed && onestep
#    @time begin
#    global (irrEU,combinations_uprima) = matdiag_distributed_mat( 
#                multiplets_block , 
#                multiplets_shell ,
#                irrEU , 
#                hop , 
#                cg_o_fullmatint , 
#                cg_s_fullmatint ,
#                pcg , 
#                pcgmat, 
#                qq_a , 
#                combinations_uprima , 
#                oindex2dimensions ;
#                verbose=false );
#    end
#else
#    @time begin
#    global (irrEU,combinations_uprima) = matdiag_serial_mat( 
#                multiplets_block , 
#                multiplets_shell ,
#                irrEU , 
#                hop , 
#                cg_o_fullmatint , 
#                cg_s_fullmatint ,
#                pcg , 
#                pcgmat, 
#                qq_a , 
#                combinations_uprima , 
#                oindex2dimensions ;
#                verbose=false );
#    end
#end
#minE = minimum([e for (E,U) in values(irrEU) for e in E])
#irrEU = Dict( G=>(E.-minE,U) for (G,(E,U)) in irrEU )
#energies = sort([e for (E,U) in values(irrEU) for e in E])
#@show energies
#@show length(energies)
#
#mm_i = imp_mults( irrEU ,
#                  oindex2dimensions ,
#                  combinations_uprima ,
#                  mm_i )
#@show mm_i
#m_imp = mult_thermo( irrEU ,
#                     betabar ,
#                     oindex2dimensions ,
#                     mm_i )
#
#
##   ~~~~~~~~~~~~~~~~~~~   #
##%% transformation of M %%#
##   ~~~~~~~~~~~~~~~~~~~   #
#@show combinations_uprima
#
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
#
#
##   =============   #
##%% NRG PROCEDURE %%#
##   =============   #
#hopchannels = collect(keys( hop ))
#nrg = GRG_mat_spectral( 
#           iterations , 
#           cutoff_type , 
#           cutoff_magnitude , 
#           L , 
#           hopchannels , 
#           irrEU , 
#           multiplets_shell ,
#           cg_o_fullmatint , 
#           pcg ,
#           pcgmat ,
#           qq_a ,
#           combinations_uprima ,
#           betabar ,
#           oindex2dimensions ,
#           mm_i ,
#           M ,
#           AA ;
#           spinarray=true ,
#           cg_s_fullmatint=cg_s_fullmatint ,
#           distributed=distributed ,
#           z=z );
#println()
#
#
##   ===========   #
##%% PERFORMANCE %%#
##   ===========   #
#
#println( "===========" )
#println( "PERFORMANCE" )
#println( "===========" )
#print_performance_onestep(nrg)
#println()
#
#
###   ==============   #
###%% SAVING TO FILE %%#
###   ==============   #
##
### impurity properties 
##if calculation=="IMP" 
##    write_impurity_info( nrg , omults , mult2index , orbital , z )
##end
##
### thermodata for this given value of z
##write_thermodata_onez( nrg , calculation , orbital , z )
##
### thermo diff
##if calculation=="IMP"
##    write_thermodiff( orbital , z )
##end
#
#
##   ~~~~~~~~~~~~~~~~~   #
##%% SPECTRAL FUNCTION %%#
##   ~~~~~~~~~~~~~~~~~   #
#
#println( "# ########################## #" )
#println( "# COMPUTING SPECTRAL FUNCION #" )
#println( "# ########################## #" )
#println()
#
#function is_in_interval( omega , emin , emax ) 
#    return (omega>emin && omega<emax)
#end
#
#@show α
#
#widthfac = L
#etafac = 1.0
#q_a = qq_a[1]
#omegas = [sign*L^(-(x-2)/2.0) for sign in [-1.0,1.0] for x in 2:0.1:(iterations-2)]
##omegas = collect(-1.0:0.1:1.0)
#spectral = [0.0 for o in omegas]
#sort!( omegas )
#@show q_a 
#@show omegas
#println( "COMPUTING SPECTRAL FUNCTION" )
#
#alternation = "even"
#alternation=="even" && (Ninit=3)
#alternation=="odd" && (Ninit=2)
#
#for (i,o) in enumerate(omegas)
#
#    omega = -o
#
#    println( "ω = $omega" )
#
#    for N in Ninit:2:iterations 
#        omegaN = Float64(α * L^(-(N-2)/2.0) )
#        emin   = omegaN#/widthfac
#        #emax  = omegaN*nrg.maxes[N-1]
#        emax   = omegaN*widthfac
#        p0     = nrg.partitions0[N-1]
#        A      = AA[N][q_a] 
#        eta    = etafac*omegaN
#
#        #positive interval 
#        if is_in_interval( omega , emin , emax )
#            println( "omega is in interval $(N): [$emin,$emax]" )
#            @show omegaN
#            #@show nrg.maxes[N-1]
#            @show p0
#            @show eta
#            println( "A run" )
#            for (m,(e,coeffs)) in A
#                @show m,e
#                Δ = (omega-e*omegaN)
#                @show Δ
#                c = p0^-1*
#                    coeffs[2]*
#                    P(Δ,etafac*omegaN)
#                @show c
#                spectral[i] += p0^-1*
#                               coeffs[2]*
#                               P(Δ,eta)
#            end
#            println()
#        end
#
#        # negative interval 
#        if is_in_interval( omega , -emax , -emin )
#            println( "omega is in interval $(N): [$(-emax),$(-emin)]" )
#            @show omegaN
#            @show nrg.maxes[N-1]
#            @show p0
#            @show eta
#            println( "A run" )
#            for (m,(e,coeffs)) in A
#                @show m,e
#                Δ = omega+e*omegaN
#                @show Δ
#                c = p0^-1*
#                    coeffs[2]*
#                    P(Δ,eta)
#                @show c
#                spectral[i] += p0^-1 *
#                               coeffs[1] *
#                               P(Δ,eta)
#            end
#            println()
#        end
#
#    end
#end
##for (i,o) in enumerate(omegas)
##    o>0 && (omegas[i]=log(o))
##    o<0 && (omegas[i]=-log(-o))
##end
#writedlm( "spectral/spectral.dat" , [omegas spectral] )
#
#
#
