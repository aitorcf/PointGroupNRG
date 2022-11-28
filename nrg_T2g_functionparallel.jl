#!/usr/bin/env julia

using Distributed
using DelimitedFiles
using ProfileVega
using Profile

include( "modules/symbols.jl" )
include( "modules/numericals.jl" )
include( "modules/compoundoperators.jl" )
include( "modules/shell.jl" )
include( "modules/thermo.jl" )
include( "modules/discretization.jl" )

orbital = "T2g"

# ~~~~~~~~~~~~~~~~~~~~ #
# Z-AVERAGE PARAMETERS #
# ~~~~~~~~~~~~~~~~~~~~ #

Nz = 2
distworkers = 2*Nz
Z = [ (0.0 + Float64(m)/Nz) for m=0:(Nz-1) ]


# ~~~~~~~~~~~~~~~~~ #
# DISTRIBUTED SETUP #
# ~~~~~~~~~~~~~~~~~ #

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

# load modules 
@everywhere begin 
    using ProgressMeter
    using PartialWaveFunctions
    using StaticArrays
    include( "modules/symbols.jl" )
    include( "modules/numericals.jl" )
    include( "modules/compoundoperators.jl" )
    include( "modules/shell.jl" )
    include( "modules/thermo.jl" )
    include( "modules/discretization.jl" )
    include( "modules/symmetry.jl" )
    include( "modules/diagonalization.jl" )
end


# ~~~~~~~~~~~~~~ #
# ONE-Z FUNCTION #
# ~~~~~~~~~~~~~~ #

#@everywhere function nrg1z( z , calculation ) 
function nrg1z( z , calculation ) 
    
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
    # numerical parameters
    L = 10.0
    betabar = 1.0
    iterations = 20
    # cutoff
    cutoff_type = "multiplet" 
    cutoff_magnitude = 10
    minmult = 0
    # one-particle parameters
    eps  = -0.1 
    gam  =  0.01
    # relative coulomb parameters 
    u_11 = 1.0
    u_12 = 1.0
    j    = 0.0
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
    println( "Computing PCG..." )
    @time pcg = get_pseudoCG( symstates_0 , 
                              basis_0 , 
                              hiztegia , 
                              oirreps2indices )
    println()
    #println( "PCG" )
    #print_dict( pcg ) 
    #println()
    
    # PCG matrix in int format
    irrmult_0 = get_irreps( multiplets_0 ; multiplicity=true )
    println( "Computing PCG matrix..." )
    @time pcgmat = get_pseudoCG_mat( pcg , 
                                     irrmult_0 , 
                                     oindex2dimensions , 
                                     oirreps2dimensions );
    println()
    
    
    #   -------------------   #
    #%% shell cops and qq_a %%#
    #   -------------------   #
    println( "Computing shell creation operators..." )
    @time shell_cops = shell_coperators( basis_0 , hiztegia )
    println()
    qq_a = collect( convert_to_int(q_a,oirreps2indices) 
                    for q_a in keys(shell_cops) )
    println( "ONE-PARTICLE HOPPERS" )
    @show qq_a;
    println()
    
    
    #   -----------------   #
    #%% hopping parameter %%#
    #   -----------------   #
    hop = Dict( 
        (oirreps2indices["T2g"],1) => ComplexF64(gam)
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
    
    #   ------------------------   #
    #%% impurity quantum numbers %%#
    #   ------------------------   #
    
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
    println( "Solving atom + innermost shell..." )
    if distributed
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
    
    mm_i = imp_mults( irrEU ,
                      oindex2dimensions ,
                      combinations_uprima ,
                      mm_i )
    m_imp = mult_thermo( irrEU ,
                         betabar ,
                         oindex2dimensions ,
                         mm_i )
    println()
    
    #println( "irrEU after adding the first shell to $calculation" )
    #println( "**********************************" )
    #for (G,(E,U)) in irrEU 
    #    irrEU[G] = (E./sqrt(L),U)
    #    @show G, E
    #end
    #println()
    
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
               z=z ,
               discretization=discretization ,
               distributed=parallel )

    # impurity properties 
    if calculation=="IMP" 
        write_impurity_info( nrg , omults , mult2index , orbital , z )
    end

    return Dict( (z,calculation)=>nrg )

end


#- %%%%%%%%%%%%%%% -#
#-                 -#
#- MAIN Z-AVG LOOP -#
#-                 -#
#- %%%%%%%%%%%%%%% -#
calcs = [(z,c) for z in Z, c in ["CLEAN","IMP"]]
#nrg_zavg = @sync @distributed (merge) for calc in calcs
#    nrg1z( calc... )
#end
nrg_zavg = Dict( nrg1z(calc...) for calc in calcs )


#   ==============   #
#%% SAVING TO FILE %%#
#   ==============   #

for calculation in ["CLEAN","IMP"], z in Z

    nrg = nrg_zavg[(z,calculation)]

    # thermodata for this given value of z
    write_thermodata_onez( nrg , calculation , orbital , z )

    # thermo diff
    if calculation=="IMP"
        write_thermodiff( orbital , z )
    end
end

#   ===============================   #
#%% POST-CALCULATION THERMODYNAMICS %%#
#   ===============================   #

# differences
for z in Z
    write_thermodiff( orbital , z )
end

# averages
println( "Averaging over values of z..." )
th_tot  = Dict()
for z in Z
    @show z
    th_z   = readdlm( "thermodata/th_diff_$(orbital)_z$z.dat" )
    t = th_z[:,1] 
    th_tot[z] = Dict( round(t[i],sigdigits=3)=>th_z[i,:] for i in 1:length(t) )
end
th_zavg = Dict()
T = sort(collect(keys(th_tot[Z[end]])))
T = T[Nz:(length(T)-(Nz-1))]
@show sort(collect(T))
for t in T
    th_zavg[t]  = sum( th_tot[z][t] for z in Z )/length(Z)
end
println()
th_zavg_vec = [th_zavg[t] for t in T]
writedlm( "thermodata/th_zavg_$(orbital).dat" , th_zavg_vec )
