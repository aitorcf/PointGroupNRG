#!/usr/bin/env julia

#using Pkg
#Pkg.activate("/scratch/aitorcalvo/juliaenv")
#Pkg.instantiate()

#   ======================   #
#%% CALCULATION PARAMETERS %%g
#   ======================   #   

orbital = "B2DIMER"
norbitals = 2

#run = "spectral" 
run = "thermo" 

# construction and diagonalization in a single step
onestep = true
# use precomputed spin-CG coefficients
spinarray = true
# parallel computation
distributed = parallel = false
# parallel method: distfor or async
method = "distfor"
# discretization ("standard" or "co2005")
discretization = "standard"

# clean system or with impurity added
calculation = "IMP"

# twisting parameter
z = 0.0

# numerical parameters
L = 10.0
betabar = 1.0
etafac = 1.0

# cutoff
cutoff_type = "multiplet" 
cutoff_magnitude = 100
minmult = 0

# free molecule spectrum
# neutral
E_0  =  0.0
E_1  =  32e-3 # discard?
# positive charge
E_0p = 467e-3
E_1p = 507e-3
E_2p = 534e-3
E_3p = 634e-3 # discarded
# negative charge
E_0n = -47e-3
E_1n = -41e-3
E_2n = +306e-3 # discarded 
E_3n = +380e-3 # discarded

# adsorbed molecule spectrum
mu = -200e-3
E_0p += mu
E_1p += mu
E_2p += mu
E_3p += mu
E_0n -= mu 
E_1n -= mu 
E_2n -= mu
E_3n -= mu

# hybridization
gam = 0.0
gam = sqrt(2*gam/pi)

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

iterations = 80

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
@show cutoff_type
@show cutoff_magnitude
@show iterations
@show max_spin2
@show betabar
@show L
@show E_0
@show E_1 
@show E_0p
@show E_1p
@show E_0n
@show E_1n
@show gam  
println()


#   =======================   #
#%% MODULES AND DISTRIBUTED %%#
#   =======================   #

using DelimitedFiles
using Profile

moduledir = "/home/aitor/Bulegoa/PointGroupNRG/modules" 

include( "$(moduledir)/symbols.jl" )
include( "$(moduledir)/numericals.jl" )
include( "$(moduledir)/compoundoperators.jl" )
include( "$(moduledir)/shell.jl" )
include( "$(moduledir)/thermo.jl" )
include( "$(moduledir)/diagonalization.jl" )
include( "$(moduledir)/spectral.jl" )
include( "$(moduledir)/reddiag.jl" ) 
include( "$(moduledir)/automatization.jl" )

if parallel 

    using Distributed 

    #1 kill current processes
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
        include( "$(moduledir)/symmetry.jl" )
        include( "$(moduledir)/diagonalization.jl" )
    end

else 

    println( "SERIAL CALCULATION" )

    using ProgressMeter
    using PartialWaveFunctions
    using StaticArrays

    include( "$(moduledir)/symmetry.jl" )
    include( "$(moduledir)/diagonalization.jl" )

end
println()


#   ==============   #
#%% CLEBSCH-GORDAN %%#
#   ==============   #
#
#CG_PATH = "/scratch/aitorcalvo/Bulegoa/ClebschGordan/Oh/cg_symbolic/"
#ASYM_PATH = "/scratch/aitorcalvo/Bulegoa/AntiSymmetricPart/Oh/";
CG_PATH = "/home/aitor/Bulegoa/ClebschGordan/Oh/cg_symbolic/"
ASYM_PATH = "/home/aitor/Bulegoa/AntiSymmetricPart/Oh/";

print( "Obtaining Clebsch-Gordan coefficients... " )
@time begin 
(oirreps,
 oirreps2indices,
 oirreps2dimensions,
 oindex2dimensions,
 cg_o_fullmatint) = get_cg_o_info( CG_PATH , ["A1g"]  )
cg_s_fullmatint = get_cg_s_fullmatint( max_spin2 );
end; println()


#   ==================================   #
#&& ONE-SHELL SYMSTATES AND MULTIPLETS &&#
#   ==================================   #

#   --------------------- #
#%% atomic symstate basis #
#   --------------------- #
identityrep = "A1g"
hiztegia = Dict( 
    "o1" => "A1g",
    "o2" => "A1g",
    "u" =>  0.5 ,
    "d" => -0.5
)

# orbital 1
o1 = (0,"o1")                                                                       
o1_up = ( o1... , 1 , "u" )                                                             
o1_do = ( o1... , 1 , "d" )                                                             
hilbert_o1 = HilbertSpace([ o1_up , o1_do ])
symstates_o1_nor = oneirrep_symstates( 
                    hilbert_o1 , 
                    hiztegia ,
                    identityrep ,
                    "$(ASYM_PATH)$(identityrep)_julia/" ) 
symstates_o1 = Dict( (q[1:5]...,1)=>s 
                    for (q,s) in symstates_o1_nor )
# orbital 2 (b)
o2 = (0,"o2")                                                                        
o2_up = ( o2... , 1 , "u" )                                                             
o2_do = ( o2... , 1 , "d" )                                                             
hilbert_o2 = HilbertSpace([ o2_up , o2_do ])
symstates_o2_nor = oneirrep_symstates( 
                    hilbert_o2 , 
                    hiztegia ,
                    identityrep ,
                    "$(ASYM_PATH)$(identityrep)_julia/" ) 
symstates_o2 = Dict( (q[1:5]...,2)=>s 
                    for (q,s) in symstates_o2_nor )

# mix 
symstates_0 = cg_reduce_product_states( symstates_o1 ,
                                         symstates_o2 ,
                                         CG_PATH )

# basis 
basis_0 = collect(values(symstates_0))[1].basis 

# multiplets 
multiplets_0 = get_multiplets( symstates_0 )

# printing
println( "BASIS FOR ATOMIC SHELL" )                                               
println( basis_0 )
println()
#println( "SYMSTATES" )
#print_symstates_Nordered( symstates_0 )
#println( "ATOMIC MULTIPLETS" )
#print_multiplets_Nordered( multiplets_0 )
#println()

##   ------------------ #
##%% selected symstates #
##   ------------------ #

# atomic symstates
symstates_atom = Dict{Tuple{Int64,String,Float64,Int64,Float64,Int64},State}()

# vacuum (seed) states
vac = State( basis_0.states[1] , basis_0 )

# creation operators 
cre_1u = Operator( SymbolCreationOperator(o1_up) , basis_0 )
cre_1d = Operator( SymbolCreationOperator(o1_do) , basis_0 )
cre_2u = Operator( SymbolCreationOperator(o2_up) , basis_0 )
cre_2d = Operator( SymbolCreationOperator(o2_do) , basis_0 )

# MAYBE IT IS MORE CORRECT TO REMOVE NORMALIZATIONS

# 0 triplet
i_0 = (4,"A1g",1.0)
m_0 = (i_0...,1)
# m=+1 
slat_uu = cre_1u*cre_2u*vac
s_0_p1 = 0.7*sqrt(2)*slat_uu
# m=0 
slat_ud = cre_1u*cre_2d*vac
slat_du = cre_1d*cre_2u*vac
s_0_0 = 0.7*( slat_ud + slat_du )
# m=-1 
slat_dd = cre_1d*cre_2d*vac
s_0_n1 = 0.7*sqrt(2)*slat_dd 
# symstates
symstates_m0 = Dict{Tuple{Int64,String,Float64,Int64,Float64,Int64},State}( 
                (i_0...,1, 1.0,1)=>s_0_p1 ,
                (i_0...,1, 0.0,1)=>s_0_0 ,
                (i_0...,1,-1.0,1)=>s_0_n1 )

# 1 singlet
i_1 = (4,"A1g",0.0)
m_1 = (i_1...,1)
# m=0 
slat_20 = cre_1u*cre_1d*vac 
slat_02 = cre_2u*cre_2d*vac
slat_du = cre_1d*cre_2u*vac
slat_ud = cre_1u*cre_2d*vac
s_1_0 = 0.7*(slat_20-slat_02) + 0.02*(slat_ud-slat_du)
# symstates
symstates_m1 = Dict{Tuple{Int64,String,Float64,Int64,Float64,Int64},State}( 
                (i_1...,1,0.0,1)=>s_1_0 )

# 0+ doublet
i_0p = (5,"A1g",0.5)
m_0p = (i_0p...,1)
# m=+1/2
slat_2u = cre_1u*cre_1d*cre_2u*vac
slat_u2 = cre_1u*cre_2u*cre_2d*vac
s_0p_p05 = -0.59*slat_2u - 0.51*slat_u2
# m=-1/2
slat_2d = cre_1u*cre_1d*cre_2d*vac
slat_d2 = cre_1d*cre_2u*cre_2d*vac
s_0p_n05 = -0.59*slat_2d - 0.51*slat_d2
# symstates
symstates_m0p = Dict{Tuple{Int64,String,Float64,Int64,Float64,Int64},State}( 
                 (i_0p...,1, 0.5,1)=>s_0p_p05,
                 (i_0p...,1,-0.5,1)=>s_0p_n05 )

## 1+ quadruplet
#i_1p = (5,"A1g",1.5)
#m_1p = (i_1p...,1)
## m=+3/2
#slat_2uuu = cre_1u*cre_1d*cre_2u*cre_3u*cre_4u*vac 
#slat_uu2u = cre_1u*cre_2u*cre_3u*cre_3d*cre_4u*vac
#slat_uuu2 = cre_1u*cre_2u*cre_3u*cre_4u*cre_4d*vac 
#s_1p_p15 = -0.99*slat_2uuu + 0.09*slat_uu2u + 0.05*slat_uuu2 
## m=+1/2 
#slat_2uud = cre_1u*cre_1d*cre_2u*cre_3u*cre_4d*vac
#slat_2udu = cre_1u*cre_1d*cre_2u*cre_3d*cre_4u*vac 
#slat_2duu = cre_1u*cre_1d*cre_2d*cre_3u*cre_4u*vac
#slat_uu2d = cre_1u*cre_2u*cre_3u*cre_3d*cre_4d*vac
#slat_ud2u = cre_1u*cre_2d*cre_3u*cre_3d*cre_4u*vac 
#slat_du2u = cre_1d*cre_2u*cre_3u*cre_3d*cre_4u*vac
#slat_uud2 = cre_1u*cre_2u*cre_3d*cre_4u*cre_4d*vac
#slat_udu2 = cre_1u*cre_2d*cre_3u*cre_4u*cre_4d*vac
#slat_duu2 = cre_1d*cre_2u*cre_3u*cre_4u*cre_4d*vac
#s_1p_p05 = -0.57*( slat_2uud + slat_2udu + slat_2duu ) +
#            0.05*( slat_uu2d + slat_ud2u + slat_du2u ) +
#            0.03*( slat_uud2 + slat_udu2 + slat_duu2 )
## m=-1/2 
#slat_2ddu = cre_1u*cre_1d*cre_2d*cre_3d*cre_4u*vac
#slat_2dud = cre_1u*cre_1d*cre_2d*cre_3u*cre_4d*vac 
#slat_2udd = cre_1u*cre_1d*cre_2u*cre_3d*cre_4d*vac
#slat_dd2u = cre_1d*cre_2d*cre_3u*cre_3d*cre_4u*vac
#slat_du2d = cre_1d*cre_2u*cre_3u*cre_3d*cre_4d*vac 
#slat_ud2d = cre_1u*cre_2d*cre_3u*cre_3d*cre_4d*vac
#slat_ddu2 = cre_1d*cre_2d*cre_3u*cre_4u*cre_4d*vac
#slat_dud2 = cre_1d*cre_2u*cre_3d*cre_4u*cre_4d*vac
#slat_udd2 = cre_1u*cre_2d*cre_3d*cre_4u*cre_4d*vac
#s_1p_n05 = -0.57*( slat_2ddu + slat_2dud + slat_2udd ) +
#            0.05*( slat_dd2u + slat_du2d + slat_ud2d ) +
#            0.03*( slat_ddu2 + slat_dud2 + slat_udd2 )
## m=-3/2
#slat_2ddd = cre_1u*cre_1d*cre_2d*cre_3d*cre_4d*vac 
#slat_dd2d = cre_1d*cre_2d*cre_3u*cre_3d*cre_4d*vac
#slat_ddd2 = cre_1d*cre_2d*cre_3d*cre_4u*cre_4d*vac 
#s_1p_n15 = -0.99*slat_2ddd + 0.09*slat_dd2d + 0.05*slat_ddd2 
## symstates 
#symstates_m1p = Dict{Tuple{Int64,String,Float64,Int64,Float64,Int64},State}( 
#                 (i_1p...,1,+1.5,1)=>s_1p_p15,
#                 (i_1p...,1,+0.5,1)=>s_1p_p05,
#                 (i_1p...,1,-0.5,1)=>s_1p_n05,
#                 (i_1p...,1,-1.5,1)=>s_1p_n15 )

# 2+ doublet 
i_2p = (5,"A1g",0.5) 
m_2p = (i_2p...,2)
# m=+1/2
slat_u2 = cre_1u*cre_2u*cre_2d*vac
slat_2u = cre_1u*cre_1d*cre_2u*vac
s_2p_p05 = -0.61*slat_u2 + 0.43*slat_2u
# m=-1/2
slat_d2 = cre_1d*cre_2u*cre_2d*vac
slat_2d = cre_1u*cre_1d*cre_2d*vac
s_2p_n05 = -0.61*slat_d2 + 0.43*slat_2d
# symstates 
symstates_m2p = Dict{Tuple{Int64,String,Float64,Int64,Float64,Int64},State}( 
                 (i_2p...,1,+0.5,2)=>s_2p_p05,
                 (i_2p...,1,-0.5,2)=>s_2p_n05 )

# 0- doublet
i_0n = (3,"A1g",0.5)
m_0n = (i_0n...,1)
# m=+1/2
slat_0u = cre_2u*vac
slat_u0 = cre_1u*vac
s_0n_p05 = -0.98*slat_0u + 0.1*slat_u0
# m=-1/2
slat_0d = cre_2d*vac 
slat_d0 = cre_1d*vac
s_0n_n05 = -0.99*slat_0d + 0.1*slat_d0
# symstates 
symstates_m0n = Dict{Tuple{Int64,String,Float64,Int64,Float64,Int64},State}( 
                 (i_0n...,1,+0.5,1)=>s_0n_p05,
                 (i_0n...,1,-0.5,1)=>s_0n_n05 )


# 1- doublet
i_1n = (3,"A1g",0.5)
m_1n = (i_0n...,2)
# m=+1/2
s_1n_p05 = -0.99*slat_u0 + 0.1*slat_0u
# m=-1/2
s_1n_n05 = -0.99*slat_d0 + 0.1*slat_0d
# symstates 
symstates_m1n = Dict{Tuple{Int64,String,Float64,Int64,Float64,Int64},State}( 
                 (i_1n...,1,+0.5,2)=>s_1n_p05,
                 (i_1n...,1,-0.5,2)=>s_1n_n05 )

symstates_clean = Dict{Tuple{Int64,String,Float64,Int64,Float64,Int64},State}(
                    (0,"A1g",0.0,1,0.0,1)=>vac ) 
symstates_imp = merge( symstates_m0  , 
                       symstates_m1  ,
                       symstates_m0n ,
                       symstates_m1n ,
                       symstates_m0p ,
                       #symstates_m1p ,
                       symstates_m2p )
multiplets_imp = get_multiplets(symstates_imp)

symstates_atom = calculation=="IMP" ? symstates_imp : symstates_clean

# atomic symstates and multiplets
multiplets_atom = get_multiplets( symstates_atom )

println( "***********************************************************************")
println( "ATOMIC SELECTED MULTIPLETS" )
print_multiplets_Nordered( multiplets_atom )
#println( "ATOMIC SELECTED SYMSTATES" )
#print_symstates_Nordered( symstates_atom )
#println( "***********************************************************************")
println()

##   =================== #
##&& INITIAL HAMILTONIAN #
##   =================== #

print( "Preparing initial Hamiltonian... " )
@time begin

# parameter renormalization
function rescale( 
            num::N ,
            L::Float64 ,
            z::Float64 ,
            discretization::String ) where {N<:Number}
    a::Float64 = compute_ebar0_z(
                     z,
                     L;
                     discretization=discretization)
    return num/(a*sqrt(L))
end
α = compute_ebar0_z( z,
                     L;
                     discretization=discretization)
E_0  = rescale(E_0,  L, z, discretization) 
E_1  = rescale(E_1,  L, z, discretization)
E_0p = rescale(E_0p, L, z, discretization)
E_1p = rescale(E_1p, L, z, discretization)
E_2p = rescale(E_2p, L, z, discretization)
E_0n = rescale(E_0n, L, z, discretization)
E_1n = rescale(E_1n, L, z, discretization)
gam  = rescale(gam,  L, z, discretization)


#   ----- #
#%% irreu #
#   ----- #
irrEU_clean = get_irrEU_clean( "A1g" )
irrEU_imp   = Dict{Tuple{Int64, String, Float64}, Tuple{Vector{Float64}, Matrix{ComplexF64}}}( 
                   i_0  => ([E_0], [ComplexF64(1.0)][:,:]) ,
                   i_1  => ([E_1], [ComplexF64(1.0)][:,:]) ,
                   i_0p => ([E_0p,E_2p],ComplexF64.(diagm([1.0,1.0]))) ,#i_0p=i_2p
                   #i_1p => ([E_1p],[ComplexF64(1.0)][:,:]) ,
                   i_0n => ([E_0n,E_1n],ComplexF64.(diagm([1.0,1.0]))) )#i_0n=i_1n
irrEU = calculation=="IMP" ? irrEU_imp : irrEU_clean
multiplets_block = calculation=="IMP" ? multiplets_atom :
                                        Set([ (0,"A1g",0.0,1) ])    
end #timing
println()

print_spectrum( irrEU )


##   ================================================   #
##&& NUMERICAL ADDITION OF INNERMOST CONDUCTION SHELL &&#
##   ================================================   #

#   ---------------   #
#%% combinations u' %%#
#   ---------------   #
print( "Computing combinations_uprima... " )
@time begin

combinations_uprima = 
    Dict{ Tuple{Int64,Int64,Int64,Int64} , NTuple{2,Tuple{Int64,Int64,Int64,Int64}} }()
m_vac = (0,oirreps2indices["A1g"],0,1)
if calculation=="IMP"
    for m_mu in multiplets_atom
        mint_mu = convert_to_int( m_mu , oirreps2indices )
        push!( combinations_uprima , mint_mu=>(mint_mu,m_vac) )
    end 
elseif calculation=="CLEAN" 
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
end #timing
println()

#   ---------------------------   #
#%% atom pseudo-CG coefficients %%#
#   ---------------------------   #

print( "Constructing atomic pcg... " )
@time begin

# hoppers
pcg_atom = Dict{NTuple{3,Tuple{Int64,String,Float64,Int64,Float64,Int64}},ComplexF64}()
q_1up = (1,"A1g",0.5,1, 0.5,1)
q_1do = (1,"A1g",0.5,1,-0.5,1)
q_2up = (1,"A1g",0.5,1, 0.5,2)
q_2do = (1,"A1g",0.5,1,-0.5,2)
qq_a  = [q_1up,q_1do,q_2up,q_2do] 

# pcg
for (q1,s1) in symstates_atom,
    (q2,s2) in symstates_atom

    cup_1 = s1 * cre_1u * s2 
    cdo_1 = s1 * cre_1d * s2

    cup_2 = s1 * cre_2u * s2 
    cdo_2 = s1 * cre_2d * s2

    isapprox(cup_1,0.0) || (pcg_atom[(q1,q_1up,q2)]=cup_1)
    isapprox(cdo_1,0.0) || (pcg_atom[(q1,q_1do,q2)]=cdo_1)

    isapprox(cup_2,0.0) || (pcg_atom[(q1,q_2up,q2)]=cup_2)
    isapprox(cdo_2,0.0) || (pcg_atom[(q1,q_2do,q2)]=cdo_2)

end

end #timing
println()


#   ----------------------------------- #
#%% shell symstates, multiplets and pcg #
#   ----------------------------------- #

print( "Obtaining shell info... " )
@time begin

# symstates
symstates_shell = symstates_0 

# multiplets
multiplets_shell = get_multiplets( symstates_shell )
irrmult_shell = get_irreps( multiplets_shell ; multiplicity=true )

# pcg
pcg_shell = Dict{NTuple{3,Tuple{Int64,String,Float64,Int64,Float64,Int64}},ComplexF64}()
for (q1,s1) in symstates_shell,
    (q2,s2) in symstates_shell

    cup_1 = s1 * cre_1u * s2 
    cdo_1 = s1 * cre_1d * s2

    cup_2 = s1 * cre_2u * s2 
    cdo_2 = s1 * cre_2d * s2

    isapprox(cup_1,0.0) || (pcg_shell[(q1,q_1up,q2)]=cup_1)
    isapprox(cdo_1,0.0) || (pcg_shell[(q1,q_1do,q2)]=cdo_1)

    isapprox(cup_2,0.0) || (pcg_shell[(q1,q_2up,q2)]=cup_2)
    isapprox(cdo_2,0.0) || (pcg_shell[(q1,q_2do,q2)]=cdo_2)
end

end #timing 
println()


#   -----------------   #
#%% hopping parameter %%#
#   -----------------   #
hop = Dict( 
    (oirreps2indices["A1g"],1) => ComplexF64(gam) ,
    (oirreps2indices["A1g"],2) => ComplexF64(gam) 
)
# hopchannels 
hopchannels = collect(keys( hop ))


#   ------------------------ #
##% conversion to int format #
#   ------------------------ #
multiplets_block = Set([ convert_to_int(m,oirreps2indices) 
                         for m in multiplets_block ])
multiplets_imp = Set([ convert_to_int(m,oirreps2indices) 
                         for m in multiplets_imp ])
multiplets_shell = Set([ convert_to_int(m,oirreps2indices) 
                         for m in multiplets_shell ])
irrEU = Dict( convert_to_int(G,oirreps2indices)=>x 
                         for (G,x) in irrEU )
qq_a = [ convert_to_int(q,oirreps2indices) 
                         for q in qq_a ]
pcg_atom  = Dict( (convert_to_int(k[1],oirreps2indices),
                   convert_to_int(k[2],oirreps2indices),
                   convert_to_int(k[3],oirreps2indices))=>v 
                 for (k,v) in pcg_atom )
pcg_shell = Dict( (convert_to_int(k[1],oirreps2indices),
                   convert_to_int(k[2],oirreps2indices),
                   convert_to_int(k[3],oirreps2indices))=>v 
                 for (k,v) in pcg_shell )

#   --------------------   #
#%% reduced pcg matrices %%#
#   --------------------   #

print( "Computing reduced pcg... " )
@time begin

multiplets_a = Set( (q[1:3]...,q[6]) for q in qq_a )
pcgred_atom = get_redmat2(
                pcg_atom ,
                multiplets_block ,
                multiplets_a ,
                cg_o_fullmatint ,
                cg_s_fullmatint ;
                verbose=false )
pcgred_shell = get_redmat2( 
                pcg_shell ,
                multiplets_shell ,
                multiplets_a ,
                cg_o_fullmatint ,
                cg_s_fullmatint ;
                verbose=false )

end
println()
println( "PCGRED ATOM" )
print_dict( pcgred_atom )
println( "PCGRED SHELL" )
print_dict( pcgred_shell )
println( "==================" )
println( "ATOMIC EXCITATIONS" )
println( "==================" )
for ((G_i,G_a,G_j),mat) in pcgred_atom
    (G_i==(4,1,2) || G_j==(4,1,2)) || continue
    println("askldjf")
    for rrr in CartesianIndices(mat)
        r_i,r_a,r_j = Tuple(rrr)
        matel = mat[r_i,r_a,r_j]
        if !isapprox(matel,0.0im)

            m_i = (G_i...,r_i) 
            m_j = (G_j...,r_j)
            m_a = (G_a...,r_a)

            @show m_i,m_a,m_j
            @show matel 
            println()


        end
    end

end
println()


#   -----------------   #
#%% excitation matrix %%#
#   -----------------   #

print( "Computing excitation matrix... ") 
@time begin

## selected MO
#atom_creops = Dict(
#    (1,"A1g",0.5,1, 0.5,1) => cre_3u,
#    (1,"A1g",0.5,1,-0.5,1) => cre_3d)
# all MOs
atom_creops = Dict(
    (1,"A1g",0.5,1, 0.5,1) => cre_1u, 
    (1,"A1g",0.5,1,-0.5,1) => cre_1d,
    (1,"A1g",0.5,1, 0.5,2) => cre_2u,
    (1,"A1g",0.5,1,-0.5,2) => cre_2d)

#pcg_atom_full = Dict()
#for (qbra,sbra) in symstates_atom,
#    (qket,sket) in symstates_atom,
#    (qa,oa)     in atom_creops 
#
#    c = sbra*oa*sket 
#    if !isapprox( abs2(c) , 0.0 )
#        pcg_atom_full[(qbra,qa,qket)] = c 
#    end
#
#end
#pcg_atom_full = Dict( 
#    (convert_to_int(k[1],oirreps2indices),
#     convert_to_int(k[2],oirreps2indices),
#     convert_to_int(k[3],oirreps2indices))=>v
#     for (k,v) in pcg_atom_full )
pcg_atom_full = pcg_atom
#multiplets_a = Set((q[1:3]...,q[6]) for (s1,q,s2) in keys(pcg_atom_full))
multiplets_atomhop = Set((convert_to_int(q,oirreps2indices)[1:3]...,q[6]) 
                   for q in keys(atom_creops))
Mred, AA = setup_redmat_AA_orbitalresolved(
            pcg_atom_full ,
            multiplets_block ,
            multiplets_atomhop ,
            cg_o_fullmatint ,
            cg_s_fullmatint ,
            irrEU ;
            verbose=false )

end
println()

#println( "MRED" )
#print_dict( Mred )
#println()

#   ------------------- #
#%% impurity multiplets #
#   ------------------- #
print( "Impurity multiplet setup... " )
@time begin 
omults = ordered_multiplets(multiplets_block)
mult2index = Dict( m=>i for (i,m) in 
                   enumerate(omults))
mm_i,m_imp = setup_impmultinfo( 
                multiplets_block ,
                irrEU ,
                betabar ,
                oindex2dimensions )
end #timing
println()

#   -------------------   #
#%% clebsch-gordan sums %%#
#   -------------------   #
print( "Precomputing Clebsch-Gordan sums..." )
@time begin
Bsum_o_dict,Bsum_s_dict,Csum_o_dict,Csum_s_dict =
    precompute_CGsums(
            oirreps ,
            multiplets_a,
            multiplets_imp ,
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
end #timing 
println()


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

Mred, AA = update_redmat_AA_CGsummethod_orbitalresolved(
            Mred ,
            irrEU ,
            combinations_uprima ,
            collect(multiplets_atomhop) ,
            cg_o_fullmatint ,
            cg_s_fullmatint ,
            Karray_orbital ,
            Karray_spin ,
            AA ;
            verbose=false )
end #timing
println()
println( "SPECTRUM AFTER ADDING FIRST SHELL" )
print_spectrum(irrEU)

#   =============== #
#%% NRG CALCULATION #
#   =============== #

@time begin
if run=="spectral"
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
                      alpha=Float64(α) ,
                      multiplets_atomhop=collect(multiplets_atomhop) ,
                      orbitalresolved=true)
elseif run=="thermo" 
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
                      spectral=false ,
                      alpha=Float64(α) , )
end
end #timing


#   ===========   #
#%% PERFORMANCE %%#
#   ===========   #

println( "===========" )
println( "PERFORMANCE" ) println( "===========" )
print_performance_onestep(nrg)
println()

#   ==============   #
#%% SAVING TO FILE %%#
#   ==============   #

if run=="thermo" 
    
    # impurity properties 
    if calculation=="IMP" 
        write_impurity_info( nrg , omults , mult2index , orbital , z )
    end
    
    # thermodata for this given value of z
    write_thermodata_onez( nrg , calculation , orbital , z )
    
    # thermo diff
    if calculation=="IMP"
        write_thermodiff( orbital , z )
    end

elseif run=="spectral" 

    for (m_a,spectralo) in nrg.specfunc
        writedlm( "spectral/spectral_$(m_a[end]).dat" , spectralo )
    end

end
