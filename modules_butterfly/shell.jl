include( "symbols.jl" )
include( "numericals.jl" )
include( "symmetry.jl" )
include( "discretization.jl" )

using StaticArrays
using ProgressMeter


# ###########################################
# SHELL CREATION OPERATORS
#
# - Used for computing pseudo-CG coefficients
#
# ...........................................


function basis2coperators( basis::CanonicalBasis ; n=-1 )
    # given a canonical basis, it returns an array of 
    # creation operators, one for each element of the basis.
    #
    # input:
    #
    # - basis : canonical basis of the system 
    # - n : number of states in the basis --if n==-1, it computes
    #       all the elements, otherwise it selects the n 
    #       subspace.
    #
    # output:
    # 
    # - coperators : array of creation operators for the basis.
    #
    coperators::Vector{Operator} = []
    for sfs::SymbolFockState in basis.states
        if n==-1
            push!( coperators , sfs2coperator(sfs,basis) )
        else
            if sum(sfs.occ)==n
                push!( coperators , sfs2coperator(sfs,basis) )
            end
        end
    end
    return coperators 
end

function shell_coperators( shell_basis::CanonicalBasis , hiztegia )
    # same as basis2coperators (above) but stores the creation
    # operators as a dictionary.
    #
    # input
    #
    # - shell_basis : canonical basis of shell.
    # - hiztegia : irrep pseudoname => irrep name
    #
    # output
    #
    # - coperator : sym => creation operatorcoperator : sym => creation operator
    #
    coperators = Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , Operator }()
    @inbounds for (i,sfs) in enumerate(shell_basis.states)
        sum(sfs.occ)!==1 && continue
        idx = findfirst( sfs.occ )
        tup = sfs.hilbert.states[idx] 
        # sym = ( N , I , S , mu , m , r )
        sym = ( convert(Int64,1) , 
                convert(String,hiztegia[tup[2]]) , 
                convert(Float64,0.5) , 
                convert(Int64,tup[3]) , 
                convert(Float64,hiztegia[tup[4]]) , 
                convert(Int64,1) )
        cop = sfs2coperator(sfs,shell_basis)
        sym2cop = Dict( sym => cop )
        merge!( coperators , sym2cop )
    end
    return coperators 
end


# ###########################################
# PSEUDO--CLEBSCH-GORDAN COEFFICIENTS
#
# - Matrix elements of the creation operators 
#   with states from a single shell:
#           ( q_nu | c^\dagger_q_a | q_mu )
#
# ............................................

# int format method
function get_pseudoCG( 
                symstates_shell::Dict , 
                basis_shell::CanonicalBasis , 
                hiztegia::Dict , 
                oirreps2dimensions::Dict{String,Int64} ; 
                verbose=false )
    # compute the pseudo-CG coefficients for the shells.
    #
    # input
    #
    # - symstates_shell : all symstates of a single shell (any)
    # - basis_shell : basis for only that shell 
    # - hiztegia : dict( irrep name => standard notation )
    # - [format : "int" or "standard"]
    #
    # output
    #
    # - pseudo-CG coefficients as dict( (q_nu,q_a,q_mu)=>coeff )
    
    pseudoCG = Dict()

    # shell creation operators
    shell_cops = shell_coperators( basis_shell , hiztegia )

    for (q_nu,s_nu) in symstates_shell, 
        (q_mu,s_mu) in symstates_shell

        for (q_a,c_a) in shell_cops 
            q = ( q_nu , q_a , q_mu ) 
            # ( q_nu | c^\dagger_{q_a} | q_mu )
            # notice that it is not the complex conjugate!!!
            coeff = s_nu * c_a * s_mu 
            isapprox( abs(coeff) , 0.0 , atol=1e-6 ) || push!( pseudoCG , q=>coeff ) 
        end

    end

    pseudoCG = Dict( (convert_to_int(k[1],oirreps2dimensions),
                      convert_to_int(k[2],oirreps2dimensions),
                      convert_to_int(k[3],oirreps2dimensions))=>v 
                      for (k,v) in pseudoCG )

    if verbose 
        for (q,coeff) in pseudoCG 
            println( "q = $q" )
            println( "coeff = $coeff" )
            println()
        end
    end

    return pseudoCG 
end

function get_pseudoCG_mat( 
            pseudoCG::Dict{Tuple{NTuple{6, Int64}, NTuple{6, Int64}, NTuple{6, Int64}}, ComplexF64} ,
            irrmult_0::Set{Tuple{Tuple{Int64, String, Float64}, Int64}} ,
            oindex2dimensions::Vector{Int64} )
    #
    # generate pseudo-CG coefficients in matrix format
    #
    # input
    #
    # - pseudoCG : coefficients in string format 
    # - irrmult_0 : irrep=>multiplicity for the atomic shell
    # - oindex2dimensions: oindex=>dimension 
    #
    # output
    #
    # - pcgmat : coefficients in matrix format 
    #
    
    # MAYBE IT CAN BE RESTRICTED TO IRREP COMBINATIONS WHERE
    # I3 IS CONTAINED IN I1xI2
    
    pcgmat = Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,9} }()

    # irrmult in int format (isn't it already?)
    shell_irrmult = Set( (convert_to_int(G,oirreps2indices),R) 
                        for (G,R) in irrmult_0 )

    for (G_mu,R_mu) in shell_irrmult,
        (G_a ,R_a ) in shell_irrmult,
        (G_nu,R_nu) in shell_irrmult 

        (N_mu,I_mu,S_mu) = G_mu 
        (N_a ,I_a ,S_a ) = G_a 
        (N_nu,I_nu,S_nu) = G_nu

        D_Imu = oindex2dimensions[I_mu]
        D_Inu = oindex2dimensions[I_nu]
        D_Ia  = oindex2dimensions[I_a]
        D_Smu = S_mu+1
        D_Snu = S_nu+1
        D_Sa  = S_a+1

        mat = zeros( ComplexF64 , D_Imu, D_Smu , D_Inu , D_Snu ,
                                  R_mu , R_nu ,
                                  D_Ia , D_Sa , R_a ) 

        @inbounds for i_mu in 1:D_Imu,
                      s_mu in 1:D_Smu,
                      r_mu in 1:R_mu,
                      i_a  in 1:D_Ia,
                      s_a  in 1:D_Sa,
                      r_a  in 1:R_a,
                      i_nu in 1:D_Inu,
                      s_nu in 1:D_Snu,
                      r_nu in 1:R_nu

            q_mu = (N_mu,I_mu,S_mu,i_mu,(2*(s_mu-1)-S_mu),r_mu)
            q_a  = (N_a, I_a, S_a, i_a, (2*(s_a -1)-S_a ), r_a)
            q_nu = (N_nu,I_nu,S_nu,i_nu,(2*(s_nu-1)-S_nu),r_nu)

            # ordering for optimal pseudoCG_up_vp
            mat[ i_mu , s_mu , i_nu , s_nu , r_mu , r_nu , i_a , s_a , r_a ] =
                    get( pseudoCG , (q_mu,q_a,q_nu) , zero(ComplexF64) )
        end

        pcgmat[G_mu,G_a,G_nu] = mat
    end
    return pcgmat
end

# standard format method
function get_pseudoCG( 
            symstates_shell::Dict , 
            basis_shell::CanonicalBasis , 
            hiztegia::Dict ; 
            verbose=false )
    #
    # compute the pseudo-CG coefficients for the shells.
    #
    # input:
    #
    # - symstates_shell : all symstates of a single shell (any)
    # - basis_shell : basis for only that shell 
    # - hiztegia : dict( irrep name => standard notation )
    # - [format : "int" or "standard"]
    #
    # output:
    #
    # - pseudo-CG coefficients as dict( (q_nu,q_a,q_mu)=>coeff )
    #
    pseudoCG = Dict()

    # shell creation operators
    shell_cops = shell_coperators( basis_shell , hiztegia )

    for (q_nu,s_nu) in symstates_shell, (q_mu,s_mu) in symstates_shell
        for (q_a,c_a) in shell_cops 
            q = ( q_nu , q_a , q_mu ) 
            # ( q_nu | c^\dagger_{q_a} | q_mu )
            # notice that it is not the complex conjugate!!!
            coeff = s_nu * c_a * s_mu 
            isapprox( coeff , 0 ) || push!( pseudoCG , q=>coeff ) 
        end
    end

    if verbose 
        for (q,coeff) in pseudoCG 
            println( "q = $q" )
            println( "coeff = $coeff" )
            println()
        end
    end
    return pseudoCG 
end


# ################################
# MULTIPLETS u IN COMBINATION i,mu
# ................................
# string method
function get_combination_multiplets( 
            multiplets_block::Set{Tuple{Int64,String,Float64,Int64}} , 
            multiplets_shell::Set{Tuple{Int64,String,Float64,Int64}} ,
            cg_o_fullmat::Dict{ NTuple{3,Tuple{Int64,String,Float64}} , Matrix{ComplexF64} } ;
            verbose=false )

    mui2u = Dict{NTuple{2,Tuple{Int64,String,Float64,Int64}},Vector{Tuple{Int64,String,Float64,Int64}}}()
    cg_mui_s::Dict{NTuple{6, Float64}, ComplexF64} = Dict()

    G2R::Dict{Tuple{Int64,String,Float},Int64} = Dict()

    for m_i::Tuple{Int64,String,Float64,Int64} in multiplets_block, 
        m_mu::Tuple{Int64,String,Float64,Int64} in multiplets_shell 
        if verbose 
            println( "m_mu = $m_mu ; m_i = $m_i" )
            println( "============================" )
        end
        
        (N_i,I_i,S_i,r_i) = m_i
        (N_mu,I_mu,S_mu,r_mu) = m_mu 

        cg_mui_o = Dict( k=>v for (k,v) in cg_o_fullmat
                         if (k[1],k[2])==(I_mu,I_i) )
        cg_mui_s = cg_spin(    S_mu , S_i )

        N_u = N_i + N_mu
        II_u = Set( k[3] for k in keys(cg_mui_o) )::Set{String}
        SS_u = Set( k[5] for k in keys(cg_mui_s) )::Set{Float64}

        for I_u::String in II_u, S_u::Float64 in SS_u
            G_u = (N_u,I_u,S_u)
            (G_u in keys(G2R)) ? G2R[G_u]+=1 : merge!(G2R,Dict(G_u=>1))
            m_u = ( G_u... , G2R[G_u] )
            mergewith!( (x,y)->append!(x,y) , mui2u , Dict( (m_mu,m_i)=>[m_u]) )
            verbose && println( m_u )
        end
        verbose && println()
    end
    return mui2u 
end

# method for int 
function get_combination_multiplets( 
            multiplets_block::Set{NTuple{4,Int64}} , 
            multiplets_shell::Set{NTuple{4,Int64}} ,
            cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} } ;
            verbose=false )

    mui2u = Dict{NTuple{2,NTuple{4,Int64}},Vector{NTuple{4,Int64}}}()
    cg_mui_s::Dict{NTuple{6,Float64}, ComplexF64} = Dict()

    G2R = Dict{Tuple{Int64,Int64,Int64},Int64}()

    for m_i::NTuple{4,Int64} in multiplets_block, 
        m_mu::NTuple{4,Int64} in multiplets_shell 
        if verbose 
            println( "m_mu = $m_mu ; m_i = $m_i" )
            println( "============================" )
        end
        
        (N_i,I_i,S_i,r_i) = m_i
        (N_mu,I_mu,S_mu,r_mu) = m_mu 

        cg_mui_o = Dict( k=>v for (k,v) in cg_o_fullmatint
                         if (k[1],k[2])==(I_mu,I_i) )
        cg_mui_s = cg_spin(    S_mu/2.0 , S_i/2.0 )

        N_u = N_i + N_mu
        II_u = Set( k[3] for k in keys(cg_mui_o) )::Set{Int64}
        SS_u = Set( convert(Int64,2*k[5]) for k in keys(cg_mui_s) )::Set{Int64}

        for I_u::Int64 in II_u, S_u::Int64 in SS_u
            G_u = (N_u,I_u,S_u)
            (G_u in keys(G2R)) ? G2R[G_u]+=1 : merge!(G2R,Dict(G_u=>1))
            m_u = ( G_u... , G2R[G_u] )
            mergewith!( (x,y)->append!(x,y) , mui2u , Dict( (m_mu,m_i)=>[m_u]) )
            verbose && println( m_u )
        end
        verbose && println()
    end
    return mui2u 
end


# #############################################
# PSEUDO-CG COEFFICIENTS TRANSFORMED 
# BY DIAGONALIZATION MATRICES U
# 
# - Coefficient 
#      ( i=u'_tilde | c^\dagger_a | j=v'_tilde )
#   obtained by applying transfromation U to 
#   (n-1) in-shell pseudo-CG coefficients
# - It is computed in-place in order to avoid 
#   having to store large amounts of data,
#   since only the universal P-CG and the 
#   regular CG are needed.
# ..............................................

function get_pseudoCG_up_vp( 
        u_mu_i::Tuple{NTuple{6,Int64},NTuple{4,Int64},NTuple{4,Int64}} , 
        v_nu_j::Tuple{NTuple{6,Int64},NTuple{4,Int64},NTuple{4,Int64}} , 
        q_a::NTuple{6,Int64} , 
        oindex2dimensions::Vector{Int64} ,
        pcg::Dict{ NTuple{3,NTuple{6,Int64}} , ComplexF64 } , 
        pcgmat::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,9} } ,
        cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} } ;
        verbose=false )
    # compute the clebsch-gordan series of pseudo-CG coefficients 
    # of the block states before diagonalization (hence the prima,
    # and not the tilde). every letter (u,v,i,j,mu,nu) in the function 
    # refers to the corresponding primed index.
    # input:
    # - u_mu_i : (q_u,m_mu,m_i)
    # - v_nu_j : (q_v,m_nu,m_j)
    # - q_a  
    # - oirreps2dimensions : Dict( irrep => dimension )
    # - pcg : universal pseudo-CG coefficients 
    #          ( q_nu | c^\dagger_{q_a} | q_mu )
    # - pcgmat : the same but in the format 
    #          ( G_nu , G_a , G_mu ) => coefficient matrix
    # output:
    # - coefficient
    #          ( q_u' | c^\dagger | q_v')
    
    # type annotations
    coeff::ComplexF64 = 0.0im
    cg_muiu::ComplexF64 = 0.0im
    cg_nujv::ComplexF64 = 0.0im
    cg_total::ComplexF64 = 0.0im

    # we need all quantum numbers for u,v
    # and the multiplets for mu,nu,i,j
    (q_u,m_mu,m_i) = u_mu_i
    (q_v,m_nu,m_j) = v_nu_j

    # Informative but expensive. Uncomment for testing.
    #if verbose 
    #    println( "COMPUTING (N-1) PSG" )
    #    println( "q_u' = $q_u ; q_v' = $q_v" )
    #    println( "m_i' = $m_i ; m_j' = $m_j" )
    #    println( "m_mu' = $m_mu ; m_nu' = $m_nu" ) 
    #    println( "q_a = $q_a" )
    #end

    # quantum numbers 
    # q_u, q_v 
    (N_u,I_u,S_u,i_u,s_u,r_u) = q_u
    (N_v,I_v,S_v,i_v,s_v,r_v) = q_v
    # m_mu, m_nu 
    (N_mu,I_mu,S_mu,r_mu) = m_mu
    (N_nu,I_nu,S_nu,r_nu) = m_nu
    # m_i, m_j
    (N_i,I_i,S_i,r_i) = (N_ij,I_ij,S_ij,r_ij) = m_i
    (N_j,I_j,S_j,r_j) = m_j

    # pseudo-cg matrix 
    G_mu = (N_mu,I_mu,S_mu)
    G_nu = (N_nu,I_nu,S_nu)
    G_a  = (N_a,I_a,S_a) = (q_a[1],q_a[2],q_a[3])
    (i_a,si_a,r_a) = (q_a[4],Int64((q_a[5]+S_a)/2+1),q_a[6]) 
    ((G_mu,G_a,G_nu) in keys(pcgmat)) || (return zero(ComplexF64))
    pcg_local = @view pcgmat[G_mu,G_a,G_nu][:,:,:,:,r_mu,r_nu,i_a,si_a,r_a]

    # clebsch-gordan matrices
    !( ((I_mu,I_ij,I_u) in keys(cg_o_fullmatint)) &&
       ((I_nu,I_ij,I_v) in keys(cg_o_fullmatint)) ) && (return zero(ComplexF64))
    @inbounds begin
        cgomat_muiu = @view cg_o_fullmatint[I_mu,I_ij,I_u][:,:,i_u]
        cgomat_nujv = @view cg_o_fullmatint[I_nu,I_ij,I_v][:,:,i_v]
    end

    # dimension of orbital irreps -> info about partners
    d_Iij = oindex2dimensions[I_i] # I_i = I_j
    d_Imu = oindex2dimensions[I_mu]
    d_Inu = oindex2dimensions[I_nu]

    # main iteration 
    @inbounds for s_ij in (-S_ij):2:S_ij,
                  i_ij in 1:d_Iij,
                  s_nu in (-S_nu):2:S_nu,
                  i_nu in 1:d_Inu,
                  s_mu in (-S_mu):2:S_mu,
                  i_mu in 1:d_Imu

        i_i = i_j = i_ij 
        s_i = s_j = s_ij        

        if verbose 
            println( "g_mu' = ( $i_mu , $(s_mu/2.0) )" )
            println( "g_nu' = ( $i_nu , $(s_nu/2.0) )" )
        end

        # coefficient
        si_mu = Int64((s_mu+S_mu)/2+1)
        si_nu = Int64((s_nu+S_nu)/2+1)
        pcgcoeff = pcg_local[i_mu,si_mu,i_nu,si_nu]
        #q_mu = (N_mu,I_mu,S_mu,i_mu,s_mu,r_mu)
        #q_nu = (N_nu,I_nu,S_nu,i_nu,s_nu,r_nu)
        #pcgcoeff = get( pcg , (q_mu,q_a,q_nu) , zero(ComplexF64) )
        pcgcoeff==zero(pcgcoeff) && continue
        verbose && println( "pcgcoeff = $pcgcoeff")

        # cg coefficients
        cg_muiu = cgomat_muiu[i_mu,i_ij] * 
                PartialWaveFunctions.CG_doublearg(S_mu,s_mu,S_i,s_i,S_u,s_u)
        cg_nujv = cgomat_nujv[i_nu,i_ij] * 
                PartialWaveFunctions.CG_doublearg(S_nu,s_nu,S_j,s_j,S_v,s_v) 
        cg_total = conj( cg_muiu ) * cg_nujv 

        if verbose 
            println( "cg_muiu = $cg_muiu" )
            println( "cg_nujv = $cg_nujv" )
            println( "cg_total = $cg_total" )
        end

        # putting all together
        coeff += cg_total * pcgcoeff
        verbose && println( "coeff = $coeff" )
    end
    verbose && println( "final coeff = $coeff" )
    verbose && println()
    return coeff 
end

function get_pseudoCG_up_vp_spinarray( 
        u_mu_i::Tuple{NTuple{6,Int64},NTuple{4,Int64},NTuple{4,Int64}} , 
        v_nu_j::Tuple{NTuple{6,Int64},NTuple{4,Int64},NTuple{4,Int64}} , 
        q_a::NTuple{6,Int64} , 
        oindex2dimensions::Vector{Int64} ,
        pcg::Dict{ NTuple{3,NTuple{6,Int64}} , ComplexF64 } , 
        pcgmat::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,9} } ,
        cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} } ,
        cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} } ;
        verbose=false )
    # compute the clebsch-gordan series of pseudo-CG coefficients 
    # of the block states before diagonalization (hence the prima,
    # and not the tilde). every letter (u,v,i,j,mu,nu) in the function 
    # refers to the corresponding primed index.
    # input:
    # - u_mu_i : (q_u,m_mu,m_i)
    # - v_nu_j : (q_v,m_nu,m_j)
    # - q_a  
    # - oirreps2dimensions : Dict( irrep => dimension )
    # - pcg : universal pseudo-CG coefficients 
    #          ( q_nu | c^\dagger_{q_a} | q_mu )
    # - pcgmat : the same but in the format 
    #          ( G_nu , G_a , G_mu ) => coefficient matrix
    # output:
    # - coefficient
    #          ( q_u' | c^\dagger | q_v')
    
    # type annotations
    coeff::ComplexF64 = 0.0im
    cg_muiu::ComplexF64 = 0.0im
    cg_nujv::ComplexF64 = 0.0im
    cg_total::ComplexF64 = 0.0im

    # we need all quantum numbers for u,v
    # and the multiplets for mu,nu,i,j
    (q_u,m_mu,m_i) = u_mu_i
    (q_v,m_nu,m_j) = v_nu_j

    # Informative but expensive. Uncomment for testing.
    #if verbose 
    #    println( "COMPUTING (N-1) PSG" )
    #    println( "q_u' = $q_u ; q_v' = $q_v" )
    #    println( "m_i' = $m_i ; m_j' = $m_j" )
    #    println( "m_mu' = $m_mu ; m_nu' = $m_nu" ) 
    #    println( "q_a = $q_a" )
    #end

    # quantum numbers 
    # q_u, q_v 
    (N_u,I_u,S_u,i_u,s_u,r_u) = q_u
    (N_v,I_v,S_v,i_v,s_v,r_v) = q_v
    si_u = Int64((s_u+S_u)/2+1)
    si_v = Int64((s_v+S_v)/2+1)
    # m_mu, m_nu 
    (N_mu,I_mu,S_mu,r_mu) = m_mu
    (N_nu,I_nu,S_nu,r_nu) = m_nu
    # m_i, m_j
    (N_i,I_i,S_i,r_i) = (N_ij,I_ij,S_ij,r_ij) = m_i
    (N_j,I_j,S_j,r_j) = m_j

    # pseudo-cg matrix 
    G_mu = (N_mu,I_mu,S_mu)
    G_nu = (N_nu,I_nu,S_nu)
    G_a  = (N_a,I_a,S_a) = (q_a[1],q_a[2],q_a[3])
    (i_a,si_a,r_a) = (q_a[4],Int64((q_a[5]+S_a)/2+1),q_a[6]) 
    ((G_mu,G_a,G_nu) in keys(pcgmat)) || (return zero(ComplexF64))
    pcg_local = @view pcgmat[G_mu,G_a,G_nu][:,:,:,:,r_mu,r_nu,i_a,si_a,r_a]

    # clebsch-gordan matrices
    !( ((I_mu,I_ij,I_u) in keys(cg_o_fullmatint)) &&
       ((I_nu,I_ij,I_v) in keys(cg_o_fullmatint)) ) && (return zero(ComplexF64))
    @inbounds begin
        cgomat_muiu = @view cg_o_fullmatint[I_mu,I_ij,I_u][:,:,i_u]
        cgsmat_muiu = @view cg_s_fullmatint[S_mu,S_ij,S_u][:,:,si_u]
        cgomat_nujv = @view cg_o_fullmatint[I_nu,I_ij,I_v][:,:,i_v]
        cgsmat_nujv = @view cg_s_fullmatint[S_nu,S_ij,S_v][:,:,si_v]
    end

    # dimension of orbital irreps -> info about partners
    d_Iij = oindex2dimensions[I_i] # I_i = I_j
    d_Imu = oindex2dimensions[I_mu]
    d_Inu = oindex2dimensions[I_nu]

    # main iteration 
    @inbounds for (si_ij,s_ij) in enumerate((-S_ij):2:S_ij),
                  i_ij in 1:d_Iij,
                  (si_nu,s_nu) in enumerate((-S_nu):2:S_nu),
                  i_nu in 1:d_Inu,
                  (si_mu,s_mu) in enumerate((-S_mu):2:S_mu),
                  i_mu in 1:d_Imu

        i_i = i_j = i_ij 
        s_i = s_j = s_ij        

        if verbose 
            println( "g_mu' = ( $i_mu , $(s_mu/2.0) )" )
            println( "g_nu' = ( $i_nu , $(s_nu/2.0) )" )
        end

        # coefficient
        si_mu = Int64((s_mu+S_mu)/2+1)
        si_nu = Int64((s_nu+S_nu)/2+1)
        pcgcoeff = pcg_local[i_mu,si_mu,i_nu,si_nu]
        #q_mu = (N_mu,I_mu,S_mu,i_mu,s_mu,r_mu)
        #q_nu = (N_nu,I_nu,S_nu,i_nu,s_nu,r_nu)
        #pcgcoeff = get( pcg , (q_mu,q_a,q_nu) , zero(ComplexF64) )
        pcgcoeff==zero(pcgcoeff) && continue
        verbose && println( "pcgcoeff = $pcgcoeff")

        # cg coefficients
        cg_muiu = cgomat_muiu[i_mu,i_ij] * cgsmat_muiu[si_mu,si_ij]
        cg_nujv = cgomat_nujv[i_nu,i_ij] * cgsmat_nujv[si_nu,si_ij]
        cg_total = conj( cg_muiu ) * cg_nujv 

        if verbose 
            println( "cg_muiu = $cg_muiu" )
            println( "cg_nujv = $cg_nujv" )
            println( "cg_total = $cg_total" )
        end

        # putting all together
        coeff += cg_total * pcgcoeff
        verbose && println( "coeff = $coeff" )
    end
    verbose && println( "final coeff = $coeff" )
    verbose && println()
    return coeff 
end


# #####################################
# MATRIX ELEMENTS ( m_u | H_1 | m_v ) 
#
# - Compute the matrix elements for a 
#   step given the necessary info about 
#   the previous one 
# ..................................... 

function compute_step_matels( 
        multiplets_block::Set{NTuple{4,Int64}}, 
        multiplets_shell::Set{NTuple{4,Int64}},
        irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
        hop::Dict{ Tuple{Int64,Int64} , ComplexF64 },
        cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        pcg::Dict{ NTuple{3,NTuple{6,Int64}} , ComplexF64 },
        pcgmat::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,9} },
        qq_a::Vector{NTuple{6,Int64}}, 
        combinations_uprima::Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} },
        oindex2dimensions::Vector{ Int64 } ;
        verbose=false )
    # computes matrix elements ( m_u | H_1 | m_v ) for 
    # the next step.
    # input:
    # - multiplets_block : { q_i }
    # - multiplets_shell : { q_mu }
    # - irrEU : G => E, U (main result of previous step)
    # - hop : m_a => xi for n-th step
    # - cg_o_fullmatin : orbital CG coefficients for this problem (int format)
    # - pcg : pseudo-CG coefficients
    # - pcgmat : pseudo-CG coefficients in matrix form
    # - qq_a : set of all q_a
    # - combinations_uprima : m_u' => m_mu', m_i'
    # - oirreps2dimensions : orbital_irrep => dimension 
    # output:
    # - u_H1_v : ( m_u | H_1 | m_v )
    # - combinations_uprima_new : m_u => m_mu, m_i

    # type annotations
    m_up::Tuple{Int64,Int64,Int64,Int64}  = (0,0,0,0)
    m_vp::Tuple{Int64,Int64,Int64,Int64}  = (0,0,0,0)
    m_ip::Tuple{Int64,Int64,Int64,Int64}  = (0,0,0,0)
    m_jp::Tuple{Int64,Int64,Int64,Int64}  = (0,0,0,0)
    m_mup::Tuple{Int64,Int64,Int64,Int64} = (0,0,0,0)
    m_nup::Tuple{Int64,Int64,Int64,Int64} = (0,0,0,0)
    iu::Int64 = 0
    iv::Int64 = 0
    uu::ComplexF64      = 0
    hoparam::ComplexF64 = 0
    pseudo::ComplexF64  = 0
    iaj::ComplexF64     = 0
    

    # block-shell combination multiplets 
    multiplets_mui2u::Dict{NTuple{2,NTuple{4,Int64}},Vector{NTuple{4,Int64}}} = 
            get_combination_multiplets( multiplets_block , 
                                        multiplets_shell , 
                                        cg_o_fullmatint ;
                                        verbose=false )
    println( "COMBINATION MULTIPLETS:" )
    for (k,v) in multiplets_mui2u 
        @show k, v
    end

    # new block-shell combinations 
    GG_u = Set( m_u[1:3] for mm_u in values(multiplets_mui2u) 
                         for m_u in mm_u )
    combinations_uprima_new = Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} }(
                G => NTuple{3,NTuple{4,Int64}}[ (m_u,m_mu,m_i) 
                                                for ((m_mu,m_i),mm_u) in multiplets_mui2u 
                                                for m_u in mm_u
                                                if m_u[1:3]==G 
                                              ] 
                for G in GG_u
    )

    # Hamiltonian matrix
    newmults = Set{NTuple{4,Int64}}( m for mm in values(multiplets_mui2u) for m in mm )
    newirrepmults = get_irreps( newmults ; multiplicity=true )
    u_H1_v_blocks = Dict{ NTuple{3,Int64} , Array{ComplexF64,3} }(
            G => zeros(ComplexF64,R,R,2) for (G,R) in newirrepmults
    )
    
    # iteration over multiplets in i,mu and j,nu, 
    # and their combinations u,v
    for ((m_mu,m_i),mm_u) in multiplets_mui2u, 
        ((m_nu,m_j),mm_v) in multiplets_mui2u

        if verbose
            println( "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" )
            println( "m_mu = $m_mu, m_i = $m_i" )
            println( "m_nu = $m_nu, m_j = $m_j" )
            println( "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" )
        end

        # block and shell quantum numbers 
        (N_i, I_i, S_i, r_i ) = m_i
        (N_j, I_j, S_j, r_j ) = m_j
        G_i = (N_i, I_i, S_i)
        G_j = (N_j, I_j, S_j)
        (N_mu,I_mu,S_mu,r_mu) = m_mu
        (N_nu,I_nu,S_nu,r_nu) = m_nu
        
        (N_i+N_mu)==(N_j+N_nu) || continue # block-diagonality in N

        # U matrices (not transposed for faster loop)
        U_i = irrEU[G_i][2]
        U_j = irrEU[G_j][2]

        if (N_nu==(N_mu+1) && N_i==(N_j+1)) #hopping possible
            # (mu',i'=>u'), (nu',j'=>v') combinations 
            @inbounds begin
                combs_uprima_local = @view combinations_uprima[G_i][:]
                combs_vprima_local = @view combinations_uprima[G_j][:]
            end                          
        end

        # iterate over irreps of multiplets u,v in 
        # combinations (i,mu)=>u and (j,nu)=>v
        #for m_u::NTuple{4,Int64} in mm_u
        @showprogress for m_u::NTuple{4,Int64} in mm_u
            for m_v::NTuple{4,Int64} in mm_v

            r_u = m_u[4]
            r_v = m_v[4]

            if verbose 
                println( "..........................." )
                println( "( $m_u | H_1 | $m_v )" )
                println( "..........................." )
                println()
            end
            
            # irreps and block-diagonality
            G_u = m_u[1:3]
            G_v = m_v[1:3]
            G_u==G_v || continue # block-diagonality in irrep
            G = G_u 
            verbose && println( "G_u = G_v = $G" )

            # select H block view
            hblockview = @view u_H1_v_blocks[G][:,:,:]

            # chosen partner: orbital 1, max z-spin
            p = (N,I,S,i,s) = ( G... , 1 , G[3] )
            verbose && println( "p_u = p_v = $((N,I,S,i,s))" )

            # clebsch-gordan matrices
            @inbounds begin
                cgomat_muiu = @view cg_o_fullmatint[I_mu,I_i,I][:,:,i]
                cgomat_nujv = @view cg_o_fullmatint[I_nu,I_j,I][:,:,i]
            end

            # clebsch-gordan series
            verbose && println( "expanding CG series" )
            @inbounds for i_j::Int64  = 1:(oindex2dimensions[I_j]::Int64),
                          i_i::Int64  = 1:(oindex2dimensions[I_i]::Int64), 
                          i_nu::Int64 = 1:(oindex2dimensions[I_nu]::Int64),
                          i_mu::Int64 = 1:(oindex2dimensions[I_mu]::Int64), 
                          s_j::Int64  = (-S_j):2:S_j,
                          s_i::Int64  = (-S_i):2:S_i,
                          s_nu::Int64 = (-S_nu):2:S_nu, 
                          s_mu::Int64 = (-S_mu):2:S_mu
                
                # qnums
                q_i  = (N_i,I_i,S_i,i_i,s_i,r_i)
                q_j  = (N_j,I_j,S_j,i_j,s_j,r_j)
                q_mu = (N_mu,I_mu,S_mu,i_mu,s_mu,r_mu)
                q_nu = (N_nu,I_nu,S_nu,i_nu,s_nu,r_nu)
                G_i  = (N_i,I_i,S_i)
                G_j  = (N_j,I_j,S_j)

                # cg coefficients
                cg_muiu = cgomat_muiu[i_mu,i_i] * 
                          PartialWaveFunctions.CG_doublearg(S_mu,s_mu,S_i,s_i,S,s)
                cg_nujv = cgomat_nujv[i_nu,i_j] * 
                          PartialWaveFunctions.CG_doublearg(S_nu,s_nu,S_j,s_j,S,s)
                cg_tot  = conj( cg_muiu ) * cg_nujv
                isapprox( cg_tot , 0.0im ) && continue

                if verbose
                    println( "q_mu=$q_mu ; q_i=$q_i" )
                    println( "q_nu=$q_nu ; q_j=$q_j" )
                end

                # diagonal block energy
                E = (q_i==q_j && q_mu==q_nu) ? irrEU[G_i][1][r_i] : 
                                               zero(ComplexF64)
                verbose && println( "E = $E" )

                # hopping part iteration
                hopart = 0.0im
                
                # hopping number matching
                if (N_nu==(N_mu+1) && N_i==(N_j+1)) #hopping possible

                    verbose && println( "N_i=$N_i ; N_j=$N_j" )

                    verbose && println( "hopping iteration" )
                    for q_a in qq_a

                        verbose && println( "q_a = $q_a" )

                        # hopping parameter (orbital multiplet dependent)
                        hoparam = hop[(q_a[2],q_a[end])]
                        hoparam==0.0im && continue

                        # universal pseudo-CG coefficient (needs to be conjugated)
                        # MIGHT BE IMPROVED WITH PCGMAT 
                        nuamu = get( pcg , (q_nu,q_a,q_mu) , zero(ComplexF64) )
                        nuamu==zero(nuamu) && continue

                        # U-transformed coefficient
                        iaj = 0.0im
                        verbose && println( "u-transformed coefficients. iterating..." )
                        @inbounds for iu::Int64 in 1:length(combs_uprima_local)    

                            # i -- U --> u' 
                            (m_up,m_mup,m_ip) = combs_uprima_local[iu]::NTuple{3,NTuple{4,Int64}}
                            r_up = m_up[4]::Int64

                            q_up = ( q_i[1], q_i[2], q_i[3], q_i[4], q_i[5], r_up )
                            u_mu_i = ( q_up , m_mup , m_ip )

                            @inbounds for iv::Int64 in 1:length(combs_vprima_local)

                                # j -- U --> v'
                                (m_vp,m_nup,m_jp) = combs_vprima_local[iv]::NTuple{3,NTuple{4,Int64}}

                                # hopping only changes last shell states
                                m_ip==m_jp || continue

                                r_vp = m_vp[4]::Int64

                                verbose && println( "r_up=$r_up , r_vp=$r_vp" )

                                # multiplication of U matrix elements
                                uu = conj(U_i[r_up,r_i]) * U_j[r_vp,r_j]
                                uu==zero(uu) && continue
                                verbose && println( "uu = $uu" )

                                q_vp = ( q_j[1], q_j[2], q_j[3], q_j[4], q_j[5], r_vp )
                                v_nu_j = ( q_vp , m_nup , m_jp )

                                # pseudo-CG coefficient
                                pseudo = get_pseudoCG_up_vp( u_mu_i,
                                                             v_nu_j,
                                                             q_a,
                                                             oindex2dimensions ,
                                                             pcg , 
                                                             pcgmat ,
                                                             cg_o_fullmatint ;
                                                             verbose=false )::ComplexF64
                                verbose && println( "pseudo = $pseudo" )

                                iaj += hoparam * uu * pseudo
                                # Uncomment for testing.
                                verbose && println( "-> iaj += $(hoparam * uu * pseudo)" )
                            end
                        end
                        verbose && println( "q_a contribution (no sign) = $iaj" )
                        hopart += iaj * conj(nuamu)
                    end
                    hopart *= (-1)^N_mu
                end

                verbose && println( "hopart = $hopart" )

                # putting everything together
                hblockview[r_u,r_v,:] .+= cg_tot .* SVector{2,ComplexF64}( E , hopart ) 
            end
            verbose && println()
            verbose && println( "matel = $(hblockview[r_u,r_v,:])" )
            verbose && println()
        end
    end
    end

    # add h.c. 
    u_H1_v_new = Dict{ NTuple{3,Int64} , Matrix{ComplexF64} }(
        G => ComplexF64[ (sum(mat[r_u,r_v,:]) + conj(mat[r_v,r_u,2]))
                         for r_u in 1:size(mat,1), r_v in 1:size(mat,2) ] 
        for (G,mat) in u_H1_v_blocks
    )

    return ( u_H1_v_new, combinations_uprima_new )
end

function compute_step_matels_spinarray( 
        multiplets_block::Set{NTuple{4,Int64}}, 
        multiplets_shell::Set{NTuple{4,Int64}},
        irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
        hop::Dict{ Tuple{Int64,Int64} , ComplexF64 },
        cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        pcg::Dict{ NTuple{3,NTuple{6,Int64}} , ComplexF64 },
        pcgmat::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,9} },
        qq_a::Vector{NTuple{6,Int64}}, 
        combinations_uprima::Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} },
        oindex2dimensions::Vector{ Int64 } ;
        verbose=false )
    # computes matrix elements ( m_u | H_1 | m_v ) for 
    # the next step.
    # input:
    # - multiplets_block : { q_i }
    # - multiplets_shell : { q_mu }
    # - irrEU : G => E, U (main result of previous step)
    # - hop : m_a => xi for n-th step
    # - cg_o_fullmatin : orbital CG coefficients for this problem (int format)
    # - pcg : pseudo-CG coefficients
    # - pcgmat : pseudo-CG coefficients in matrix form
    # - qq_a : set of all q_a
    # - combinations_uprima : m_u' => m_mu', m_i'
    # - oirreps2dimensions : orbital_irrep => dimension 
    # output:
    # - u_H1_v : ( m_u | H_1 | m_v )
    # - combinations_uprima_new : m_u => m_mu, m_i

    # type annotations
    m_up::Tuple{Int64,Int64,Int64,Int64}  = (0,0,0,0)
    m_vp::Tuple{Int64,Int64,Int64,Int64}  = (0,0,0,0)
    m_ip::Tuple{Int64,Int64,Int64,Int64}  = (0,0,0,0)
    m_jp::Tuple{Int64,Int64,Int64,Int64}  = (0,0,0,0)
    m_mup::Tuple{Int64,Int64,Int64,Int64} = (0,0,0,0)
    m_nup::Tuple{Int64,Int64,Int64,Int64} = (0,0,0,0)
    iu::Int64 = 0
    iv::Int64 = 0
    uu::ComplexF64      = 0
    hoparam::ComplexF64 = 0
    pseudo::ComplexF64  = 0
    iaj::ComplexF64     = 0
    
    # block-shell combination multiplets 
    multiplets_mui2u::Dict{NTuple{2,NTuple{4,Int64}},Vector{NTuple{4,Int64}}} = 
            get_combination_multiplets( multiplets_block , 
                                        multiplets_shell , 
                                        cg_o_fullmatint ;
                                        verbose=false )

    # new block-shell combinations 
    GG_u = Set( m_u[1:3] for mm_u in values(multiplets_mui2u) 
                         for m_u in mm_u )
    combinations_uprima_new = Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} }(
                G => NTuple{3,NTuple{4,Int64}}[ (m_u,m_mu,m_i) 
                                                for ((m_mu,m_i),mm_u) in multiplets_mui2u 
                                                for m_u in mm_u
                                                if m_u[1:3]==G 
                                              ] 
                for G in GG_u
    )

    # Hamiltonian matrix
    newmults = Set{NTuple{4,Int64}}( m for mm in values(multiplets_mui2u) for m in mm )
    newirrepmults = get_irreps( newmults ; multiplicity=true )
    u_H1_v_blocks = Dict{ NTuple{3,Int64} , Array{ComplexF64,3} }(
            G => zeros(ComplexF64,R,R,2) for (G,R) in newirrepmults
    )
    
    # iteration over multiplets in i,mu and j,nu, 
    # and their combinations u,v
    for ((m_mu,m_i),mm_u) in multiplets_mui2u, 
        ((m_nu,m_j),mm_v) in multiplets_mui2u

        if verbose
            println( "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" )
            println( "m_mu = $m_mu, m_i = $m_i" )
            println( "m_nu = $m_nu, m_j = $m_j" )
            println( "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" )
        end

        # block and shell quantum numbers 
        (N_i, I_i, S_i, r_i ) = m_i
        (N_j, I_j, S_j, r_j ) = m_j
        G_i = (N_i, I_i, S_i)
        G_j = (N_j, I_j, S_j)
        (N_mu,I_mu,S_mu,r_mu) = m_mu
        (N_nu,I_nu,S_nu,r_nu) = m_nu
        
        (N_i+N_mu)==(N_j+N_nu) || continue # block-diagonality in N

        # U matrices (not transposed for faster loop)
        U_i = irrEU[G_i][2]
        U_j = irrEU[G_j][2]

        if (N_nu==(N_mu+1) && N_i==(N_j+1)) #hopping possible
            # (mu',i'=>u'), (nu',j'=>v') combinations 
            @inbounds begin
                combs_uprima_local = @view combinations_uprima[G_i][:]
                combs_vprima_local = @view combinations_uprima[G_j][:]
            end                          
        end

        # iterate over irreps of multiplets u,v in 
        # combinations (i,mu)=>u and (j,nu)=>v
        #@showprogress for m_u::NTuple{4,Int64} in mm_u
        for m_u::NTuple{4,Int64} in mm_u
            for m_v::NTuple{4,Int64} in mm_v

            r_u = m_u[4]
            r_v = m_v[4]

            if verbose 
                println( "..........................." )
                println( "( $m_u | H_1 | $m_v )" )
                println( "..........................." )
                println()
            end
            
            # irreps and block-diagonality
            G_u = m_u[1:3]
            G_v = m_v[1:3]
            G_u==G_v || continue # block-diagonality in irrep
            G = G_u 
            verbose && println( "G_u = G_v = $G" )

            # select H block view
            hblockview = @view u_H1_v_blocks[G][:,:,:]

            # chosen partner: orbital 1, max z-spin
            p = (N,I,S,i,s) = ( G... , 1 , G[3] )
            si = Int64( (s+S)/2.0 + 1 )
            verbose && println( "p_u = p_v = $((N,I,S,i,s))" )

            # clebsch-gordan matrices
            @inbounds begin
                cgomat_muiu = @view cg_o_fullmatint[I_mu,I_i,I][:,:,i]
                cgsmat_muiu = @view cg_s_fullmatint[S_mu,S_i,S][:,:,si]
                cgomat_nujv = @view cg_o_fullmatint[I_nu,I_j,I][:,:,i]
                cgsmat_nujv = @view cg_s_fullmatint[S_nu,S_j,S][:,:,si]
            end

            # clebsch-gordan series
            verbose && println( "expanding CG series" )
            @inbounds for i_j::Int64  = 1:(oindex2dimensions[I_j]::Int64),
                          i_i::Int64  = 1:(oindex2dimensions[I_i]::Int64), 
                          i_nu::Int64 = 1:(oindex2dimensions[I_nu]::Int64),
                          i_mu::Int64 = 1:(oindex2dimensions[I_mu]::Int64), 
                          (si_j::Int64,s_j::Int64) in enumerate((-S_j):2:S_j),
                          (si_i::Int64,s_i::Int64) in enumerate((-S_i):2:S_i),
                          (si_nu::Int64,s_nu::Int64) in enumerate((-S_nu):2:S_nu),
                          (si_mu::Int64,s_mu::Int64) in enumerate((-S_mu):2:S_mu)
                
                # qnums
                q_i  = (N_i,I_i,S_i,i_i,s_i,r_i)
                q_j  = (N_j,I_j,S_j,i_j,s_j,r_j)
                q_mu = (N_mu,I_mu,S_mu,i_mu,s_mu,r_mu)
                q_nu = (N_nu,I_nu,S_nu,i_nu,s_nu,r_nu)
                G_i  = (N_i,I_i,S_i)
                G_j  = (N_j,I_j,S_j)

                # cg coefficients
                cg_muiu = cgomat_muiu[i_mu,i_i] * cgsmat_muiu[si_mu,si_i]
                cg_nujv = cgomat_nujv[i_nu,i_j] * cgsmat_nujv[si_nu,si_j]
                cg_tot  = conj( cg_muiu ) * cg_nujv
                isapprox( cg_tot , 0.0im ) && continue

                if verbose
                    println( "q_mu=$q_mu ; q_i=$q_i" )
                    println( "q_nu=$q_nu ; q_j=$q_j" )
                end

                # diagonal block energy
                E = (q_i==q_j && q_mu==q_nu) ? irrEU[G_i][1][r_i] : 
                                               zero(ComplexF64)
                verbose && println( "E = $E" )

                # hopping part iteration
                hopart = 0.0im
                
                # hopping number matching
                if (N_nu==(N_mu+1) && N_i==(N_j+1)) #hopping possible

                    verbose && println( "N_i=$N_i ; N_j=$N_j" )

                    verbose && println( "hopping iteration" )
                    for q_a in qq_a

                        verbose && println( "q_a = $q_a" )

                        # hopping parameter (orbital multiplet dependent)
                        hoparam = hop[(q_a[2],q_a[end])]
                        hoparam==0.0im && continue

                        # universal pseudo-CG coefficient (needs to be conjugated)
                        # MIGHT BE IMPROVED WITH PCGMAT 
                        nuamu = get( pcg , (q_nu,q_a,q_mu) , zero(ComplexF64) )
                        nuamu==zero(nuamu) && continue

                        # U-transformed coefficient
                        iaj = 0.0im
                        verbose && println( "u-transformed coefficients. iterating..." )
                        @inbounds for iu::Int64 in 1:length(combs_uprima_local)    

                            # i -- U --> u' 
                            (m_up,m_mup,m_ip) = combs_uprima_local[iu]::NTuple{3,NTuple{4,Int64}}
                            r_up = m_up[4]::Int64

                            q_up = ( q_i[1], q_i[2], q_i[3], q_i[4], q_i[5], r_up )
                            u_mu_i = ( q_up , m_mup , m_ip )

                            @inbounds for iv::Int64 in 1:length(combs_vprima_local)

                                # j -- U --> v'
                                (m_vp,m_nup,m_jp) = combs_vprima_local[iv]::NTuple{3,NTuple{4,Int64}}

                                # hopping only changes last shell states
                                m_ip==m_jp || continue

                                r_vp = m_vp[4]::Int64

                                verbose && println( "r_up=$r_up , r_vp=$r_vp" )

                                # multiplication of U matrix elements
                                uu = conj(U_i[r_up,r_i]) * U_j[r_vp,r_j]
                                uu==zero(uu) && continue
                                verbose && println( "uu = $uu" )

                                q_vp = ( q_j[1], q_j[2], q_j[3], q_j[4], q_j[5], r_vp )
                                v_nu_j = ( q_vp , m_nup , m_jp )

                                # pseudo-CG coefficient
                                pseudo = get_pseudoCG_up_vp_spinarray( u_mu_i,
                                                                       v_nu_j,
                                                                       q_a,
                                                                       oindex2dimensions ,
                                                                       pcg , 
                                                                       pcgmat ,
                                                                       cg_o_fullmatint ,
                                                                       cg_s_fullmatint ;
                                                                       verbose=false )::ComplexF64
                                verbose && println( "pseudo = $pseudo" )

                                iaj += hoparam * uu * pseudo
                                # Uncomment for testing.
                                verbose && println( "-> iaj += $(hoparam * uu * pseudo)" )
                            end
                        end
                        verbose && println( "q_a contribution (no sign) = $iaj" )
                        hopart += iaj * conj(nuamu)
                    end
                    hopart *= (-1)^N_mu
                end

                verbose && println( "hopart = $hopart" )

                # putting everything together
                hblockview[r_u,r_v,:] .+= cg_tot .* SVector{2,ComplexF64}( E , hopart ) 
            end
            verbose && println()
            verbose && println( "matel = $(hblockview[r_u,r_v,:])" )
            verbose && println()
        end
    end
    end

    # add h.c. 
    u_H1_v_new = Dict{ NTuple{3,Int64} , Matrix{ComplexF64} }(
        G => ComplexF64[ (sum(mat[r_u,r_v,:]) + conj(mat[r_v,r_u,2]))
                         for r_u in 1:size(mat,1), r_v in 1:size(mat,2) ] 
        for (G,mat) in u_H1_v_blocks
    )

    return ( u_H1_v_new, combinations_uprima_new )
end


# ####################################
# GOMBINATIONS U'
#
#       m_i', G_u' => [ m_mu' , m_u' ]
# ....................................

function get_Gup2mip2mmupmup( combinations_mui2u )

    combinations_mui2u_flat = [ (m_mup,m_ip,m_up) 
                                for ((m_ip,m_mup),mm_up) in combinations_mui2u 
                                for m_up in mm_up ]
    Gup2mip = Dict{ Tuple{Int64,String,Float64} ,
                    Dict{Tuple{Int64,String,Float64,Int64},
                         Vector{NTuple{2,Tuple{Int64,String,Float64,Int64}}}}}()
    for (m_mup,m_ip,m_up) in combinations_mui2u_flat 
        mergewith!( (a,b)->mergewith!( vcat , a , b ) ,
                    Gup2mip ,
                    Dict( m_up[1:3]=>Dict(m_ip=>[(m_mup,m_up)]) ) )
    end

    #for (Gup,mipdict) in Gup2mip
    #    for (m_mup,m_ip,m_up) in combinations_mui2u_flat 
    #        m_up[1:3]==Gup || continue
    #        push!( Gup2mip[Gup] , mip=>(m_mup,m_up) )
    #    end
    #end

    return Gup2mip
end


# ###############################
# DIAGONALIZE ( m_u | H_1 | m_v ) 
#
# - Diagonalize symmetry-adapted 
#   Hamiltonian matrix by blocks
# ...............................

function diagonalize_step_blockform( 
        u_H1_v::Dict{ NTuple{3,Int64} , Matrix{ComplexF64} },
        oindex2dimensions::Vector{Int64} ;
        verbose=false )

    verbose && println( "DIAGONALIZATION" )
    # eigenvalues and transformation matrices
    irrEU = Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} }()

    # iterate through blocks
    for (G,H) in u_H1_v 
        
        verbose && println( "G = $G" )

        # dimensionality 
        Do = oindex2dimensions[G[2]]
        Ds = 2*G[3] + 1 
        D = Do*Ds

        # multiplicity
        R = size(H,1)

        # diagonalize matrix
        F = eigen( H )
        (e,u) = ( real(F.values) , F.vectors )
        push!( irrEU , G=>(e,u) )

        if verbose 
            println( "E = $e" )
            println()
        end
    end

    # subtract ground energy 
    minE = minimum(collect( minimum(v[1]) for v in values(irrEU) ))
    irrEU = Dict( G=>(E.-minE,U) for (G,(E,U)) in irrEU )

    return irrEU
end


function diagonalize_step( 
            u_H1_v::Dict , 
            oirreps2dimensions::Dict{String,Int64} 
            ; verbose=false )

    verbose && println( "DIAGONALIZATION" )
    # blocks = irreps
    GG = Set( k[1][1:3] for k in keys(u_H1_v) )

    # multiplets 
    mm = Set( (k[1][1:3]...,k[1][end]) for k in keys(u_H1_v) )

    # eigenvalues and transformation matrices
    irrEU = Dict{ Tuple{Int64,String,Float64} , Tuple{Vector{Float64},Matrix{ComplexF64}} }()

    # iterate through blocks
    for G in GG 
        
        verbose && println("G = $G" )

        # dimensionality 
        Do = oirreps2dimensions[G[2]]
        Ds = 2*G[3] + 1 
        D = Do*Ds

        # multiplicity
        R = get_multiplicity( mm , G )

        # iterate through block to construct matrix
        mat = zeros( ComplexF64 , R , R )
        for r_u=1:R, r_v=1:R
            mat[r_u,r_v] = u_H1_v[((G...,r_u),(G...,r_v))]
        end

        if verbose 
            println( "H_sub" )
            for i=1:size(mat,1)
                println( mat[i,:] )
            end
        end

        # diagonalize matrix
        F = eigen( mat )
        (e,u) = ( real(F.values) , F.vectors )
        push!( irrEU , G=>(e,u) )

        if verbose 
            println( "E = $e" )
            println()
        end
    end

    # subtract ground energy 
    minE = minimum(collect( minimum(v[1]) for v in values(irrEU) ))
    irrEU = Dict( G=>(E.-minE,U) for (G,(E,U)) in irrEU )

    return irrEU
end

function diagonalize_step_blockform( 
            u_H1_v::Dict , 
            oirreps2dimensions::Dict{String,Int64} 
            ; verbose=false )

    verbose && println( "DIAGONALIZATION" )
    # eigenvalues and transformation matrices
    irrEU = Dict{ Tuple{Int64,String,Float64} , Tuple{Vector{Float64},Matrix{ComplexF64}} }()

    # iterate through blocks
    for (G,H) in u_H1_v 
        
        verbose && println( "G = $G" )

        # dimensionality 
        Do = oirreps2dimensions[G[2]]
        Ds = 2*G[3] + 1 
        D = Do*Ds

        # multiplicity
        R = size(H,1)

        # diagonalize matrix
        @show H
        F = eigen( H )
        (e,u) = ( real(F.values) , F.vectors )
        push!( irrEU , G=>(e,u) )

        if verbose 
            println( "E = $e" )
            println()
        end
    end

    # subtract ground energy 
    minE = minimum(collect( minimum(v[1]) for v in values(irrEU) ))
    irrEU = Dict( G=>(E.-minE,U) for (G,(E,U)) in irrEU )

    return irrEU
end

# ###############################
# CUTOFF 
# 
# - Introduce energy/state cutoff 
#   in irrEU 
# ...............................

function cut_off!( 
            irrEU::Dict{ Tuple{Int64,Int64,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} };
            type::String="multiplet" , 
            cutoff::T=200 , 
            safeguard::Bool=true , 
            minmult::Int64=0 , 
            verbose::Bool=true ) where {T<:Number}
    # input 
    # - irrEU : G => (E,U)
    # - type : "multiplet" --> keep a limited number of lowest-energy multiplets
    #          "energy"    --> upper bound to energy
    # - cutoff : upper bound in multiplets/energy
    # output: 
    # - [irrEU : altered, not returned]
    # - { kept } : kept block multiplets
    # - discarded : amount of multiplets discarded

    # order multiplets in terms of their energy 
    mm::Vector{Tuple{NTuple{4,Int64},Float64}} = collect( ((G...,r),E[r]) for (G,(E,U)) in irrEU for r=1:length(E) )
    discarded::Vector{NTuple{4,Int64}} = []
    kept::Vector{NTuple{4,Int64}} = []
    sg::Vector{NTuple{4,Int64}} = []
    sort!( mm , by=x->x[2] )

    # multiplet cutoff
    if type=="multiplet"
        if cutoff>=length(mm)
            kept = map( x->x[1] , mm )
            discarded = []
        else
            kept = map( x->x[1] , mm[1:cutoff]::Vector{Tuple{NTuple{4,Int64},Float64}} )::Vector{NTuple{4,Int64}}
            if safeguard
                sg = map( x->x[1] ,
                        [ m for m in mm[(cutoff+1):end]::Vector{Tuple{NTuple{4,Int64},Float64}} 
                          if isapprox(m[2],mm[cutoff][2];atol=1e-1)] )::Vector{NTuple{4,Int64}}
                append!( kept , sg )
            end
            discarded = map( x->x[1] , 
                             mm[(cutoff+length(safeguard)+1):end]::Vector{Tuple{NTuple{4,Int64},Float64}} )::Vector{NTuple{4,Int64}}
        end
    # energy cutoff
    elseif type=="energy" 
        kept      = collect( pair[1] for pair in mm if pair[2]<=cutoff  )
        cutoff    = length(kept)
        (cutoff<minmult && minmult<length(mm)) && (cutoff=minmult)
        kept      = collect( pair[1] for pair in mm[1:cutoff] )
        if cutoff<length(mm) 
            if safeguard
                safeguard = map( x->x[1] ,
                                [ m for m in mm[(cutoff+1):end] 
                                  if isapprox(m[2],mm[cutoff][2];atol=1e-1)] )
                append!( kept , safeguard )
            end
            cutoff = length(kept)
            discarded = map( x->x[1] , mm[(cutoff+1):end] )
        else 
            discarded = []
        end
    else 
        println( "error: type must be 'multiplet' or 'energy'" )
        return "error: type must be 'multiplet' or 'energy'" 
    end

    # apply cutoff to irrEU
    for (G,(E,U)) in irrEU 
        N = length(collect( k for k in kept if k[1:3]==G ))
        if N==0 
            pop!( irrEU , G )
            continue
        end
        # U does not need cutoff (multiplets run over)
        irrEU[G] = ( E[1:N] , U )
    end

    return ( Set(kept) , discarded )::Tuple{Set{NTuple{4,Int64}},Vector{NTuple{4,Int64}}}
end



# #############
# MISCELLANEOUS
# .............

function print_irrEU_lowmults( irrEU , nmults::Int64 ; grandparents=Set() ) 
    mults = collect( (G,e) for (G,(E,U)) in irrEU for e in E )
    sort!( mults , by=x->x[2] ) 
    println( "Low-energy multiplets:" )
    for (i,(G,e)) in enumerate(mults)
        if grandparents==Set()
            println( "$G => $e" )
        else 
            println( "$G => $e || $(collect(grandparents[G]))" )
        end
        i==nmults && break
    end
end

function print_irrEU( irrEU::Dict{ Tuple{Int64,Int64,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} }, 
                      oirreps2dimensions , oirreps )
    for (G,(E,U)) in irrEU
        println( "G = $G" )
        println( "size = $(convert( Int64 , oirreps2dimensions[oirreps[G[2]]] * (G[3]+1) ))")
        println( "E = $E" )
        println( "U = $U" )
        println()
    end
end

function print_irrEU( irrEU::Dict{ Tuple{Int64,String,Float64} , Tuple{Vector{Float64},Matrix{ComplexF64}} }, 
                      oirreps2dimensions , oirreps )
    for (G,(E,U)) in irrEU
        println( "G = $G" )
        println( "size = $(convert(Int64,oirreps2dimensions[G[2]] * (2*G[3]+1) ))")
        println( "E = $E" )
        println( "U = $U" )
        println()
    end
end



# #############################
# CONVERSION TO INT KEYS FORMAT
# 
# - Convert key-type elements,
#   usually containing Float 
#   and String elements, to 
#   Int-only elements
# .............................

function convert_to_int( multiplet::Tuple{Int64,String,Float64,Int64} , oirreps2indices )
    ( N , I , S , r ) = multiplet
    return ( N , oirreps2indices[I] , Int64(2*S) , r )
end

function convert_to_int( irrep::Tuple{Int64,String,Float64} , oirreps2indices )
    ( N , I , S ) = irrep 
    return ( N , oirreps2indices[I] , Int64(2*S) )
end

function convert_to_int( irrep::Tuple{Int64,String,Int64} , oirreps2indices )
    ( N , I , S ) = irrep 
    return ( N , oirreps2indices[I] , Int64(2*S) )
end

function convert_to_int( q::Tuple{Int64,String,Number,Int64,Number,Int64} , oirreps2indices )
    ( N , I , S , i , s , r ) = q
    Ii = oirreps2indices[I]
    Si = Int64(2*S)
    si = Int64(2*s)
    return (N,Ii,Si,i,si,r)
end

function convert_to_int( multiplet::Tuple{Int64,AbstractString,Float64,Int64} , oirreps2indices )
    ( N , I , S , r ) = multiplet
    return ( N , oirreps2indices[I] , Int64(2*S) , r )
end

function convert_to_int( irrep::Tuple{Int64,AbstractString,Float64} , oirreps2indices )
    ( N , I , S ) = irrep 
    return ( N , oirreps2indices[I] , Int64(2*S) )
end

function convert_to_int( irrep::Tuple{Int64,AbstractString,Int64} , oirreps2indices )
    ( N , I , S ) = irrep 
    return ( N , oirreps2indices[I] , Int64(2*S) )
end

function convert_to_int( q::Tuple{Int64,AbstractString,Number,Int64,Number,Int64} , oirreps2indices )
    ( N , I , S , i , s , r ) = q
    Ii = oirreps2indices[I]
    Si = Int64(2*S)
    si = Int64(2*s)
    return (N,Ii,Si,i,si,r)
end


# #############
# NRG PROCEDURE
# .............

function NRG_onestep_test( iterations::Int64, 
                      cutoff_type::String, 
                      cutoff_magnitude::Number,
                      L::Float64,
                      hopchannels,
                      irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
                      cg_o_fullmatint::Dict{Tuple{Int64, Int64, Int64}, Array{ComplexF64, 3}},
                      pcg::Dict{Tuple{NTuple{6, Int64}, NTuple{6, Int64}, NTuple{6, Int64}}, ComplexF64}, 
                      pcgmat::Dict{Tuple{Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64, Int64}}, Array{ComplexF64, 9}}, 
                      qq_a::Vector{NTuple{6,Int64}}, 
                      combinations_uprima::Dict{NTuple{3,Int64}, Vector{NTuple{3,NTuple{4,Int64}}}},
                      betabar::Float64 ,
                      oindex2dimensions::Vector{Int64} ,
                      #ss_i , nn_i , 
                      mm_i ;
                      spinarray=false ,
                      cg_s_fullmatint::Dict{Tuple{Int64, Int64, Int64}, Array{ComplexF64, 3}}=Dict() ,
                      distributed=false ,
                      method="distfor" ,
                      minmult=0 , 
                      z::Float64=0.0 ,
                      discretization="standard" )

    println( "=============" )
    println( "NRG PROCEDURE" )
    println( "=============" )
    println()

    temperatures   = []
    magnetizations = []
    energies       = []
    numbers        = []
    partitions     = []
    entropies      = []
    impspins       = []
    impnums        = []
    impmults       = []

    performance = []
    maxes = []
    eigenvals5::Matrix{Float64} = zeros(Float64,iterations-1,5)

    ξ = compute_xi_vector( iterations , z , L ; discretization=discretization )

    for n=2:iterations

        println( "ITERATION n=$n" )

        # cut off and block block multiplets
        print( "Applying cutoff to obtain block multiplets... " )
        @time (multiplets_block, discarded) = cut_off!( irrEU ; 
                                                        type=cutoff_type , 
                                                        cutoff=cutoff_magnitude , 
                                                        safeguard=true ,
                                                        minmult=minmult ,
                                                        verbose=true )
        equivstates = Int64( sum( (m[3]+1) for m in multiplets_block ) )
        println( "$(length(multiplets_block)) multiplets kept ($equivstates states), $(length(discarded)) multiplets discarded" )
        proportion = length(multiplets_block)/Float64( length(multiplets_block) + length(discarded) ) * 100
        maxe = maximum(collect( e for (G,(E,U)) in irrEU for e in E ))
        push!( maxes , maxe )
        maxs = maximum(collect( G[3] for (G,(E,U)) in irrEU ))
        println( "proportion: $(proportion)%. max energy: $maxe. max spin: $maxs" )

        # renormalize by √Λ
        println( "Renormalizing eigenvalues...")
        for (G,(E,U)) in irrEU 
            irrEU[G] = ( E.*sqrt(L) , U )
        end

        # hopping parameter
        hop = Dict( hopchannel=>ComplexF64(ξ[n-1]) # at n=2, we want ξ[1]=ξ_0
                    for hopchannel in hopchannels )
        println( "shell hopping = $(ξ[n-1])" )

        # construct ( m_u | H_1 | m_v )
        println( "Diagonalizing Hamiltonian..." )
        #if distributed
            # spinarray must be true
        #    ppp = @timed (irrEU_dist,combinations_uprima_dist) = 
        #                matdiag_distributed( 
        #                            multiplets_block , 
        #                            multiplets_shell ,
        #                            irrEU , 
        #                            hop , 
        #                            cg_o_fullmatint , 
        #                            cg_s_fullmatint ,
        #                            pcg , 
        #                            pcgmat , 
        #                            qq_a , 
        #                            combinations_uprima , 
        #                            oindex2dimensions ;
        #                            method=method ,
        #                            verbose=false )
        #else 
        ppp = @timed (irrEU,combinations_uprima) = 
                                matdiag_serial( 
                                    multiplets_block , 
                                    multiplets_shell ,
                                    irrEU , 
                                    hop , 
                                    cg_o_fullmatint , 
                                    cg_s_fullmatint ,
                                    pcg , 
                                    pcgmat , 
                                    qq_a , 
                                    combinations_uprima , 
                                    oindex2dimensions ;
                                    verbose=false )
        #end
        @show ppp.time, ppp.bytes*10^-6, ppp.gctime
        push!( performance , ppp )

        # impurity properties 
        #ss_i,nn_i = imp_qnums( irrEU ,
        #                       oindex2dimensions ,
        #                       combinations_uprima ,
        #                       ss_i ,
        #                       nn_i )
        #s_imp = impspin( irrEU ,
        #                 betabar , 
        #                 oindex2dimensions ,
        #                 ss_i )
        #n_imp = impnum( irrEU ,
        #                betabar , 
        #                oindex2dimensions ,
        #                nn_i )
        #@show s_imp 
        #@show n_imp
        mm_i = imp_mults( irrEU ,
                          oindex2dimensions ,
                          combinations_uprima ,
                          mm_i )
        m_imp = mult_thermo( irrEU ,
                             betabar ,
                             oindex2dimensions ,
                             mm_i )
        @show m_imp
        #push!( impspins , s_imp ) 
        #push!( impnums  , n_imp )
        push!( impmults , m_imp )

        # thermodynamics 
        println( "THERMODYNAMICS" )
        t = temperature( n , L , betabar ; z=z )
        ρ = partition( irrEU , betabar , oindex2dimensions )
        entr= entropy( irrEU , betabar , oindex2dimensions )
        mag = magsusc( irrEU , betabar , oindex2dimensions )
        N = number(    irrEU , betabar , oindex2dimensions )
        en = energy(   irrEU , betabar , oindex2dimensions )
        @show t 
        @show ρ 
        @show entr
        @show mag
        @show N
        @show en
        println()
        push!( magnetizations , mag )
        push!( temperatures , t )
        push!( energies , en )
        push!( partitions , ρ )
        push!( numbers , N )
        push!( entropies , entr )

        eigs = sort([ e for (E,U) in values(irrEU) for e in E ])[1:5]
        eigenvals5[(n-1),:] = [(e-eigs[1]) for e in eigs[1:5]]

    end

    writedlm( "eigenvals_z$z.dat" , eigenvals5 )

    maxe_avg = sum(maxes)/length(maxes)
    @show maxe_avg

    return ( t=temperatures , 
             m=magnetizations , 
             e=energies , 
             p=partitions ,
             n=numbers , 
             entr=entropies,
             perf=performance ,
             #impspins=impspins ,
             #impnums=impnums ,
             impmults=impmults )
end
function NRG_onestep( iterations::Int64, 
                      cutoff_type::String, 
                      cutoff_magnitude::Number,
                      L::Float64,
                      hopchannels,
                      irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
                      multiplets_shell::Set{NTuple{4,Int64}}, # NEW ADDITION!!!!
                      cg_o_fullmatint::Dict{Tuple{Int64, Int64, Int64}, Array{ComplexF64, 3}},
                      pcg::Dict{Tuple{NTuple{6, Int64}, NTuple{6, Int64}, NTuple{6, Int64}}, ComplexF64}, 
                      pcgmat::Dict{Tuple{Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64, Int64}}, Array{ComplexF64, 9}}, 
                      qq_a::Vector{NTuple{6,Int64}}, 
                      combinations_uprima::Dict{NTuple{3,Int64}, Vector{NTuple{3,NTuple{4,Int64}}}},
                      betabar::Float64 ,
                      oindex2dimensions::Vector{Int64} ,
                      mm_i ;
                      spinarray=false ,
                      cg_s_fullmatint::Dict{Tuple{Int64, Int64, Int64}, Array{ComplexF64, 3}}=Dict() ,
                      distributed=false ,
                      method="distfor" ,
                      minmult=0 , 
                      z::Float64=0.0 ,
                      discretization="standard" )

    println( "=============" )
    println( "NRG PROCEDURE" )
    println( "=============" )
    println()

    temperatures   = []
    magnetizations = []
    energies       = []
    numbers        = []
    partitions     = []
    entropies      = []
    impspins       = []
    impnums        = []
    impmults       = []

    performance = []
    maxes = []
    eigenvals5::Matrix{Float64} = zeros(Float64,iterations-1,5)

    ξ = compute_xi_vector( iterations , z , L ; discretization=discretization )

    for n=2:iterations

        println( "ITERATION n=$n" )

        # cut off and block block multiplets
        print( "Applying cutoff to obtain block multiplets... " )
        @time (multiplets_block, discarded) = cut_off!( irrEU ; 
                                                        type=cutoff_type , 
                                                        cutoff=cutoff_magnitude , 
                                                        safeguard=true ,
                                                        minmult=minmult ,
                                                        verbose=true )
        equivstates = Int64( sum( (m[3]+1) for m in multiplets_block ) )
        println( "$(length(multiplets_block)) multiplets kept ($equivstates states), $(length(discarded)) multiplets discarded" )
        proportion = length(multiplets_block)/Float64( length(multiplets_block) + length(discarded) ) * 100
        maxe = maximum(collect( e for (G,(E,U)) in irrEU for e in E ))
        push!( maxes , maxe )
        maxs = maximum(collect( G[3] for (G,(E,U)) in irrEU ))
        println( "proportion: $(proportion)%. max energy: $maxe. max spin: $maxs" )

        # renormalize by √Λ
        println( "Renormalizing eigenvalues...")
        for (G,(E,U)) in irrEU 
            irrEU[G] = ( E.*sqrt(L) , U )
        end

        # hopping parameter
        hop = Dict( hopchannel=>ComplexF64(ξ[n-1]) # at n=2, we want ξ[1]=ξ_0
                    for hopchannel in hopchannels )
        @show hop
        println( "shell hopping = $(ξ[n-1])" )

        # construct ( m_u | H_1 | m_v )
        println( "Diagonalizing Hamiltonian..." )
        if distributed
            # spinarray must be true
            #ppp = @timed (irrEU,combinations_uprima) = 
            #            matdiag_distributed( 
            #                        multiplets_block , 
            #                        multiplets_shell ,
            #                        irrEU , 
            #                        hop , 
            #                        cg_o_fullmatint , 
            #                        cg_s_fullmatint ,
            #                        pcg , 
            #                        pcgmat , 
            #                        qq_a , 
            #                        combinations_uprima , 
            #                        oindex2dimensions ;
            #                        method=method ,
            #                        verbose=false )
            ppp = @timed (irrEU,combinations_uprima) = 
                        matdiag_distributed_H0( 
                                    multiplets_block , 
                                    multiplets_shell ,
                                    irrEU , 
                                    hop , 
                                    cg_o_fullmatint , 
                                    cg_s_fullmatint ,
                                    pcg , 
                                    pcg , 
                                    qq_a , 
                                    combinations_uprima , 
                                    oindex2dimensions ;
                                    method=method ,
                                    verbose=false )
        else 
            #if n==10
            #    @profile begin 
            global ppp = @timed (irrEU,combinations_uprima) = 
                                    matdiag_serial_H0( 
                                        multiplets_block , 
                                        multiplets_shell ,
                                        irrEU , 
                                        hop , 
                                        cg_o_fullmatint , 
                                        cg_s_fullmatint ,
                                        pcg , 
                                        pcg , 
                                        qq_a , 
                                        combinations_uprima , 
                                        oindex2dimensions ;
                                        verbose=false )
#                    println( "irrEU at N=$n" )
#                    for (G,(E,U)) in irrEU 
#                        @show G 
#                        @show E 
#                    end
            #    end
            #else
            #    global ppp = @timed (irrEU,combinations_uprima) = 
            #                    matdiag_serial( 
            #                            multiplets_block , 
            #                            multiplets_shell ,
            #                            irrEU , 
            #                            hop , 
            #                            cg_o_fullmatint , 
            #                            cg_s_fullmatint ,
            #                            pcg , 
            #                            pcgmat , 
            #                            qq_a , 
            #                            combinations_uprima , 
            #                            oindex2dimensions ;
            #                            verbose=false )
            #end 
#        else
#            if !spinarray
#                ppp = @timed (irrEU,combinations_uprima) = 
#                            matdiag( 
#                                        multiplets_block , 
#                                        multiplets_shell ,
#                                        irrEU , 
#                                        hop , 
#                                        cg_o_fullmatint , 
#                                        pcg , 
#                                        pcgmat , 
#                                        qq_a , 
#                                        combinations_uprima , 
#                                        oindex2dimensions ;
#                                        verbose=false )
#            else
#                ppp = @timed (irrEU,combinations_uprima) = 
#                            matdiag_spinarray( 
#                                        multiplets_block , 
#                                        multiplets_shell ,
#                                        irrEU , 
#                                        hop , 
#                                        cg_o_fullmatint , 
#                                        cg_s_fullmatint ,
#                                        pcg , 
#                                        pcgmat , 
#                                        qq_a , 
#                                        combinations_uprima , 
#                                        oindex2dimensions ;
#                                        verbose=false )
#            end
        end
        @show ppp.time, ppp.bytes*10^-6, ppp.gctime
        push!( performance , ppp )

        mm_i = imp_mults( irrEU ,
                          oindex2dimensions ,
                          combinations_uprima ,
                          mm_i )
        m_imp = mult_thermo( irrEU ,
                             betabar ,
                             oindex2dimensions ,
                             mm_i )
        @show m_imp
        #push!( impspins , s_imp ) 
        #push!( impnums  , n_imp )
        push!( impmults , m_imp )

        # thermodynamics 
        println( "THERMODYNAMICS" )
        t = temperature( n , L , betabar ; z=z )
        ρ = partition( irrEU , betabar , oindex2dimensions )
        entr= entropy( irrEU , betabar , oindex2dimensions )
        mag = magsusc( irrEU , betabar , oindex2dimensions )
        N = number(    irrEU , betabar , oindex2dimensions )
        en = energy(   irrEU , betabar , oindex2dimensions )
        @show t 
        @show ρ 
        @show entr
        @show mag
        @show N
        @show en
        println()
        push!( magnetizations , mag )
        push!( temperatures , t )
        push!( energies , en )
        push!( partitions , ρ )
        push!( numbers , N )
        push!( entropies , entr )

        eigs = sort([ e for (E,U) in values(irrEU) for e in E ])[1:5]
        eigenvals5[(n-1),:] = [(e-eigs[1]) for e in eigs[1:5]]

    end

    writedlm( "eigenvals_z$z.dat" , eigenvals5 )

    maxe_avg = sum(maxes)/length(maxes)
    @show maxe_avg

    return ( t=temperatures , 
             m=magnetizations , 
             e=energies , 
             p=partitions ,
             n=numbers , 
             entr=entropies,
             perf=performance ,
             #impspins=impspins ,
             #impnums=impnums ,
             impmults=impmults )
end


function NRG( iterations::Int64, 
              cutoff_type::String, 
              cutoff_magnitude::Number,
              L::Float64,
              hopchannels,
              irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
              cg_o_fullmatint::Dict{Tuple{Int64, Int64, Int64}, Array{ComplexF64, 3}},
              pcg::Dict{Tuple{NTuple{6, Int64}, NTuple{6, Int64}, NTuple{6, Int64}}, ComplexF64}, 
              pcgmat::Dict{Tuple{Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64, Int64}}, Array{ComplexF64, 9}}, 
              qq_a::Vector{NTuple{6,Int64}}, 
              combinations_uprima::Dict{NTuple{3,Int64}, Vector{NTuple{3,NTuple{4,Int64}}}},
              oindex2dimensions::Vector{Int64} ;
              spinarray=false ,
              cg_s_fullmatint::Dict{Tuple{Int64, Int64, Int64}, Array{ComplexF64, 3}}=Dict() ,
              distributed=false ,
              method="distfor" )

    println( "NRG PROCEDURE" )

    temperatures   = []
    magnetizations = []
    energies       = []
    numbers        = []
    partitions     = []

    performance_mat = []
    performance_dia = []

    for n=2:iterations

        println( "ITERATION n=$n" )

        # cut off and block block multiplets
        print( "Applying cutoff to obtain block multiplets... " )
        @time (multiplets_block, discarded) = cut_off!( irrEU ; 
                                                        type=cutoff_type , 
                                                        cutoff=cutoff_magnitude , 
                                                        safeguard=false ,
                                                        verbose=true )
        equivstates = Int64( sum( (m[3]+1) for m in multiplets_block ) )
        println( "$(length(multiplets_block)) multiplets kept ($equivstates states), $(length(discarded)) multiplets discarded" )
        proportion = length(multiplets_block)/Float64( length(multiplets_block) + length(discarded) ) * 100
        maxe = maximum(collect( e for (G,(E,U)) in irrEU for e in E ))
        maxs = maximum(collect( G[3] for (G,(E,U)) in irrEU ))
        println( "proportion: $(proportion)%. max energy: $maxe. max spin: $maxs" )

        # renormalize by L^{1/2}
        println( "Renormalizing eigenvalues...")
        for (G,(E,U)) in irrEU 
            irrEU[G] = ( E.*sqrt(L) , U )
        end

        # hopping parameter
        hop = Dict( hopchannel=>ComplexF64(xi(n,L)) 
                    for hopchannel in hopchannels )
        println( "shell hopping = $(xi(n,L))" )

        # construct ( m_u | H_1 | m_v )
        println( "Constructing Hamiltonian matrix..." )
        if distributed
            # spinarray must be true
            ppp = @timed (u_H1_v,combinations_uprima) = 
                        compute_step_matels_distributed( 
                                    multiplets_block , 
                                    multiplets_shell ,
                                    irrEU , 
                                    hop , 
                                    cg_o_fullmatint , 
                                    cg_s_fullmatint ,
                                    pcg , 
                                    pcgmat , 
                                    qq_a , 
                                    combinations_uprima , 
                                    oindex2dimensions ;
                                    method=method ,
                                    verbose=false )
        else
            if !spinarray
                ppp = @timed (u_H1_v,combinations_uprima) = 
                            compute_step_matels( 
                                        multiplets_block , 
                                        multiplets_shell ,
                                        irrEU , 
                                        hop , 
                                        cg_o_fullmatint , 
                                        pcg , 
                                        pcgmat , 
                                        qq_a , 
                                        combinations_uprima , 
                                        oindex2dimensions ;
                                        verbose=false )
            else
                ppp = @timed (u_H1_v,combinations_uprima) = 
                            compute_step_matels_spinarray( 
                                        multiplets_block , 
                                        multiplets_shell ,
                                        irrEU , 
                                        hop , 
                                        cg_o_fullmatint , 
                                        cg_s_fullmatint ,
                                        pcg , 
                                        pcgmat , 
                                        qq_a , 
                                        combinations_uprima , 
                                        oindex2dimensions ;
                                        verbose=false )
            end

        end
        @show ppp.time, ppp.bytes, ppp.gctime
        push!( performance_mat , ppp )

        # diagonalize 
        println( "Diagonalizing..." )
        ppp = @timed irrEU = diagonalize_step_blockform( u_H1_v , 
                                                        oindex2dimensions )
        @show ppp.time, ppp.bytes, ppp.gctime
        push!( performance_dia , ppp )

        print_irrEU_lowmults( irrEU , 10 )

        # thermodynamics 
        println( "THERMODYNAMICS" )
        t = temperature( n , L , betabar )
        Z = partition(irrEU,betabar,oindex2dimensions)
        mag = magsusc( irrEU , betabar , oindex2dimensions )
        N = number( irrEU , betabar , oindex2dimensions  )
        en = energy( irrEU , betabar , oindex2dimensions )
        @show t 
        @show Z
        @show mag
        @show N
        @show en
        println()
        push!( magnetizations , mag )
        push!( temperatures , t )
        push!( energies , en )
        push!( partitions , Z )
        push!( numbers , N )
    end

    return ( t=temperatures , 
             m=magnetizations , 
             e=energies , 
             p=partitions ,
             n=numbers , 
             p_mat=performance_mat ,
             p_dia=performance_dia )
end


# #######################
# PERFORMANCE INFORMATION
# .......................

function print_performance_onestep( nrg )
    performance = nrg.perf

    times  = [ p.time for p in performance ]
    t_min  = minimum( times )
    t_max  = maximum( times )
    t_mean = sum(times)/length(times)
    println( "TIME (s)" )
    println( "minimum time: $t_min" )
    println( "maximum time: $t_max" )
    println( "mean time:    $t_mean" )

    println()

    bytes  = [ p.bytes for p in performance]
    b_min  = minimum( bytes )/10^9
    b_max  = maximum( bytes )/10^9
    b_mean = sum(bytes)/length(bytes)/10^9
    println( "minimum gigabytes: $b_min" )
    println( "maximum gigabytes: $b_max" )
    println( "mean gigabytes:    $b_mean" )

    println()

    gc = [ p.gctime for p in performance]
    gc_min = minimum( gc )
    gc_max = maximum( gc )
    gc_mean = sum(gc)/length(gc)
    println( "minimum gc time: $gc_min" )
    println( "maximum gc time: $gc_max" )
    println( "mean gc time:    $gc_mean" )
end

function print_performance_mat( nrg )
    performance_mat = nrg.p_mat

    println( "MATRIX CONSTRUCTION PERFORMANCE" )
    println()

    times_mat = [ p.time for p in performance_mat ]
    t_min_mat = minimum( times_mat )
    t_max_mat = maximum( times_mat )
    t_mean_mat = sum(times_mat)/length(times_mat)
    println( "TIME (s)" )
    println( "minimum time: $t_min_mat" )
    println( "maximum time: $t_max_mat" )
    println( "mean time:    $t_mean_mat" )

    println()

    bytes_mat = [ p.bytes for p in performance_mat ]
    b_min_mat = minimum( bytes_mat )/10^9
    b_max_mat = maximum( bytes_mat )/10^9
    b_mean_mat = sum(bytes_mat)/length(bytes_mat)/10^9
    println( "minimum gigabytes: $b_min_mat" )
    println( "maximum gigabytes: $b_max_mat" )
    println( "mean gigabytes:    $b_mean_mat" )

    println()

    gc_mat = [ p.gctime for p in performance_mat ]
    gc_min_mat = minimum( gc_mat )
    gc_max_mat = maximum( gc_mat )
    gc_mean_mat = sum(gc_mat)/length(gc_mat)
    println( "minimum gc time: $gc_min_mat" )
    println( "maximum gc time: $gc_max_mat" )
    println( "mean gc time:    $gc_mean_mat" )
end

function print_performance_dia( nrg )
    performance_dia = nrg.p_dia

    println( "MATRIX DIAGONALIZATION PERFORMANCE" )
    println()

    times_dia = [ p.time for p in performance_dia ]
    t_min_dia = minimum( times_dia )
    t_max_dia = maximum( times_dia )
    t_mean_dia = sum(times_dia)/length(times_dia)
    println( "TIME (s)" )
    println( "minimum time: $t_min_dia" )
    println( "maximum time: $t_max_dia" )
    println( "mean time:    $t_mean_dia" )

    println()

    bytes_dia = [ p.bytes for p in performance_dia ]
    b_min_dia = minimum( bytes_dia )/10^9
    b_max_dia = maximum( bytes_dia )/10^9
    b_mean_dia = sum(bytes_dia)/length(bytes_dia)/10^9
    println( "minimum gigabytes: $b_min_dia" )
    println( "maximum gigabytes: $b_max_dia" )
    println( "mean gigabytes:    $b_mean_dia" )

    println()

    gc_dia = [ p.gctime for p in performance_dia ]
    gc_min_dia = minimum( gc_dia )
    gc_max_dia = maximum( gc_dia )
    gc_mean_dia = sum(gc_dia)/length(gc_dia)
    println( "minimum gc time: $gc_min_dia" )
    println( "maximum gc time: $gc_max_dia" )
    println( "mean gc time:    $gc_mean_dia" )
end

# =============================== #
# ASYMMETRIC ORBITAL CONSTRUCTION #
# =============================== #
