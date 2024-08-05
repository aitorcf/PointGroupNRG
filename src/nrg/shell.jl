using Printf
using Interpolations

# ###########################################
# SHELL CREATION OPERATORS
#
# - Used for computing pseudo-CG coefficients
#
# ...........................................

function basis2coperators( basis::CB ; 
                           n::Int64=-1 ) where {CB<:CanonicalBasis}
    #
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
    coperators::Vector{Operator{CB}} = []
    for sfs in basis.states
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

function shell_coperators( shell_basis::CB , 
                           hiztegia::D )::Dict{ ClearQNums , Operator{CB} } where {CB<:CanonicalBasis,D<:Dict}
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
    coperators = Dict{ ClearQNums , Operator{CB} }()
    @inbounds for (i,sfs) in enumerate(shell_basis.states)
        sum(sfs.occ)!==1 && continue
        idx = findfirst( sfs.occ )
        tup = sfs.hilbert.states[idx] 
        # sym = ( N , I , S , mu , m , r )
        sym::ClearQNums = ( convert(Int64,1) , 
                            convert(String,hiztegia[tup[2]]) ,
                            convert(Float64,tup[4]=="-" ? 0.0 : 0.5) ,
                            convert(Int64,tup[3]) ,
                            convert(Float64,hiztegia[tup[4]]) ,
                            convert(Int64,tup[1]) )
        cop::Operator{CB} = sfs2coperator(sfs,shell_basis)
        sym2cop::Dict{ ClearQNums , Operator{CB} } = Dict{ClearQNums,Operator{CB}}( sym => cop )
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
                symstates_shell::Dict{ClearQNums,S} , 
                basis_shell::CB , 
                hiztegia::D , 
                oirreps2dimensions::Dict{String,Int64} ; 
                verbose=false )::IntQPCG where {CB<:CanonicalBasis,S<:State,D<:Dict}
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
    
    pseudoCG_clear::ClearQPCG = ClearQPCG()

    # shell creation operators
    shell_cops::Dict{ ClearQNums , Operator{CB} } = shell_coperators( basis_shell , hiztegia )

    for (q_nu::ClearQNums,s_nu::S) in symstates_shell, 
        (q_mu::ClearQNums,s_mu::S) in symstates_shell

        for (q_a::ClearQNums,c_a::Operator{CB}) in shell_cops 
            q::ClearTripleQ = ( q_nu , q_a , q_mu ) 
            # ( q_nu | c^\dagger_{q_a} | q_mu )
            # notice that it is not the complex conjugate!!!
            coeff::ComplexF64 = (s_nu * c_a * s_mu)::ComplexF64
            isapprox( abs(coeff) , 0.0 , atol=1e-6 ) || push!( pseudoCG_clear , q=>coeff ) 
        end

    end

    pseudoCG_int::IntQPCG = IntQPCG( (convert_to_int(k[1],oirreps2dimensions),
                                      convert_to_int(k[2],oirreps2dimensions),
                                      convert_to_int(k[3],oirreps2dimensions))=>v 
                                      for (k::ClearTripleQ,v::ComplexF64) in pseudoCG_clear )

    if verbose 
        for (q,coeff) in pseudoCG_int
            println( "q = $q" )
            println( "coeff = $coeff" )
            println()
        end
    end

    return pseudoCG_int::IntQPCG
end

function get_pseudoCG_mat( 
            pseudoCG::Dict{Tuple{NTuple{6, Int64}, NTuple{6, Int64}, NTuple{6, Int64}}, ComplexF64} ,
            irrmult_0::Set{Tuple{Tuple{Int64, String, Float64}, Int64}} ,
            oindex2dimensions::Vector{Int64} ,
            oirreps2indices::Dict{String,Int64} )
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
#function get_pseudoCG( 
#            symstates_shell::Dict , 
#            basis_shell::CanonicalBasis , 
#            hiztegia::Dict ; 
#            verbose=false )
#    #
#    # compute the pseudo-CG coefficients for the shells.
#    #
#    # input:
#    #
#    # - symstates_shell : all symstates of a single shell (any)
#    # - basis_shell : basis for only that shell 
#    # - hiztegia : dict( irrep name => standard notation )
#    # - [format : "int" or "standard"]
#    #
#    # output:
#    #
#    # - pseudo-CG coefficients as dict( (q_nu,q_a,q_mu)=>coeff )
#    #
#    pseudoCG = Dict()
#
#    # shell creation operators
#    shell_cops = shell_coperators( basis_shell , hiztegia )
#
#    for (q_nu,s_nu) in symstates_shell, (q_mu,s_mu) in symstates_shell
#        for (q_a,c_a) in shell_cops 
#            q = ( q_nu , q_a , q_mu ) 
#            # ( q_nu | c^\dagger_{q_a} | q_mu )
#            # notice that it is not the complex conjugate!!!
#            coeff = s_nu * c_a * s_mu 
#            isapprox( coeff , 0 ) || push!( pseudoCG , q=>coeff ) 
#        end
#    end
#
#    if verbose 
#        for (q,coeff) in pseudoCG 
#            println( "q = $q" )
#            println( "coeff = $coeff" )
#            println()
#        end
#    end
#    return pseudoCG 
#end


# ################################
# MULTIPLETS u IN COMBINATION i,mu
# ................................
function get_combination_multiplets( 
            multiplets_block::Set{NTuple{4,Int64}} , 
            multiplets_shell::Set{NTuple{4,Int64}} ,
            keys_as_dict_o::Dict{NTuple{2,Int64},Vector{Int64}} ,
            keys_as_dict_s::Dict{NTuple{2,Int64},Vector{Int64}} ;
            verbose=false )

    combinations_u::Dict{IntIrrep,Vector{NTuple{3,IntMultiplet}}} = Dict{IntIrrep,Vector{NTuple{3,IntMultiplet}}}()

    G2R::Dict{ NTuple{3,Int64} , Int64 } = Dict{ NTuple{3,Int64} , Int64 }()

    m_u::IntMultiplet = (0,0,0,0)
    
    for m_i::IntMultiplet  in multiplets_block, 
        m_mu::IntMultiplet in multiplets_shell 

        if verbose 
            println( "m_mu = $m_mu ; m_i = $m_i" )
            println( "============================" )
        end
        
        (N_i::Int64,I_i::Int64,S_i::Int64,r_i) = m_i
        (N_mu::Int64,I_mu::Int64,S_mu::Int64,r_mu) = m_mu 

        N_u::Int64 = N_i + N_mu
        @views II_u::Vector{Int64} = keys_as_dict_o[(I_mu,I_i)]
        @views SS_u::Vector{Int64} = keys_as_dict_s[(S_mu,S_i)]

        for I_u::Int64 in II_u, 
            S_u::Int64 in SS_u

            G_u::IntIrrep = (N_u,I_u,S_u)
            if G_u in keys(G2R)
                G2R[G_u] += 1 
                m_u = ( G_u... , G2R[G_u] )
                push!( combinations_u[G_u] , (m_u,m_mu,m_i)::NTuple{3,IntMultiplet} )
            else 
                G2R[G_u] = 1
                m_u = ( G_u... , G2R[G_u] )
                combinations_u[G_u] = [(m_u,m_mu,m_i)]
            end
            #mergewith!( (x,y)->append!(x,y) , 
            #             combinations_u , 
            #             Dict{IntIrrep,Vector{NTuple{3,IntMultiplet}}}(G_u=>[(m_u,m_mu,m_i)]) )
            verbose && println( m_u )
        end
        verbose && println()
    end
    return combinations_u 
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
# CUTOFF 
# 
# - Introduce energy/state cutoff 
#   in irrEU 
# ...............................

function normalize_irrEU( irrEU )
    minE = minimum([e for (G,(E,U)) in irrEU for e in E])
    return Dict( G=>(E.-minE,U) for (G,(E,U)) in irrEU ) 
end

function cut_off!( 
            irrEU::Dict{ Tuple{Int64,Int64,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} };
            type::String="multiplet" , 
            cutoff::T=200 , 
            safeguard::Bool=true , 
            safeguard_tol::Float64=1e-2 ,
            safeguard_max::Int64=200 ,
            minmult::Int64=0 , 
            mine::Float64=-1.0 ,
            verbose::Bool=true ,
            M::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} }=Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} }()) where {T<:Number}
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
    mm::Vector{Tuple{IntMultiplet,Float64}} = collect( ((G...,r),E[r]) for (G,(E,U)) in irrEU for r=1:length(E) )
    sort!( mm , by=x->x[2] )

    # discard, kept and safeguarded
    discarded::Vector{IntMultiplet} = IntMultiplet[]
    kept::Vector{IntMultiplet} = IntMultiplet[]
    sg::Vector{IntMultiplet} = IntMultiplet[]

    # multiplet cutoff
    if type=="multiplet"

        # excessive cutoff, no discard
        if cutoff>=length(mm)

            kept = map( x->x[1] , mm )
            discarded = []

        # discard needed
        else

            mine_idx::Int64 = mm[end][2]<mine ? cutoff : findfirst( x->x[2]>=mine , mm )
            cutoff = maximum([ mine_idx , cutoff ])
            cutoff_energy = mm[cutoff][2]
            #kept = map( x->x[1] , mm[1:cutoff] )

            mm_kept = filter( x->(x[2]<=cutoff_energy || (safeguard && isapprox(x[2],cutoff_energy;atol=safeguard_tol))) , mm ) 
            length(mm_kept)>(cutoff+safeguard_max) && (mm_kept=mm_kept[1:(cutoff+safeguard_max)])
            mm_discarded = mm[(length(mm_kept)+1):end]

            kept      = map( x->x[1] , mm_kept )
            discarded = map( x->x[1] , mm_discarded )
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
        error( "type must be 'multiplet' or 'energy'" )
    end

    # apply cutoff to irrEU
    for (G,(E,U)) in irrEU 
        N = length(collect( k for k in kept if k[1:3]==G ))
        if N==0 
            pop!( irrEU , G )
            continue
        end
        # U does not need cutoff (multiplets run over)
        irrEU[G] = ( E[1:N] , U[:,1:N] )
    end

    if length(M)!==0
        for ((G_u,G_a,G_v),uavmat) in M 
            if !(haskey(irrEU,G_u) && haskey(irrEU,G_v)) 
                pop!( M , (G_u,G_a,G_v) )
                continue
            end
            R_u = length(irrEU[G_u][1]) 
            R_v = length(irrEU[G_v][1]) 
            M[G_u,G_a,G_v] = uavmat[1:R_u,:,1:R_v]
        end
    end
    return ( Set(kept) , discarded )::Tuple{Set{NTuple{4,Int64}},Vector{NTuple{4,Int64}}}
end


# #############
# MISCELLANEOUS
# .............

function print_symstates_ordered( symstates ) 
    K = sort( collect(keys(symstates)) ) 
    for k in K 
        println( k )
        println( symstates[k] )
        println()
    end
end

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

function convert_to_int( multiplet::ClearMultiplet , oirreps2indices::Dict{String,Int64} )::IntMultiplet
    ( N , I , S , r ) = multiplet
    return ( N , oirreps2indices[I] , Int64(2*S) , r )
end
function convert_to_int( irrep::ClearIrrep , oirreps2indices::Dict{String,Int64} )
    ( N , I , S ) = irrep 
    return ( N , oirreps2indices[I] , Int64(2*S) )
end
function convert_to_int( q::ClearQNums , oirreps2indices::Dict{String,Int64} )
    ( N , I , S , i , s , r ) = q
    Ii = oirreps2indices[I]
    Si = Int64(2*S)
    si = Int64(2*s)
    return (N,Ii,Si,i,si,r)
end
function convert_to_int( s::String , oirreps2indices::Dict{String,Int64} )
    return oirreps2indices[s]
end
function convert_to_int( i::Int64 , oirreps2indices::Dict{String,Int64} )
    return i
end
function convert_to_int( f::Float64 , oirreps2indices::Dict{String,Int64} )
    return Int64(2*f)
end


# ============== #
# NRG ITERATIONS #
# ============== #

function NRG( label::String ,
              calculation::String ,
              iterations::Int64, 
              cutoff_type::String, 
              cutoff_magnitude::Number,
              L::Float64,
              hop_symparams::Dict{ Int64 , Matrix{ComplexF64} },
              irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
              multiplets_shell::Set{NTuple{4,Int64}}, 
              cg_o_fullmatint::Dict{Tuple{Int64, Int64, Int64}, Array{ComplexF64, 3}},
              cg_s_fullmatint::Dict{Tuple{Int64, Int64, Int64}, Array{ComplexF64, 3}},
              keys_as_dict_o::Dict{NTuple{2,Int64},Vector{Int64}} ,
              keys_as_dict_s::Dict{NTuple{2,Int64},Vector{Int64}} ,
              Csum_o_array::Array{ComplexF64,6} ,
              Csum_s_array::Array{ComplexF64,6} ,
              Bsum_o_array::Array{ComplexF64,6} ,
              Bsum_s_array::Array{ComplexF64,6} ,
              pcgred_shell::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} },
              multiplets_a::Vector{NTuple{4,Int64}} , 
              combinations_uprima::Dict{NTuple{3,Int64}, Vector{NTuple{3,NTuple{4,Int64}}}},
              betabar::Float64 ,
              oindex2dimensions::Vector{Int64} ,
              channels_codiagonals::Vector{Dict{Int64,Vector{Float64}}},
              max_spin2::Int64;
              verbose::Bool=false ,
              distributed::Bool=false ,
              minmult::Int64=0 , 
              mine::Float64=0.0 ,
              z::Float64=0.0 ,
              spectral::Bool=false ,
              spectral_functions::Dict{String,Dict{IntMultiplet,Matrix{Float64}}}=Dict{String,Dict{IntMultiplet,Matrix{Float64}}}() ,
              K_factor::Float64=2.0 ,
              spectral_method::String="sakai1989",
              spectral_broadening::Float64=1.0 ,
              broadening_distribution::String="gauss",
              orbitalresolved::Bool=false ,
              dmnrg::Bool=false ,
              dmnrg_run::Int64=1 ,
              shell_dimension::Int64=0 ,
              density_matrices::Vector{Dict{IntIrrep,Matrix{Float64}}}=Dict{IntIrrep,Matrix{Float64}}[] ,
              half_weight_idx::Int64=0 ,
              half_weight_energy::Float64=0.0 ,
              M::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} }=Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} }() ,
              AA::Vector{T}=[] ,
              Karray_orbital::Array{ComplexF64,6}=Array{ComplexF64,6}(undef,0,0,0,0,0,0) , 
              Karray_spin::Array{ComplexF64,6}=Array{ComplexF64,6}(undef,0,0,0,0,0,0) ,
              multiplets_atomhop::Vector{NTuple{4,Int64}}=NTuple{4,Int64}[] ,
              scale::Float64=1.0 ,
              precompute_iaj::Bool=true ,
              compute_impmults::Bool=false ,
              mult2index::Dict{ClearMultiplet,Int64}=Dict{ClearMultiplet,Int64}() ,
              orbital_multiplets::Vector{ClearMultiplet}=ClearMultiplet[] ,
              mm_i::Dict{NTuple{4,Int64},Vector{Float64}}=Dict{NTuple{4,Int64},Vector{Float64}}() ,
              write_spectrum::Bool=false ,
              channels_diagonals::Vector{Dict{IntIrrep,Vector{Float64}}}=[] ,
              compute_selfenergy::Bool=false ,
              Mred_se::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} }=Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} }() ,
              AA_se=[] ,
              impurity_operators::Dict{String,Dict{IntTripleG,Array{ComplexF64,3}}}=Dict{String,Dict{IntTripleG,Array{ComplexF64,3}}}() ,
              spectral_temperature::Float64=0.0 ,
              extra_iterations::Int64=0 ) where {T}

    println( "=============" )
    println( "NRG PROCEDURE" )
    println( "=============" )
    println()

    # thermodynamic information from even and odd iterations
    thermo_even = zeros( Float64 , 0 , 8 )
    thermo_odd  = zeros( Float64 , 0 , 8 )
    thermo_evenodd  = zeros( Float64 , 0 , 8 )

    impspins       = []
    impnums        = []
    impmults       = []

    # spectrum after diagonalization in each step
    spectrum_even = Vector{Dict{IntIrrep,Vector{Float64}}}()
    spectrum_odd  = Vector{Dict{IntIrrep,Vector{Float64}}}()

    # performance at each step
    diagonalization_performances = []
    if spectral
        spectral_performances = []
    end

    # maximum energies and 2S (twice total spin)
    cutoff_eigenenergies = Float64[]
    cutoff_eigenenergy_all_iterations = 0.0
    maximum_spin2s = Int64[]
    maximum_spin2_all_iterations = 0

    # kept state ratios
    kept_discarded_ratios = Float64[]

    # create spectral directory
    isdir("spectral") || mkdir("spectral")

    # dmnrg
    diagonalizers::Vector{Dict{IntIrrep,Matrix{ComplexF64}}} = Dict{IntIrrep,Matrix{ComplexF64}}[]
    combinations_uprima_train = []
    multiplets_kept_train::Vector{Set{IntMultiplet}} = Set{IntMultiplet}[]
    multiplets_disc_train::Vector{Set{IntMultiplet}} = Set{IntMultiplet}[]
    spectrum_train = Dict{IntIrrep,Vector{Float64}}[]
    partition_dmnrg::Float64 = 0.0

    # NRG iterations
    nrg_performance = @timed for n in 2:iterations

        println( "-"^17 )
        @printf "ITERATION n = %3i\n" n
        println( "-"^17 )

        # cutoff
        #
        # apply cutoff
        irrEU_uncut = copy(irrEU)
        (multiplets_block, discarded) = cut_off!( irrEU ; 
                                                  type=cutoff_type , 
                                                  cutoff=cutoff_magnitude , 
                                                  safeguard=true ,
                                                  minmult=minmult ,
                                                  mine=mine ,
                                                  verbose=false ,
                                                  M=!iszero(length(impurity_operators)) ? impurity_operators["particle"] : Dict{IntTripleG,Array{ComplexF64,3}}() )

        # print info
        number_kept_multiplets = length(multiplets_block)
        number_kept_states = sum( (oindex2dimensions[m[2]]*(m[3]+1)) for m in multiplets_block )
        number_discarded_multiplets = length(discarded)
        number_total_multiplets = number_kept_multiplets+number_discarded_multiplets
        kept_discarded_ratio = number_kept_multiplets/number_total_multiplets
        cutoff_eigenenergy = maximum(collect( e for (G,(E,U)) in irrEU for e in E ))
        cutoff_eigenenergy_all_iterations = maximum([ cutoff_eigenenergy , cutoff_eigenenergy_all_iterations ])
        push!( cutoff_eigenenergies , cutoff_eigenenergy )
        println( "CUTOFF" )
        println( "  kept:      $number_kept_multiplets multiplets ($number_kept_states states)" )
        println( "  discarded: $number_discarded_multiplets multiplets" )
        println( "  ratio kept/total: $kept_discarded_ratio, $(number_kept_multiplets)/$(number_total_multiplets) " )
        println( "  cutoff eigenenergy: $cutoff_eigenenergy" )
        push!( kept_discarded_ratios , kept_discarded_ratio )

        if dmnrg && dmnrg_run==1 && !iszero(spectral_temperature) && !iszero(length(discarded))
            nn = n-1
            partition_dmnrg += sum( 
            oindex2dimensions[I]*(S+1)*exp(-irrEU_uncut[N,I,S][1][r]*iterscale(scale,L,nn)/spectral_temperature)*shell_dimension^(iterations-nn)
                for (N,I,S,r) in discarded
            )
        end

        # renormalize by √Λ
        for (G,(E,U)) in irrEU 
            irrEU[G] = ( E.*sqrt(L) , U )
        end

        # hopping parameter
        hop_symparams = Dict( 
            I_o => diagm(ComplexF64.(channels_codiagonals[n-1][I_o])) 
            for I_o in keys(channels_codiagonals[1])
        )
        println( "CONDUCTION COUPLING TERMS" )
        println( "Codiagonals (hoppings)")
        @printf "  %-8s  %-10s  %-s\n" "orbital" "multiplet" "amplitude"
        for (I,hop_matrix) in hop_symparams,
            cartesian_index in CartesianIndices(hop_matrix)

            i,j = Tuple(cartesian_index)
            isapprox(hop_matrix[i,j],0.0) && continue

            @printf "  %-8s  %-3i => %-3i  %-.3f\n" I i j hop_matrix[i,j]

        end
        println( "Diagonals (occupation energies)")
        @printf "  %-8s  %-10s  %-s\n" "irrep" "multiplet" "amplitude"
        for (I_o,I_o_multiplets_diagonals) in channels_diagonals[n],
            (r_o,r_o_multiplet_diagonal) in enumerate(I_o_multiplets_diagonals)

            @printf "  %-8s  %-3i => %-.3f\n" I_o r_o r_o_multiplet_diagonal

        end

        # diagonalization
        # 
        # construct and diagonalize ( m_u | H_1 | m_v )
        diagonalization_performance = @timed (irrEU,combinations_uprima) = matdiag_redmat( 
                multiplets_block , 
                multiplets_shell ,
                irrEU , 
                hop_symparams , 
                keys_as_dict_o ,
                keys_as_dict_s ,
                Csum_o_array ,
                Csum_s_array ,
                Bsum_o_array ,
                Bsum_s_array ,
                pcgred_shell ,
                pcgred_shell ,
                multiplets_a , 
                multiplets_a ,
                combinations_uprima ;
                verbose=verbose ,
                distributed=distributed ,
                precompute_iaj=precompute_iaj ,
                conduction_diagonals=channels_diagonals[n] );

        # information
        maximum_spin2 = maximum(collect( G[3] for (G,(E,U)) in irrEU ))
        if maximum_spin2>max_spin2
            error( """
                ERROR: max_spin2 too low!
                max_spin2 = $max_spin2
                maximum spin2 found in calculation = $maximum_spin2
            """)
        end
        maximum_spin2_all_iterations = maximum([ maximum_spin2 , maximum_spin2_all_iterations ])
        maximum_eigenenergy_after_diagonalization = maximum([ maximum(E) for (_,(E,_)) in irrEU ])
        println( "DIAGONALIZATION" )
        println( "  time: $(diagonalization_performance.time)s" )
        println( "  memory: $(diagonalization_performance.bytes*10^-6)Mb" )
        println( "  garbage collection time: $(diagonalization_performance.gctime)s" )
        println( "  maximum 2S needed: $maximum_spin2")
        println( "  maximum eigenenergy obtained: $maximum_eigenenergy_after_diagonalization" )

        # store performance
        push!( diagonalization_performances , diagonalization_performance )

        # dmnrg information
        if dmnrg && dmnrg_run==1

            # info gathered from cutoff to previous state
            push!( multiplets_kept_train , multiplets_block )
            push!( multiplets_disc_train , Set(discarded) )

            # info gathered from diagonalization in this step
            push!( combinations_uprima_train , combinations_uprima )
            push!( diagonalizers , Dict( G=>U for (G,(_,U)) in irrEU ) )
            push!( spectrum_train , Dict( G=>E for (G,(E,_)) in irrEU ) )
        end

        # spectrum 
        if write_spectrum
            (n%2==0)  && save_spectrum!( spectrum_even , irrEU )
            (n%2!==0) && save_spectrum!( spectrum_odd  , irrEU )
        end

        # impurity info
        if compute_impmults
            mm_i = imp_mults( irrEU ,
                              combinations_uprima ,
                              mm_i )
            m_imp::Vector{Float64} = mult_thermo( irrEU ,
                                 betabar ,
                                 oindex2dimensions ,
                                 mm_i )
            push!( impmults , m_imp )
        end

        # thermodynamics 
        #
        # iteration temperature
        temperature = compute_temperature_newdiscretization(n,L,betabar,scale)
        # compute thermodynamic variables
        magnetic_susceptibility = compute_magnetic_susceptibility( irrEU , betabar , oindex2dimensions )
        entropy = compute_entropy( irrEU , betabar , oindex2dimensions )
        heat_capacity = compute_heat_capacity( irrEU , betabar , oindex2dimensions )
        free_energy = compute_free_energy( irrEU , betabar , oindex2dimensions )
        number_particles = compute_average_number_of_particles( irrEU , betabar , oindex2dimensions )
        energy = compute_energy( irrEU , betabar , oindex2dimensions )
        partition_function = compute_partition_function( irrEU , betabar , oindex2dimensions )
        thermodynamic_matrix = [temperature magnetic_susceptibility entropy heat_capacity free_energy number_particles energy partition_function]
        # store even/odd therodynamic results
        if n%2==0
            thermo_even = vcat( thermo_even , thermodynamic_matrix )
        else
            thermo_odd  = vcat( thermo_odd , thermodynamic_matrix )
        end
        thermo_evenodd = vcat( thermo_evenodd , thermodynamic_matrix )
        # information
        println( "THERMODYNAMICS" )
        @printf "  %s = %.3e\n" "temperature" temperature
        @printf "  %s = %.3f\n" "magnetic susceptibility" magnetic_susceptibility
        @printf "  %s = %.3f\n" "entropy" entropy
        @printf "  %s = %.3f\n" "heat capacity" heat_capacity
        @printf "  %s = %.3f\n" "free energy" free_energy
        @printf "  %s = %i\n"   "average number of particles" number_particles
        @printf "  %s = %.3f\n" "energy" energy
        @printf "  %s = %.3e\n" "Z" partition_function

        # spectral function calculation
        if spectral 

            spectral_performance = @timed begin
                impurity_operators["particle"] = update_operator( impurity_operators["particle"], 
                                                                  collect(multiplets_atomhop) ,
                                                                  Karray_orbital ,
                                                                  Karray_spin ,
                                                                  combinations_uprima ,
                                                                  irrEU )
                multiplets_kept = Set{IntMultiplet}()
                multiplets_discarded = Set{IntMultiplet}()
                if dmnrg_run==2
                    irrEU_copy = copy(irrEU)
                    cut = n==iterations ? sum( (iszero(e) ? 1 : 0) for (_,(E,_)) in irrEU for e in E ) : cutoff_magnitude
                    (multiplets_kept,multiplets_discarded) = cut_off!( irrEU ; 
                                                                       type=cutoff_type , 
                                                                       cutoff=cut , 
                                                                       safeguard=(n!==iterations) ,
                                                                       minmult=minmult ,
                                                                       mine=mine )
                    irrEU = irrEU_copy
                end
                add_correlation_contribution!(
                    spectral_functions["spectral"],
                    impurity_operators["particle"],
                    impurity_operators["particle"],
                    oindex2dimensions,
                    irrEU,
                    n ,
                    broadening_distribution ,
                    spectral_broadening ,
                    iterscale(scale,L,n) ,
                    K_factor ; 
                    correlation_type="spectral",
                    T=spectral_temperature ,
                    limit_shell = n==iterations ,
                    extra_iterations=extra_iterations ,
                    density_matrix = dmnrg_run==2 ? density_matrices[n] : Dict{IntIrrep,Matrix{Float64}}() ,
                    multiplets_kept=collect(multiplets_kept) ,
                    multiplets_discarded=collect(multiplets_discarded) ,
                    L=L,
                    scale=scale,
                    half_weight_idx=half_weight_idx,
                    half_weight_energy=half_weight_energy
                )
                #compute_correlation_peaks(
                #    impurity_operators["particle"],
                #    impurity_operators["particle"],
                #    oindex2dimensions,
                #    irrEU,
                #    z,
                #    n;
                #    correlation_type="spectral",
                #    T=spectral_temperature,
                #    iteration_scale=iterscale(scale,L,n),
                #    use_density_matrix=(dmnrg && dmnrg_run==2)
                #)
            end


            #if orbitalresolved

            #    spectral_performance = @timed M, AA = update_redmat_AA_CGsummethod_orbitalresolved(
            #                M ,
            #                irrEU ,
            #                combinations_uprima ,
            #                collect(multiplets_atomhop) ,
            #                cg_o_fullmatint ,
            #                cg_s_fullmatint ,
            #                Karray_orbital ,
            #                Karray_spin ,
            #                AA ,
            #                oindex2dimensions ;
            #                verbose=false )

            #else 

            #    if !dmnrg

            #        spectral_performance = @timed M, AA = update_redmat_AA_CGsummethod(
            #                    M ,
            #                    irrEU ,
            #                    combinations_uprima ,
            #                    collect(multiplets_atomhop) ,
            #                    cg_o_fullmatint ,
            #                    cg_s_fullmatint ,
            #                    Karray_orbital ,
            #                    Karray_spin ,
            #                    AA ,
            #                    oindex2dimensions ;
            #                    verbose=false )

            #    elseif dmnrg && dmnrg_run==2

            #        irrEU_copy = deepcopy(irrEU)
            #        multiplets_kept,_ = cut_off!( irrEU ; 
            #                                      type=cutoff_type , 
            #                                      cutoff= n==iterations ? 1 : cutoff_magnitude , 
            #                                      safeguard= n<true ,
            #                                      minmult=minmult ,
            #                                      mine=mine )

            #        spectral_performance = @timed M, AA = update_redmat_AA_CGsummethod(
            #                    M ,
            #                    irrEU_copy ,
            #                    combinations_uprima ,
            #                    collect(multiplets_atomhop) ,
            #                    cg_o_fullmatint ,
            #                    cg_s_fullmatint ,
            #                    Karray_orbital ,
            #                    Karray_spin ,
            #                    AA ,
            #                    oindex2dimensions ;
            #                    verbose=false ,
            #                    density_matrix=density_matrices[n-1] ,
            #                    multiplets_kept=multiplets_kept )

            #    end


            #    if compute_selfenergy
            #        Mred_se, AA_se = update_selfenergy_CGsummethod(
            #                    M ,
            #                    Mred_se ,
            #                    irrEU ,
            #                    combinations_uprima ,
            #                    collect(multiplets_atomhop) ,
            #                    cg_o_fullmatint ,
            #                    cg_s_fullmatint ,
            #                    Karray_orbital ,
            #                    Karray_spin ,
            #                    AA_se ,
            #                    oindex2dimensions ;
            #                    verbose=false )

            #    end
            #end

            # information
            maximum_irrep_spin2 = irreps -> maximum((irreps[1][3],irreps[2][3],irreps[3][3]))
            maximum_spin2 = maximum([ maximum_irrep_spin2(irreps) for irreps in keys(impurity_operators["particle"]) ])
            maximum_spin2_all_iterations = maximum([ maximum_spin2_all_iterations , maximum_spin2 ])
            println( "EXCITATION MATRIX CALCULATION" )
            println( "  time: $(spectral_performance.time)s" )
            println( "  memory: $(spectral_performance.bytes*10^-6)Mb" )
            println( "  garbage collection time: $(spectral_performance.gctime)s" )
            println( "  maximum 2S needed: $maximum_spin2")
            push!( spectral_performances , spectral_performance )

        end

        println()
    end # NRG iterations finished
    println()


    # average even and odd thermodynamic results
    #
    # reverse thermodynamic matrices for interpolation
    reverse!( thermo_even , dims=1 )
    reverse!( thermo_odd  , dims=1 )
    reverse!( thermo_evenodd , dims=1 )
    # collect temperatures from even and odd
    temperatures_even = thermo_even[:,1]
    temperatures_odd  = thermo_odd[:,1]
    temperatures_evenodd = sort(vcat(temperatures_even,temperatures_odd))
    # interpolate even and odd data to new temperatures
    thermo_even_interpolated = interpolate_thermo_matrix( thermo_even , temperatures_evenodd )
    thermo_odd_interpolated  = interpolate_thermo_matrix( thermo_odd  , temperatures_evenodd )
    # average even and odd
    thermo_average = copy(thermo_even_interpolated)
    thermo_average[:,2:end] = 0.5*( thermo_even_interpolated[:,2:end] + thermo_odd_interpolated[:,2:end] )

    if spectral 

        save_correlation_spectral_decomposition(spectral_functions,label,z)
        #compute_correlation_spectral_decomposition(
        #    L ,
        #    iterations ,
        #    scale ,
        #    spectral_broadening ,
        #    label ,
        #    z ;
        #    broadening_distribution=broadening_distribution,
        #    K_factor=K_factor ,
        #    correlation_type="spectral" ,
        #    T=spectral_temperature ,
        #    betabar=betabar 
        #)
        #compute_spectral_function(
        #    dmnrg ? AA[1:end-1] : AA ,
        #    L ,
        #    dmnrg ? iterations-1 : iterations ,
        #    scale ;
        #    spectral_broadening=spectral_broadening,
        #    method=spectral_method,
        #    label=label ,
        #    z=z ,
        #    K_factor=K_factor ,
        #    orbitalresolved=orbitalresolved ,
        #    broadening_distribution=broadening_distribution
        #)

        #if compute_selfenergy
        #    compute_spectral_function(
        #        AA_se ,
        #        L ,
        #        iterations ,
        #        scale ;
        #        spectral_broadening=spectral_broadening,
        #        method=spectral_method,
        #        label=label ,
        #        z=z ,
        #        K_factor=K_factor ,
        #        orbitalresolved=orbitalresolved ,
        #        broadening_distribution=broadening_distribution,
        #        triple_excitation=true
        #    )
        #end
    end

    # reduced density matrices for all steps
    if dmnrg && dmnrg_run==1

        if iszero(spectral_temperature)

            Z = sum( (iszero(e) ? oindex2dimensions[I]*(S+1) : 0.0)
                     for ((_,I,S),(E,_)) in irrEU
                     for e in E )
            push!( 
                density_matrices ,
                Dict(
                    (N,I,S)=>diagm([ ( iszero(e) ? 1.0/Z : 0.0 ) for e in E ])
                    for ((N,I,S),(E,_)) in irrEU
                    if iszero(E[1])
                )
            )
            for i in reverse(eachindex(diagonalizers))
                diagonalizer = diagonalizers[i]
                combinations_uprima = combinations_uprima_train[i]
                multiplets_kept = multiplets_kept_train[i]
                insert!( 
                    density_matrices ,
                    1 ,
                    backwards_dm_reduced_T0(density_matrices[1],
                                            diagonalizer,
                                            combinations_uprima,
                                            oindex2dimensions) 
                )
            end
            #for (i,dm) in enumerate(density_matrices)
            #    println( "DM for iteration $(i) with (i+1)=$(i+1)" )
            #    @show keys(dm)
            #    trace = sum( tr(d)*oindex2dimensions[I]*(S+1) for ((_,I,S),d) in dm )
            #    tr_plus = sum(  (N>(i+1) ? tr(d) : 0.0) for ((N,_,_),d) in dm )
            #    tr_minus = sum( (N<(i+1) ? tr(d) : 0.0) for ((N,_,_),d) in dm )
            #    @show tr_plus,tr_minus
            #    @show trace
            #    println()
            #end
        else # T!==0

            # add last contribution to partition function
            partition_dmnrg += sum( 
                oindex2dimensions[I]*(S+1)*exp(-e*iterscale(scale,L,iterations)/spectral_temperature)
                for ((_,I,S),(E,_)) in irrEU
                for e in E
            )
            Z_last = sum( 
                oindex2dimensions[I]*(S+1)*exp(-e*iterscale(scale,L,iterations)/spectral_temperature)
                for ((_,I,S),(E,_)) in irrEU
                for e in E
            )
            Z = partition_dmnrg
            @show Z
            @show Z_last
            @show 

            # set first density matrix for n=iterations
            push!( 
                density_matrices ,
                Dict(
                    (N,I,S)=>diagm([ exp(-e*iterscale(scale,L,iterations)/spectral_temperature)/Z for e in E ])
                    for ((N,I,S),(E,_)) in irrEU
                )
            )

            # weights of density matrices for each iteration
            partial_dm_weight_train = Float64[]
            for i in eachindex(spectrum_train)
                n = i+1
                if n!==iterations
                    multiplets_disc = multiplets_disc_train[i+1]
                    if length(multiplets_disc)==0
                        push!(partial_dm_weight_train,0.0)
                        continue
                    end
                    spectrum = spectrum_train[i]
                    #@show multiplets_disc
                    #print_dict(spectrum)
                    push!( 
                        partial_dm_weight_train ,
                        shell_dimension^(iterations-n)/
                        Z*
                        sum(oindex2dimensions[I]*(S+1)*exp(-spectrum[N,I,S][r]*iterscale(scale,L,n)/spectral_temperature) for (N,I,S,r) in multiplets_disc)
                    )
                else
                    spectrum = spectrum_train[end]
                    push!( 
                        partial_dm_weight_train ,
                        sum(oindex2dimensions[I]*(S+1)*exp(-e*iterscale(scale,L,n)/spectral_temperature) for ((_,I,S),E) in spectrum for e in E)/Z
                    )
                end
            end
            @show partial_dm_weight_train
            max_weight, max_weight_idx = findmax(partial_dm_weight_train)
            _, half_weight_idx = findmin( abs2.(max_weight/2 .- partial_dm_weight_train[max_weight_idx:end]) )
            half_weight_idx += max_weight_idx
            half_weight = partial_dm_weight_train[half_weight_idx]
            @show max_weight, max_weight_idx
            @show half_weight, half_weight_idx
            @show iterscale(scale,L,max_weight_idx+1)
            @show iterscale(scale,L,half_weight_idx+1)
            @show spectral_temperature
            half_weight_idx = max_weight_idx
            half_weight_energy = iterscale(scale,L,half_weight_idx+1)

            for i in reverse(eachindex(diagonalizers))
                iteration = i+1
                diagonalizer = diagonalizers[i]
                combinations_uprima = combinations_uprima_train[i]
                multiplets_kept = multiplets_kept_train[i]
                multiplets_disc = multiplets_disc_train[i]
                spectrum = i-1<=0 ? nothing : spectrum_train[i-1]
                @show i-1,length(spectrum_train)
                insert!( 
                    density_matrices ,
                    1 ,
                    i-1<=0 ? Dict{IntIrrep,Matrix{ComplexF64}}() : backwards_dm_reduced_Tnonzero(
                        density_matrices[1],
                        diagonalizer,
                        combinations_uprima,
                        oindex2dimensions,
                        multiplets_kept,
                        multiplets_disc,
                        spectral_temperature,
                        Z,
                        shell_dimension,
                        iterations,
                        iteration,
                        iterscale(scale,L,iteration),
                        spectrum
                    ) 
                )
                #for (G,m) in density_matrices[1]
                #    @show G
                #    @show m[1,1]
                #    println()
                #end
            end
        end
    end

    # print summary information
    println( "=========================" )
    println( "SUMMARY OF NRG ITERATIONS" )
    println( "=========================" )
    # energy, kept states and spin2
    cutoff_eigenenergy_average = sum(cutoff_eigenenergies)/length(cutoff_eigenenergies)
    kept_discarded_ratio_average = sum(kept_discarded_ratios)/length(kept_discarded_ratios)
    println( "CHECK" )
    println( "  maximum 2S needed: $maximum_spin2_all_iterations" )
    println( "  average cutoff eigenenergy: $cutoff_eigenenergy_average" )
    println( "  average kept/discarded ratio: $kept_discarded_ratio_average" )
    println( "PERFORMANCE SUMMARY" )
    println( "  Diagonalization" )
    print_performance_summary( diagonalization_performances )
    if spectral
        println( "  Excitation matrix calculation" )
        print_performance_summary( spectral_performances )
    end
    println( "  Global" )
    println( "    time: $(nrg_performance.time)s" )
    println( "    memory: $(nrg_performance.bytes*10^-6)Mb" )
    println( "    garbage collection time: $(nrg_performance.gctime)s" )

    # save data to file (spectral is saved from within the function)
    #
    # thermodynamics

    # directory
    isdir("thermodata") || mkdir("thermodata")

    # thermodynamic data for this given value of z
    write_thermodata_onez( thermo_average , calculation , label , z )
    write_thermodata_onez( thermo_evenodd , calculation , label*"_evenodd" , z )

    # impurity contribution (diff)
    if calculation=="IMP"
        thermo_clean_filename = thermo_filename_one_z( label , "clean" , z )
        println()
        if length(glob(thermo_clean_filename))!==0 
            println( "Saving thermodynamic impurity contribution to $(thermo_filename_one_z(label,"diff",z))..." )
            write_thermodiff( label , z )
            write_thermodiff( label*"_evenodd" , z )
        end
    end

    # impurity properties 
    if compute_impmults

        println( "Saving impurity projections..." )
        # directory
        isdir("impurityprojections") || mkdir("impurityprojections")

        if calculation=="IMP" 
            write_impurity_info( impmults , orbital_multiplets , mult2index , label , z )
        end

    end
    # spectra at each step
    if write_spectrum
        write_nrg_spectra( spectrum_even , spectrum_odd )
    end

    return ( 
        thermo = thermo_average ,
        diagonalization_performances = diagonalization_performances ,
        impmults = compute_impmults ? impmults : nothing ,
        spectral_performances = spectral ? spectral_performances : nothing ,
        density_matrices = density_matrices ,
        half_weight_idx = half_weight_idx ,
        half_weight_energy = half_weight_energy
    )
end
# pcgred method
function NRG_old_discretization( label::String ,
              calculation::String ,
              iterations::Int64, 
              cutoff_type::String, 
              cutoff_magnitude::Number,
              L::Float64,
              hop_symparams::Dict{ Int64 , Matrix{ComplexF64} },
              irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
              multiplets_shell::Set{NTuple{4,Int64}}, 
              cg_o_fullmatint::Dict{Tuple{Int64, Int64, Int64}, Array{ComplexF64, 3}},
              cg_s_fullmatint::Dict{Tuple{Int64, Int64, Int64}, Array{ComplexF64, 3}},
              keys_as_dict_o::Dict{NTuple{2,Int64},Vector{Int64}} ,
              keys_as_dict_s::Dict{NTuple{2,Int64},Vector{Int64}} ,
              Csum_o_array::Array{ComplexF64,6} ,
              Csum_s_array::Array{ComplexF64,6} ,
              Bsum_o_array::Array{ComplexF64,6} ,
              Bsum_s_array::Array{ComplexF64,6} ,
              pcgred_shell::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} },
              multiplets_a::Vector{NTuple{4,Int64}} , 
              combinations_uprima::Dict{NTuple{3,Int64}, Vector{NTuple{3,NTuple{4,Int64}}}},
              betabar::Float64 ,
              oindex2dimensions::Vector{Int64} ,
              xi_symparams::Dict{ Int64 , Vector{Vector{ComplexF64}} } ;
              verbose::Bool=false ,
              distributed::Bool=false ,
              minmult::Int64=0 , 
              mine::Float64=0.0 ,
              z::Float64=0.0 ,
              discretization::String="standard" ,
              spectral::Bool=false ,
              K_factor::Float64=2.0 ,
              spectral_method::String="sakai1989",
              spectral_broadening::Float64=1.0 ,
              orbitalresolved::Bool=false ,
              M::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} }=Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} }() ,
              AA::Vector{T}=[] ,
              Karray_orbital::Array{ComplexF64,6}=Array{ComplexF64,6}(undef,0,0,0,0,0,0) , 
              Karray_spin::Array{ComplexF64,6}=Array{ComplexF64,6}(undef,0,0,0,0,0,0) ,
              multiplets_atomhop::Vector{NTuple{4,Int64}}=NTuple{4,Int64}[] ,
              alpha::Float64=1.0 ,
              eta::Function=x->1.0 ,
              precompute_iaj::Bool=true ,
              compute_impmults::Bool=false ,
              mult2index::Dict{ClearMultiplet,Int64}=Dict{ClearMultiplet,Int64}() ,
              orbital_multiplets::Vector{ClearMultiplet}=ClearMultiplet[] ,
              mm_i::Dict{NTuple{4,Int64},Vector{Float64}}=Dict{NTuple{4,Int64},Vector{Float64}}() ,
              write_spectrum::Bool=false ) where {T}

    println( "=============" )
    println( "NRG PROCEDURE" )
    println( "=============" )
    println()

    # thermodynamic information from even and odd iterations
    thermo_even = zeros( Float64 , 0 , 8 )
    thermo_odd  = zeros( Float64 , 0 , 8 )

    impspins       = []
    impnums        = []
    impmults       = []

    # spectrum after diagonalization in each step
    spectrum_even = Vector{Dict{IntIrrep,Vector{Float64}}}()
    spectrum_odd  = Vector{Dict{IntIrrep,Vector{Float64}}}()

    # performance at each step
    diagonalization_performances = []
    if spectral
        spectral_performances = []
    end

    # maximum energies and 2S (twice total spin)
    cutoff_eigenenergies = Float64[]
    cutoff_eigenenergy_all_iterations = 0.0
    maximum_spin2s = Int64[]
    maximum_spin2_all_iterations = 0
    
    # kept state ratios
    kept_discarded_ratios = Float64[]

    xi::Vector{Float64} = []
    if discretization=="lanczos" 
        _,_,xi = get_hoppings( iterations , L , z , eta )
    else
        xi = compute_xi_vector( iterations , z , L ; discretization=discretization )
    end

    # create spectral directory
    isdir("spectral") || mkdir("spectral")

    # NRG iterations
    nrg_performance = @timed for n in 2:iterations

        println( "-"^17 )
        @printf "ITERATION n = %3i\n" n
        println( "-"^17 )

        # cutoff
        #
        # apply cutoff
        (multiplets_block, discarded) = cut_off!( irrEU ; 
                                                  type=cutoff_type , 
                                                  cutoff=cutoff_magnitude , 
                                                  safeguard=true ,
                                                  minmult=minmult ,
                                                  mine=mine ,
                                                  verbose=false ,
                                                  M=M )
        # print info
        number_kept_multiplets = length(multiplets_block)
        number_kept_states = sum( (oindex2dimensions[m[2]]*(m[3]+1)) for m in multiplets_block )
        number_discarded_multiplets = length(discarded)
        number_total_multiplets = number_kept_multiplets+number_discarded_multiplets
        kept_discarded_ratio = number_kept_multiplets/number_total_multiplets
        cutoff_eigenenergy = maximum(collect( e for (G,(E,U)) in irrEU for e in E ))
        cutoff_eigenenergy_all_iterations = maximum([ cutoff_eigenenergy , cutoff_eigenenergy_all_iterations ])
        push!( cutoff_eigenenergies , cutoff_eigenenergy )
        println( "CUTOFF" )
        println( "  kept:      $number_kept_multiplets multiplets ($number_kept_states states)" )
        println( "  discarded: $number_discarded_multiplets multiplets" )
        println( "  ratio kept/total: $kept_discarded_ratio, $(number_kept_multiplets)/$(number_total_multiplets) " )
        println( "  cutoff eigenenergy: $cutoff_eigenenergy" )
        push!( kept_discarded_ratios , kept_discarded_ratio )


        # renormalize by √Λ
        for (G,(E,U)) in irrEU 
            irrEU[G] = ( E.*sqrt(L) , U )
        end

        # hopping parameter
        if discretization=="lanczos"
            hop_symparams = Dict( k=>diagm(ComplexF64.([xi_symparams[k][i][n-1]
                                                        for i in 1:length(xi_symparams[k])]))
                                  for (k,v) in xi_symparams )
        else
            hop_symparams = Dict( k=>ComplexF64.(xi[n-1]*Matrix(LinearAlgebra.I,size(v)...)) # at n=2, we want ξ[1]=ξ_0
                                  for (k,v) in hop_symparams )
        end
        println( "HOPPING" )
        @printf "  %-8s  %-10s  %-s\n" "orbital" "multiplets" "amplitude"
        for (I,hop_matrix) in hop_symparams,
            cartesian_index in CartesianIndices(hop_matrix)

            i,j = Tuple(cartesian_index)
            isapprox(hop_matrix[i,j],0.0) && continue

            @printf "  %-8s  %-3i => %-3i  %-.3f\n" I i j hop_matrix[i,j]

        end

        # diagonalization
        # 
        # construct and diagonalize ( m_u | H_1 | m_v )
        diagonalization_performance = @timed (irrEU,combinations_uprima) = matdiag_redmat_old_discretization( 
                multiplets_block , 
                multiplets_shell ,
                irrEU , 
                hop_symparams , 
                cg_o_fullmatint , 
                cg_s_fullmatint ,
                keys_as_dict_o ,
                keys_as_dict_s ,
                Csum_o_array ,
                Csum_s_array ,
                Bsum_o_array ,
                Bsum_s_array ,
                pcgred_shell ,
                pcgred_shell ,
                multiplets_a , 
                multiplets_a ,
                combinations_uprima , 
                oindex2dimensions ;
                verbose=verbose ,
                distributed=distributed ,
                precompute_iaj=precompute_iaj );
        # information
        maximum_spin2 = maximum(collect( G[3] for (G,(E,U)) in irrEU ))
        maximum_spin2_all_iterations = maximum([ maximum_spin2 , maximum_spin2_all_iterations ])
        maximum_eigenenergy_after_diagonalization = maximum([ maximum(E) for (_,(E,_)) in irrEU ])
        println( "DIAGONALIZATION" )
        println( "  time: $(diagonalization_performance.time)s" )
        println( "  memory: $(diagonalization_performance.bytes*10^-6)Mb" )
        println( "  garbage collection time: $(diagonalization_performance.gctime)s" )
        println( "  maximum 2S needed: $maximum_spin2")
        println( "  maximum eigenenergy obtained: $maximum_eigenenergy_after_diagonalization" )

        # store performance
        push!( diagonalization_performances , diagonalization_performance )

        # spectrum 
        if write_spectrum
            (n%2==0)  && save_spectrum!( spectrum_even , irrEU )
            (n%2!==0) && save_spectrum!( spectrum_odd  , irrEU )
        end

        # impurity info
        if compute_impmults
            mm_i = imp_mults( irrEU ,
                              oindex2dimensions ,
                              combinations_uprima ,
                              mm_i )
            m_imp::Vector{Float64} = mult_thermo( irrEU ,
                                 betabar ,
                                 oindex2dimensions ,
                                 mm_i )
            push!( impmults , m_imp )
        end

        # thermodynamics 
        #
        # iteration temperature
        if discretization!=="lanczos" 
            temperature = compute_temperature( n , L , betabar ; z=z , discretization=discretization )
        elseif discretization=="lanczos" 
            temperature = compute_temperature( n , L , betabar ; z=z , discretization=discretization , first_asymptotic_hopping_amplitude=ebar[1] )
        end
        # compute thermodynamic variables
        magnetic_susceptibility = compute_magnetic_susceptibility( irrEU , betabar , oindex2dimensions )
        entropy = compute_entropy( irrEU , betabar , oindex2dimensions )
        heat_capacity = compute_heat_capacity( irrEU , betabar , oindex2dimensions )
        free_energy = compute_free_energy( irrEU , betabar , oindex2dimensions )
        number_particles = compute_average_number_of_particles( irrEU , betabar , oindex2dimensions )
        energy = compute_energy( irrEU , betabar , oindex2dimensions )
        partition_function = compute_partition_function( irrEU , betabar , oindex2dimensions )
        thermodynamic_matrix = [temperature magnetic_susceptibility entropy heat_capacity free_energy number_particles energy partition_function]
        # store even/odd therodynamic results
        if n%2==0
            thermo_even = vcat( thermo_even , thermodynamic_matrix )
        else
            thermo_odd  = vcat( thermo_odd , thermodynamic_matrix )
        end
        # information
        println( "THERMODYNAMICS" )
        @printf "  %s = %.3e\n" "temperature" temperature
        @printf "  %s = %.3f\n" "magnetic susceptibility" magnetic_susceptibility
        @printf "  %s = %.3f\n" "entropy" entropy
        @printf "  %s = %.3f\n" "heat capacity" heat_capacity
        @printf "  %s = %.3f\n" "free energy" free_energy
        @printf "  %s = %i\n"   "average number of particles" number_particles
        @printf "  %s = %.3f\n" "energy" energy
        @printf "  %s = %.3e\n" "Z" partition_function

        # spectral function calculation
        if spectral 

            if orbitalresolved

                spectral_performance = @timed M, AA = update_redmat_AA_CGsummethod_orbitalresolved(
                            M ,
                            irrEU ,
                            combinations_uprima ,
                            collect(multiplets_atomhop) ,
                            cg_o_fullmatint ,
                            cg_s_fullmatint ,
                            Karray_orbital ,
                            Karray_spin ,
                            AA ,
                            oindex2dimensions ;
                            verbose=false )

            else 

                spectral_performance = @timed M, AA = update_redmat_AA_CGsummethod(
                            M ,
                            irrEU ,
                            combinations_uprima ,
                            collect(multiplets_atomhop) ,
                            cg_o_fullmatint ,
                            cg_s_fullmatint ,
                            Karray_orbital ,
                            Karray_spin ,
                            AA ,
                            oindex2dimensions ;
                            verbose=false )
            end

            # security check while in development
            for matrix in values(M)
                if any(isnan.(matrix))
                    error( "NaN in M")
                end
            end
            for weights in [x[2] for x in values(AA[end])]
                if any(isnan.(weights))
                    error( "NaN in A" )
                end
            end

            # information
            maximum_irrep_spin2 = irreps -> maximum((irreps[1][3],irreps[2][3],irreps[3][3]))
            maximum_spin2 = maximum([ maximum_irrep_spin2(irreps) for irreps in keys(M) ])
            maximum_spin2_all_iterations = maximum([ maximum_spin2_all_iterations , maximum_spin2 ])
            println( "EXCITATION MATRIX CALCULATION" )
            println( "  time: $(spectral_performance.time)s" )
            println( "  memory: $(spectral_performance.bytes*10^-6)Mb" )
            println( "  garbage collection time: $(spectral_performance.gctime)s" )
            println( "  maximum 2S needed: $maximum_spin2")
            push!( spectral_performances , spectral_performance )

        end

        println()
    end # NRG iterations finished
    println()


    # average even and odd thermodynamic results
    #
    # reverse thermodynamic matrices for interpolation
    reverse!( thermo_even , dims=1 )
    reverse!( thermo_odd  , dims=1 )
    # collect temperatures from even and odd
    temperatures_even = thermo_even[:,1]
    temperatures_odd  = thermo_odd[:,1]
    temperatures_evenodd = sort(vcat(temperatures_even,temperatures_odd))
    # interpolate even and odd data to new temperatures
    thermo_even_interpolated = interpolate_thermo_matrix( thermo_even , temperatures_evenodd )
    thermo_odd_interpolated  = interpolate_thermo_matrix( thermo_odd  , temperatures_evenodd )
    # average even and odd
    thermo_average = copy(thermo_even_interpolated)
    thermo_average[:,2:end] = 0.5*( thermo_even_interpolated[:,2:end] + thermo_odd_interpolated[:,2:end] )

    if spectral 

        compute_spectral_function(
            AA ,
            L ,
            iterations ,
            alpha ;
            spectral_broadening=spectral_broadening,
            method=spectral_method,
            label=label ,
            z=z ,
            K_factor=K_factor ,
            orbitalresolved=orbitalresolved
        )

    end

    # print summary information
    println( "=========================" )
    println( "SUMMARY OF NRG ITERATIONS" )
    println( "=========================" )
    # energy, kept states and spin2
    cutoff_eigenenergy_average = sum(cutoff_eigenenergies)/length(cutoff_eigenenergies)
    kept_discarded_ratio_average = sum(kept_discarded_ratios)/length(kept_discarded_ratios)
    println( "CHECK" )
    println( "  maximum 2S needed: $maximum_spin2_all_iterations" )
    println( "  average cutoff eigenenergy: $cutoff_eigenenergy_average" )
    println( "  average kept/discarded ratio: $kept_discarded_ratio_average" )
    println( "PERFORMANCE SUMMARY" )
    println( "  Diagonalization" )
    print_performance_summary( diagonalization_performances )
    if spectral
        println( "  Excitation matrix calculation" )
        print_performance_summary( spectral_performances )
    end
    println( "  Global" )
    println( "    time: $(nrg_performance.time)s" )
    println( "    memory: $(nrg_performance.bytes*10^-6)Mb" )
    println( "    garbage collection time: $(nrg_performance.gctime)s" )

    # save data to file (spectral is saved from within the function)
    #
    # thermodynamics
    if !spectral

        # directory
        isdir("thermodata") || mkdir("thermodata")

        # thermodynamic data for this given value of z
        write_thermodata_onez( thermo_average , calculation , label , z )

        # impurity contribution (diff)
        if calculation=="IMP"
            thermo_clean_filename = thermo_filename_one_z( label , "clean" , z )
            println()
            if length(glob(thermo_clean_filename))!==0 
                println( "Saving thermodynamic impurity contribution to $(thermo_filename_one_z(label,"diff",z))..." )
                write_thermodiff( label , z )
            end
        end

    end
    # impurity properties 
    if compute_impmults

        println( "Saving impurity projections..." )
        # directory
        isdir("impurityprojections") || mkdir("impurityprojections")

        if calculation=="IMP" 
            write_impurity_info( impmults , orbital_multiplets , mult2index , label , z )
        end

    end
    # spectra at each step
    if write_spectrum
        write_nrg_spectra( spectrum_even , spectrum_odd )
    end

    return ( 
        thermo = thermo_average ,
        diagonalization_performances = diagonalization_performances ,
        impmults = compute_impmults ? impmults : nothing ,
        spectral_performances = spectral ? spectral_performances : nothing
    )
end


# #######################
# PERFORMANCE INFORMATION
# .......................

function print_performance_summary( diagonalization_performances )

    times  = [ p.time for p in diagonalization_performances ]
    t_min  = minimum( times )
    t_max  = maximum( times )
    t_average = sum(times)/length(times)
    println( "    - time" )
    println( "        minimum: $(t_min)s" )
    println( "        maximum: $(t_max)s" )
    println( "        average: $(t_average)s" )

    bytes  = [ p.bytes for p in diagonalization_performances]
    b_min  = minimum( bytes )/1e6
    b_max  = maximum( bytes )/1e6
    b_average = sum(bytes)/length(bytes)/1e6
    println( "    - memory")
    println( "        minimum: $(b_min)Mb" )
    println( "        maximum: $(b_max)Mb" )
    println( "        average: $(b_average)Mb" )

    gc = [ p.gctime for p in diagonalization_performances]
    gc_min = minimum( gc )
    gc_max = maximum( gc )
    gc_average = sum(gc)/length(gc)
    println( "    - garbage collection time" )
    println( "        minimum: $(gc_min)s" )
    println( "        maximum: $(gc_max)s" )
    println( "        average: $(gc_average)s" )

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


# ~~~~~~~~~~~~~~~~~~~ #
# IMPURITY MULTIPLETS #
# ~~~~~~~~~~~~~~~~~~~ #
function ordered_multiplets( multiplets )
    max_N::Int64 = maximum(m[1] for m in multiplets) 
    omultiplets = []
    for N in 0:max_N,
        mult in multiplets

        mult[1]==N || continue
        push!(omultiplets,mult)

    end
    return omultiplets
end

# Specrtra from each iteration: Eigenvalues
function save_spectrum!( 
            spectrum::Vector{Dict{IntIrrep,Vector{Float64}}} , 
            irrEU::Dict{IntIrrep,Tuple{Vector{Float64},Matrix{ComplexF64}}} )

    push!(
        spectrum ,
        Dict(
            G=>E
            for (G,(E,U)) in irrEU
        )
    )

end
function write_iteration_spectrum(
        f::IO ,
        n::Int64 ,
        iteration_spectrum::Dict{IntIrrep,Vector{Float64}} )

    # n = NRG step
    write( f , "N = $(n)\n" )

    iteration_spectrum_sorted = sort(
        [ (K,V) for (K,V) in iteration_spectrum ],
        by = x->minimum(x[2])
    )

    # irrep iteration
    for ((N,I,S),E) in iteration_spectrum_sorted

        # line begins with irrep
        write( f , "$(N) $(I) $(S)" )

        # irrep spectrum
        for eigenenergy in E
            write( f , " $(eigenenergy)" )
        end

        write( f , "\n" )

    end

end
function write_nrg_spectra( 
            spectrum_even , 
            spectrum_odd )

    # even spectra
    open( "spectrum_even.dat" , write=true ) do f 

        for (i,iteration_spectrum) in enumerate(spectrum_even)

            # first dump is at n=2
            n = 2*i

            # write iteration spectrum
            write_iteration_spectrum( f , n , iteration_spectrum )

        end
    end

    # odd spectra
    open( "spectrum_odd.dat" , write=true ) do f 

        for (i,iteration_spectrum) in enumerate(spectrum_odd)

            # first dump is at n=3
            n = 2*i + 1

            # write iteration spectrum
            write_iteration_spectrum( f , n , iteration_spectrum )

        end
    end

end
