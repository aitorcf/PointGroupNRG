include( "shell.jl" )
using Base.Threads

function compute_step_matels_multithread( 
        multiplets_block::Set{NTuple{4,Int64}}, 
        multiplets_shell::Set{NTuple{4,Int64}},
        irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
        hop::Dict{ Tuple{Int64,Int64} , ComplexF64 },
        cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        pcg::Dict{ NTuple{3,NTuple{6,Int64}} , ComplexF64 },
        pcgmat::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,9} },
        qq_a::Set{NTuple{6,Int64}}, 
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
    1

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

    multiplets_mui2u_vec = collect( ((m_mu,m_i),mm_u) 
                                for ((m_mu,m_i),mm_u) in multiplets_mui2u )
    
    # iteration over multiplets in i,mu and j,nu, 
    # and their combinations u,v
    @threads for ((m_mu,m_i),mm_u) in multiplets_mui2u_vec
        for ((m_nu,m_j),mm_v) in multiplets_mui2u_vec

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
            #U_i = irrEU[G_i][2]
            #U_j = irrEU[G_j][2]

            #if (N_nu==(N_mu+1) && N_i==(N_j+1)) #hopping possible
            #    # (mu',i'=>u'), (nu',j'=>v') combinations 
            #    @inbounds begin
            #        combs_uprima_local = @view combinations_uprima[G_i][:]
            #        combs_vprima_local = @view combinations_uprima[G_j][:]
            #    end                          
            #end

            # iterate over irreps of multiplets u,v in 
            # combinations i,mu and j,nu
            for m_u::NTuple{4,Int64} in mm_u, 
                m_v::NTuple{4,Int64} in mm_v

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
                #hblockview = @view u_H1_v_blocks[G][:,:,:]

                # chosen partner: orbital 1, max z-spin
                (N,I,S,i,s) = ( G... , 1 , G[3] )
                verbose && println( "p_u = p_v = $((N,I,S,i,s))" )

                # clebsch-gordan matrices
                #@inbounds begin
                #    cgomat_muiu = @view cg_o_fullmatint[I_mu,I_i,I][:,:,i]
                #    cgomat_nujv = @view cg_o_fullmatint[I_nu,I_j,I][:,:,i]
                # end

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
                    cg_muiu = cg_o_fullmatint[I_mu,I_i,I][i_mu,i_i,i] *
                            PartialWaveFunctions.CG_doublearg(S_mu,s_mu,S_i,s_i,S,s)
                    cg_nujv = cg_o_fullmatint[I_nu,I_j,I][i_nu,i_j,i] *
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

                            # universal pseudo-CG coefficient (needs to be conjugated)
                            # MIGHT BE IMPROVED WITH PCGMAT 
                            nuamu = get( pcg , (q_nu,q_a,q_mu) , zero(ComplexF64) )
                            nuamu==zero(nuamu) && continue

                            # U-transformed coefficient
                            iaj = 0.0im
                            verbose && println( "u-transformed coefficients. iterating..." )
                            @inbounds for iu::Int64 in 1:length(combinations_uprima[G_i])    

                                # i -- U --> u' 
                                (m_up,m_mup,m_ip) = combinations_uprima[G_i][iu]::NTuple{3,NTuple{4,Int64}}
                                r_up = m_up[4]::Int64

                                q_up = ( q_i[1], q_i[2], q_i[3], q_i[4], q_i[5], r_up )
                                u_mu_i = ( q_up , m_mup , m_ip )

                                @inbounds for iv::Int64 in 1:length(combinations_uprima[G_j])

                                    # j -- U --> v'
                                    (m_vp,m_nup,m_jp) = combinations_uprima[G_j][iv]::NTuple{3,NTuple{4,Int64}}

                                    # hopping only changes last shell states
                                    m_ip==m_jp || continue

                                    r_vp = m_vp[4]::Int64

                                    verbose && println( "r_up=$r_up , r_vp=$r_vp" )

                                    # multiplication of U matrix elements
                                    uu = conj(irrEU[G_i][2][r_up,r_i]) * 
                                              irrEU[G_j][2][r_vp,r_j]
                                    uu==zero(uu) && continue
                                    verbose && println( "uu = $uu" )

                                    q_vp = ( q_j[1], q_j[2], q_j[3], q_j[4], q_j[5], r_vp )
                                    v_nu_j = ( q_vp , m_nup , m_jp )

                                    # pseudo-CG coefficient
                                    pseudo = get_pseudoCG_up_vp_multithread( 
                                                                u_mu_i,
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
                    lock( ReentrantLock() ) do
                        #u_H1_v_blocks[G][r_u,r_v,:] .+= cg_tot .* 
                        #                                ComplexF64[ E , hopart ]
                        u_H1_v_blocks[G][r_u,r_v,1] += cg_tot * E 
                        u_H1_v_blocks[G][r_u,r_v,2] += cg_tot * hopart 
                    end
                    #hblockview[r_u,r_v,:] .+= cg_tot .* SVector{2,ComplexF64}( E , hopart ) 
                    if any( isnan.(u_H1_v_blocks[G][r_u,r_v,1]))
                        @show u_H1_v_blocks[G][r_u,r_v,:]
                    end
                end
                verbose && println()
                verbose && println( "matel = $(u_H1_v_blocks[G][r_u,r_v,:])" )
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

@inline function get_pseudoCG_up_vp_multithread( 
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
    #pcg_local = @view pcgmat[G_mu,G_a,G_nu][:,:,:,:,r_mu,r_nu,i_a,si_a,r_a]

    # clebsch-gordan matrices
    #@inbounds begin
    #    cgomat_muiu = @view cg_o_fullmatint[I_mu,I_ij,I_u][:,:,i_u]
    #    cgomat_nujv = @view cg_o_fullmatint[I_nu,I_ij,I_v][:,:,i_v]
    #end

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

        s_i = s_j = s_ij        

        if verbose 
            println( "g_mu' = ( $i_mu , $(s_mu/2.0) )" )
            println( "g_nu' = ( $i_nu , $(s_nu/2.0) )" )
        end

        # coefficient
        si_mu = Int64((s_mu+S_mu)/2+1)
        si_nu = Int64((s_nu+S_nu)/2+1)
        #pcgcoeff = pcg_local[i_mu,si_mu,i_nu,si_nu]
        pcgcoeff = pcgmat[G_mu,G_a,G_nu][i_mu,si_mu,i_nu,si_nu,r_mu,r_nu,i_a,si_a,r_a]
        pcgcoeff==zero(pcgcoeff) && continue
        verbose && println( "pcgcoeff = $pcgcoeff")

        # cg coefficients
        cg_muiu = cg_o_fullmatint[I_mu,I_ij,I_u][i_mu,i_ij,i_u] *
                  PartialWaveFunctions.CG_doublearg(S_mu,s_mu,S_i,s_i,S_u,s_u)
        cg_nujv = cg_o_fullmatint[I_nu,I_ij,I_v][i_nu,i_ij,i_v] * 
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



#function compute_step_matels_multithread( 
#        multiplets_block::Set{NTuple{4,Int64}}, 
#        multiplets_shell::Set{NTuple{4,Int64}},
#        irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
#        hop::Dict{ Tuple{Int64,Int64} , ComplexF64 },
#        cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
#        pcg::Dict{ NTuple{3,NTuple{6,Int64}} , ComplexF64 },
#        pcgmat::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,9} },
#        qq_a::Set{NTuple{6,Int64}}, 
#        combinations_uprima::Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} },
#        oindex2dimensions::Vector{ Int64 };
#        verbose=false )
#    # computes matrix elements ( m_u | H_1 | m_v ) for 
#    # the next step.
#    # input:
#    # - multiplets_block : { q_i }
#    # - multiplets_shell : { q_mu }
#    # - irrEU : G => E, U (main result of previous step)
#    # - hop : m_a => xi for n-th step
#    # - cg_o_fullmatin : orbital CG coefficients for this problem (int format)
#    # - pcg : pseudo-CG coefficients
#    # - pcgmat : pseudo-CG coefficients in matrix form
#    # - qq_a : set of all q_a
#    # - combinations_uprima : m_u' => m_mu', m_i'
#    # - oirreps2dimensions : orbital_irrep => dimension 
#    # output:
#    # - u_H1_v : ( m_u | H_1 | m_v )
#    # - combinations_uprima_new : m_u => m_mu, m_i
#
#    # type annotations
#    m_up::Tuple{Int64,Int64,Int64,Int64}  = (0,0,0,0)
#    m_vp::Tuple{Int64,Int64,Int64,Int64}  = (0,0,0,0)
#    m_ip::Tuple{Int64,Int64,Int64,Int64}  = (0,0,0,0)
#    m_jp::Tuple{Int64,Int64,Int64,Int64}  = (0,0,0,0)
#    m_mup::Tuple{Int64,Int64,Int64,Int64} = (0,0,0,0)
#    m_nup::Tuple{Int64,Int64,Int64,Int64} = (0,0,0,0)
#    iu::Int64 = 0
#    iv::Int64 = 0
#    uu::ComplexF64      = 0
#    hoparam::ComplexF64 = 0
#    pseudo::ComplexF64  = 0
#    iaj::ComplexF64     = 0
#    
#
#    # block-shell combination multiplets 
#    multiplets_mui2u::Dict{NTuple{2,NTuple{4,Int64}},Vector{NTuple{4,Int64}}} = 
#            get_combination_multiplets( multiplets_block , 
#                                        multiplets_shell , 
#                                        cg_o_fullmatint ;
#                                        verbose=false )
#
#    # new block-shell combinations 
#    GG_u = Set( m_u[1:3] for mm_u in values(multiplets_mui2u) 
#                         for m_u in mm_u )
#    combinations_uprima_new = Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} }(
#                G => NTuple{3,NTuple{4,Int64}}[ (m_u,m_mu,m_i) 
#                                                for ((m_mu,m_i),mm_u) in multiplets_mui2u 
#                                                for m_u in mm_u
#                                                if m_u[1:3]==G 
#                                              ] 
#                for G in GG_u
#    )
#
#    # Hamiltonian matrix
#    newmults = Set{NTuple{4,Int64}}( m for mm in values(multiplets_mui2u) for m in mm )
#    newirrepmults = get_irreps( newmults ; multiplicity=true )
#    u_H1_v_blocks = Dict{ NTuple{3,Int64} , Array{ComplexF64,3} }(
#            G => SharedArray{ComplexF64}(R,R,2) for (G,R) in newirrepmults
#    )
#
#    #@everywhere oindex2dimensions = $oindex2dimensions
#    #@everywhere hop = $hop
#    #@everywhere cg_o_fullmatint = $cg_o_fullmatint
#    #@everywhere pcg = $pcg 
#    #@everywhere pcgmat = $pcgmat
#    #@everywhere qq_a = $qq_a
#    
#    # iteration over multiplets in i,mu and j,nu, 
#    # and their combinations u,v
#    for ((m_mu,m_i),mm_u) in multiplets_mui2u, 
#        ((m_nu,m_j),mm_v) in multiplets_mui2u
#
#        if verbose
#            println( "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" )
#            println( "m_mu = $m_mu, m_i = $m_i" )
#            println( "m_nu = $m_nu, m_j = $m_j" )
#            println( "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" )
#        end
#
#        # block and shell quantum numbers 
#        (N_i, I_i, S_i, r_i ) = m_i
#        (N_j, I_j, S_j, r_j ) = m_j
#        G_i = (N_i, I_i, S_i)
#        G_j = (N_j, I_j, S_j)
#        (N_mu,I_mu,S_mu,r_mu) = m_mu
#        (N_nu,I_nu,S_nu,r_nu) = m_nu
#        
#        #(N_i+N_mu)==(N_j+N_nu) || continue # block-diagonality in N
#        # SHORT-CIRCUIT + CONTINUE NOT WORKING
#        if (N_i+N_mu)==(N_j+N_nu)
#
#            # U matrices (not transposed for faster loop)
#            
#            U_i = @view $(irrEU[G_i][2])[:,:]
#            U_j = @view $(irrEU[G_j][2])[:,:]
#
#            hopallowed = false
#            if (N_nu==(N_mu+1) && N_i==(N_j+1)) #hopping possible
#                # (mu',i'=>u'), (nu',j'=>v') combinations 
#                hopallowed = true
#                @inbounds begin
#                    combs_uprima_local = @view $(combinations_uprima[G_i])[:]
#                    combs_vprima_local = @view $(combinations_uprima[G_j])[:]
#                end                          
#            end
#
#            # iterate over irreps of multiplets u,v in 
#            # combinations i,mu and j,nu
#            @sync for m_u::NTuple{4,Int64} in mm_u, 
#                      m_v::NTuple{4,Int64} in mm_v
#
#                @async hblock[m_u[4],m_v[4],:] .+= @fetch compute_blockelement( 
#                    m_u , m_v ,
#                    m_i  , m_j , 
#                    m_mu , m_nu )
#
#                r_u = m_u[4]
#                r_v = m_v[4]
#
#                if verbose 
#                    println( "..........................." )
#                    println( "( $m_u | H_1 | $m_v )" )
#                    println( "..........................." )
#                    println()
#                end
#                
#                # irreps and block-diagonality
#                G_u = m_u[1:3]
#                G_v = m_v[1:3]
#                #G_u==G_v || continue # block-diagonality in irrep
#                if G_u==G_v 
#
#                    G = G_u 
#                    verbose && println( "G_u = G_v = $G" )
#
#                    # select H block view
#                    hblock = @view u_H1_v_blocks[G][:,:,:]
#
#                    # chosen partner: orbital 1, max z-spin
#                    p = (N,I,S,i,s) = ( G... , 1 , G[3] )
#                    verbose && println( "p_u = p_v = $((N,I,S,i,s))" )
#
#
#                    # SPAWN HERE
#                    @async hblock[r_u,r_v,:] .= compute_blockelement( 
#                        m_u , m_v ,
#                        m_i  , m_j , 
#                        m_mu , m_nu )
#                end
#            end
#        end
#    end
#
#    # add h.c. 
#    u_H1_v_new = Dict{ NTuple{3,Int64} , Matrix{ComplexF64} }(
#        G => ComplexF64[ (sum(mat[r_u,r_v,:]) + conj(mat[r_v,r_u,2]))
#                         for r_u in 1:size(mat,1), r_v in 1:size(mat,2) ] 
#        for (G,mat) in u_H1_v_blocks
#    )
#
#    return ( u_H1_v_new, combinations_uprima_new )
#end
#
#function compute_blockelement( 
#        m_u::NTuple{4,Int64}  , m_v::NTuple{4,Int64} , 
#        m_i::NTuple{4,Int64}  , m_j::NTuple{4,Int64} , 
#        m_mu::NTuple{4,Int64} , m_nu::NTuple{4,Int64} , 
#        G_uv::NTuple{3,Int64} , r_u::Int64 , r_v::Int64 , 
#        verbose=false )
#
#    blockel = ComplexF64[ 0.0im , 0.0im ]
#
#    G_u = m_u[1:3]
#    G_v = m_v[1:3]
#    if G_u!==G_v 
#        return blockel
#    end
#    G = G_u
#    p = (N,I,S,i,s) = ( G... , 1 , G[3] )
#
#    r_u = m_u[4]
#    r_v = m_v[4]
#
#    (N_i,I_i,S_i,r_i) = m_i
#    (N_j,I_j,S_j,r_j) = m_j
#    (N_mu,I_mu,S_mu,r_mu) = m_mu
#    (N_nu,I_nu,S_nu,r_nu) = m_nu
#    (N,I,S) = G_uv
#
#    @inbounds begin
#        cgomat_muiu = @view cg_o_fullmatint[I_mu,I_i,I][:,:,i]
#        cgomat_nujv = @view cg_o_fullmatint[I_nu,I_j,I][:,:,i]
#    end
#
#    # clebsch-gordan series
#    verbose && println( "expanding CG series" )
#    @inbounds for i_j::Int64  = 1:(oindex2dimensions[I_j]::Int64),
#                  i_i::Int64  = 1:(oindex2dimensions[I_i]::Int64), 
#                  i_nu::Int64 = 1:(oindex2dimensions[I_nu]::Int64),
#                  i_mu::Int64 = 1:(oindex2dimensions[I_mu]::Int64), 
#                  s_j::Int64  = (-S_j):2:S_j,
#                  s_i::Int64  = (-S_i):2:S_i,
#                  s_nu::Int64 = (-S_nu):2:S_nu, 
#                  s_mu::Int64 = (-S_mu):2:S_mu
#        
#        # qnums
#        q_i  = (N_i,I_i,S_i,i_i,s_i,r_i)
#        q_j  = (N_j,I_j,S_j,i_j,s_j,r_j)
#        q_mu = (N_mu,I_mu,S_mu,i_mu,s_mu,r_mu)
#        q_nu = (N_nu,I_nu,S_nu,i_nu,s_nu,r_nu)
#        G_i  = (N_i,I_i,S_i)
#        G_j  = (N_j,I_j,S_j)
#
#        # cg coefficients
#        cg_muiu = cgomat_muiu[i_mu,i_i] * 
#                    PartialWaveFunctions.CG_doublearg(S_mu,s_mu,S_i,s_i,S,s)
#        cg_nujv = cgomat_nujv[i_nu,i_j] * 
#                    PartialWaveFunctions.CG_doublearg(S_nu,s_nu,S_j,s_j,S,s)
#        cg_tot  = conj( cg_muiu ) * cg_nujv
#        isapprox( cg_tot , 0.0im ) && continue
#
#        if verbose
#            println( "q_mu=$q_mu ; q_i=$q_i" )
#            println( "q_nu=$q_nu ; q_j=$q_j" )
#        end
#
#        # diagonal block energy
#        E = (q_i==q_j && q_mu==q_nu) ? convert(ComplexF64,E_Gi[r_i]) : 
#                                       zero(ComplexF64)
#        verbose && println( "E = $E" )
#
#        # hopping part iteration
#        hopart = 0.0im
#        
#        # hopping number matching
#        if (N_nu==(N_mu+1) && N_i==(N_j+1)) #hopping possible
#
#            verbose && println( "N_i=$N_i ; N_j=$N_j" )
#
#            verbose && println( "hopping iteration" )
#            for q_a in qq_a
#
#                verbose && println( "q_a = $q_a" )
#
#                # hopping parameter (orbital multiplet dependent)
#                hoparam = hop[(q_a[2],q_a[end])]
#
#                # universal pseudo-CG coefficient (needs to be conjugated)
#                # MIGHT BE IMPROVED WITH PCGMAT 
#                nuamu = get( pcg , (q_nu,q_a,q_mu) , zero(ComplexF64) )
#                nuamu==zero(nuamu) && continue
#
#                # OR SPAWN HERE
#
#                # U-transformed coefficient
#                iaj = 0.0im
#                verbose && println( "u-transformed coefficients. iterating..." )
#                @inbounds for iu::Int64 in 1:length(combs_uprima_local)    
#
#                    # i -- U --> u' 
#                    (m_up,m_mup,m_ip) = combs_uprima_local[iu]::NTuple{3,NTuple{4,Int64}}
#                    r_up = m_up[4]::Int64
#
#                    q_up = ( q_i[1], q_i[2], q_i[3], q_i[4], q_i[5], r_up )
#                    u_mu_i = ( q_up , m_mup , m_ip )
#
#                    @inbounds for iv::Int64 in 1:length(combs_vprima_local)
#
#                        # j -- U --> v'
#                        (m_vp,m_nup,m_jp) = combs_vprima_local[iv]::NTuple{3,NTuple{4,Int64}}
#
#                        # hopping only changes last shell states
#                        m_ip==m_jp || continue
#
#                        r_vp = m_vp[4]::Int64
#
#                        verbose && println( "r_up=$r_up , r_vp=$r_vp" )
#
#                        # multiplication of U matrix elements
#                        uu = conj(U_i[r_up,r_i]) * U_j[r_vp,r_j]
#                        uu==zero(uu) && continue
#                        verbose && println( "uu = $uu" )
#
#                        q_vp = ( q_j[1], q_j[2], q_j[3], q_j[4], q_j[5], r_vp )
#                        v_nu_j = ( q_vp , m_nup , m_jp )
#
#                        # pseudo-CG coefficient
#                        pseudo = get_pseudoCG_up_vp_distributed( 
#                                                    u_mu_i,
#                                                    v_nu_j,
#                                                    q_a,
#                                                    oindex2dimensions ,
#                                                    pcg , 
#                                                    pcgmat ,
#                                                    cg_o_fullmatint ;
#                                                    verbose=false )::ComplexF64
#                        verbose && println( "pseudo = $pseudo" )
#
#                        iaj += hoparam * uu * pseudo
#                        # Uncomment for testing.
#                        verbose && println( "-> iaj += $(hoparam * uu * pseudo)" )
#                    end
#                end
#                verbose && println( "q_a contribution (no sign) = $iaj" )
#                hopart += iaj * conj(nuamu)
#            end
#            hopart *= (-1)^N_mu
#        end
#
#        verbose && println( "hopart = $hopart" )
#
#        # putting everything together
#        #hblock[r_u,r_v,:] .+= cg_tot .* SVector{2,ComplexF64}( E , hopart ) 
#        blockel .+= cg_tot .* SVector{2,ComplexF64}( E , hopart ) 
#    end
#    return blockel
#end
#
#function get_pseudoCG_up_vp_multithread( 
#        u_mu_i::Tuple{NTuple{6,Int64},NTuple{4,Int64},NTuple{4,Int64}} , 
#        v_nu_j::Tuple{NTuple{6,Int64},NTuple{4,Int64},NTuple{4,Int64}} , 
#        q_a::NTuple{6,Int64} , 
#        oindex2dimensions::Vector{Int64} ,
#        pcg::Dict{ NTuple{3,NTuple{6,Int64}} , ComplexF64 } , 
#        pcgmat::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,9} } ,
#        cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} } ;
#        verbose=false )
#    # compute the clebsch-gordan series of pseudo-CG coefficients 
#    # of the block states before diagonalization (hence the prima,
#    # and not the tilde). every letter (u,v,i,j,mu,nu) in the function 
#    # refers to the corresponding primed index.
#    # input:
#    # - u_mu_i : (q_u,m_mu,m_i)
#    # - v_nu_j : (q_v,m_nu,m_j)
#    # - q_a  
#    # - oirreps2dimensions : Dict( irrep => dimension )
#    # - pcg : universal pseudo-CG coefficients 
#    #          ( q_nu | c^\dagger_{q_a} | q_mu )
#    # - pcgmat : the same but in the format 
#    #          ( G_nu , G_a , G_mu ) => coefficient matrix
#    # output:
#    # - coefficient
#    #          ( q_u' | c^\dagger | q_v')
#    
#    # type annotations
#    coeff::ComplexF64 = 0.0im
#    cg_muiu::ComplexF64 = 0.0im
#    cg_nujv::ComplexF64 = 0.0im
#    cg_total::ComplexF64 = 0.0im
#
#    # we need all quantum numbers for u,v
#    # and the multiplets for mu,nu,i,j
#    (q_u,m_mu,m_i) = u_mu_i
#    (q_v,m_nu,m_j) = v_nu_j
#
#    # Informative but expensive. Uncomment for testing.
#    #if verbose 
#    #    println( "COMPUTING (N-1) PSG" )
#    #    println( "q_u' = $q_u ; q_v' = $q_v" )
#    #    println( "m_i' = $m_i ; m_j' = $m_j" )
#    #    println( "m_mu' = $m_mu ; m_nu' = $m_nu" ) 
#    #    println( "q_a = $q_a" )
#    #end
#
#    # quantum numbers 
#    # q_u, q_v 
#    (N_u,I_u,S_u,i_u,s_u,r_u) = q_u
#    (N_v,I_v,S_v,i_v,s_v,r_v) = q_v
#    # m_mu, m_nu 
#    (N_mu,I_mu,S_mu,r_mu) = m_mu
#    (N_nu,I_nu,S_nu,r_nu) = m_nu
#    # m_i, m_j
#    (N_i,I_i,S_i,r_i) = (N_ij,I_ij,S_ij,r_ij) = m_i
#    (N_j,I_j,S_j,r_j) = m_j
#
#    # pseudo-cg matrix 
#    G_mu = (N_mu,I_mu,S_mu)
#    G_nu = (N_nu,I_nu,S_nu)
#    G_a  = (N_a,I_a,S_a) = (q_a[1],q_a[2],q_a[3])
#    (i_a,si_a,r_a) = (q_a[4],Int64((q_a[5]+S_a)/2+1),q_a[6]) 
#    ((G_mu,G_a,G_nu) in keys(pcgmat)) || (return zero(ComplexF64))
#    pcg_local = @view pcgmat[G_mu,G_a,G_nu][:,:,:,:,r_mu,r_nu,i_a,si_a,r_a]
#
#    # clebsch-gordan matrices
#    @inbounds begin
#        cgomat_muiu = @view cg_o_fullmatint[I_mu,I_ij,I_u][:,:,i_u]
#        cgomat_nujv = @view cg_o_fullmatint[I_nu,I_ij,I_v][:,:,i_v]
#    end
#
#    # dimension of orbital irreps -> info about partners
#    d_Iij = oindex2dimensions[I_i] # I_i = I_j
#    d_Imu = oindex2dimensions[I_mu]
#    d_Inu = oindex2dimensions[I_nu]
#
#    # main iteration 
#    @inbounds for s_ij in (-S_ij):2:S_ij,
#                  i_ij in 1:d_Iij,
#                  s_nu in (-S_nu):2:S_nu,
#                  i_nu in 1:d_Inu,
#                  s_mu in (-S_mu):2:S_mu,
#                  i_mu in 1:d_Imu
#
#        i_i = i_j = i_ij 
#        s_i = s_j = s_ij        
#
#        if verbose 
#            println( "g_mu' = ( $i_mu , $(s_mu/2.0) )" )
#            println( "g_nu' = ( $i_nu , $(s_nu/2.0) )" )
#        end
#
#        # coefficient
#        si_mu = Int64((s_mu+S_mu)/2+1)
#        si_nu = Int64((s_nu+S_nu)/2+1)
#        pcgcoeff = pcg_local[i_mu,si_mu,i_nu,si_nu]
#        pcgcoeff==zero(pcgcoeff) && continue
#        verbose && println( "pcgcoeff = $pcgcoeff")
#
#        # cg coefficients
#        cg_muiu = cgomat_muiu[i_mu,i_ij] * 
#                PartialWaveFunctions.CG_doublearg(S_mu,s_mu,S_i,s_i,S_u,s_u)
#        cg_nujv = cgomat_nujv[i_nu,i_ij] * 
#                PartialWaveFunctions.CG_doublearg(S_nu,s_nu,S_j,s_j,S_v,s_v) 
#        cg_total = conj( cg_muiu ) * cg_nujv 
#
#        if verbose 
#            println( "cg_muiu = $cg_muiu" )
#            println( "cg_nujv = $cg_nujv" )
#            println( "cg_total = $cg_total" )
#        end
#
#        # putting all together
#        coeff += cg_total * pcgcoeff
#        verbose && println( "coeff = $coeff" )
#    end
#    verbose && println( "final coeff = $coeff" )
#    verbose && println()
#    return coeff 
#end