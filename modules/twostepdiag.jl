# ))))))))))))))))))))))))) # 
#
# TWO-STEP METHOD: OBSOLETE
#
# ))))))))))))))))))))))))) # 


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

