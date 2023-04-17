# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# REDUCED MATRIX REPRESENTATION #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# correction that allows it to work for
# Butterfly, but is intended to be general 
function get_redmat3( dictmat , 
                      multiplets_atom , 
                      multiplets_operator ,
                      cg_o_fullmatint , 
                      cg_s_fullmatint ;
                      verbose=false )

    if verbose 
        println( "~~~~~~~~~~~~~~~~~~~~~~~~~~" )
        println( "REDUCED MATRIX COMPUTATION" )
        println( "~~~~~~~~~~~~~~~~~~~~~~~~~~" )
        println()
        @show multiplets_atom 
        @show multiplets_operator
        println()
    end
    
    # irrep combinations in dictmat
    Gcombs = Set( (q_i[1:3],q_a[1:3],q_j[1:3]) 
                  for (q_i,q_a,q_j) in keys(dictmat) )
    irrmult_ij = get_irreps( Set(multiplets_atom) ; multiplicity=true )
    irrmult_a  = get_irreps( Set(multiplets_operator) ; multiplicity=true )

    # valid irrep combinations
    cg_o_combs = keys(cg_o_fullmatint)
    cg_s_combs = keys(cg_s_fullmatint)
    
    if verbose
        @show irrmult_ij 
        @show irrmult_a
        @show Gcombs
        println()
    end

    # reduced matrix with zeros
    #       < m_i || f^+_{m_a} || m_j > = 0
    irr2mult_ij = Dict( G=>R for (G,R) in irrmult_ij )
    irr2mult_a  = Dict( G=>R for (G,R) in irrmult_a  )
    redmat::Dict{NTuple{3,NTuple{3,Int64}},Array{ComplexF64,3}} = 
        Dict( (G1,G2,G3)=>zeros(ComplexF64,irr2mult_ij[G1],irr2mult_a[G2],irr2mult_ij[G3]) 
              for (G1,G2,G3) in Gcombs )

    # iterate over irreps (with multiplicities)
    @inbounds for (G_i,R_i) in irrmult_ij,
                  (G_a,R_a) in irrmult_a,
                  (G_j,R_j) in irrmult_ij

        ((G_i,G_a,G_j) in Gcombs) || continue

        # irrep quantum numbers
        (N_i,I_i,S_i) = G_i
        (N_j,I_j,S_j) = G_j
        (N_a,I_a,S_a) = G_a
        
        # filter irrep combination
        N_i==(N_j+1)                  || continue 
        ((I_a,I_j,I_i) in cg_o_combs) || continue
        ((S_a,S_j,S_i) in cg_s_combs) || continue

        if verbose
            println( "IRREPS" )
            @show G_i,G_a,G_j
            println()
        end

        # clebsch-gordan matrices
        cgomat_aji = @view cg_o_fullmatint[I_a,I_j,I_i][:,:,:]
        cgsmat_aji = @view cg_s_fullmatint[S_a,S_j,S_i][:,:,:]

        # irrep dimensions
        (dI_a,dI_j,dI_i) = size(cgomat_aji)
        (dS_a,dS_j,dS_i) = size(cgsmat_aji)

        # iterate over irrep multiplets 
        @inbounds for r_i in 1:R_i,
                      r_a in 1:R_a,
                      r_j in 1:R_j 

            m_i = (G_i...,r_i)
            m_a = (G_a...,r_a)
            m_j = (G_j...,r_j)

            verbose && (@show r_i,r_a,r_j)

            # Store  < m_i || f^+_{m_a} || m_j > = 0
            # without having to search in array
            matel = zero(ComplexF64)

            # choose representative

            coinciding = Dict( (q_i,q_a,q_j)=>v
                              for ((q_i,q_a,q_j),v) in dictmat 
                              if ((q_i[1:3]...,q_i[6]),(q_a[1:3]...,q_a[6]),(q_j[1:3]...,q_j[6]))==(m_i,m_a,m_j) )
            length(coinciding)==0 && continue
            representative_key   = collect(keys(coinciding))[1]
            representative_value = coinciding[representative_key]
            q_i,q_a,q_j = representative_key
            if verbose 
                @show representative_key 
                @show dictmat[representative_key] 
            end
            (N_i,I_i,S_i,i_i,s_i,r_i) = q_i
            (N_a,I_a,S_a,i_a,s_a,r_a) = q_a
            (N_j,I_j,S_j,i_j,s_j,r_j) = q_j
            si_i = Int64((s_i+S_i)/2+1)
            si_a = Int64((s_a+S_a)/2+1)
            si_j = Int64((s_j+S_j)/2+1)

            # clebsch-gordan
            cgo = cgomat_aji[i_a,i_j,i_i] 
            cgs = cgsmat_aji[si_a,si_j,si_i]
            cg  = conj(cgo*cgs)

            # compute reduced matrix element
            matel = representative_value/cg
            verbose && (@show matel; println()) 
            # introduce value in final matrix
            redmat[(G_i,G_a,G_j)][r_i,r_a,r_j] = matel
            
        end
    end
    return redmat
end
function get_redmat2( dictmat , 
                      multiplets_atom , 
                      multiplets_operator ,
                      cg_o_fullmatint , 
                      cg_s_fullmatint ;
                      verbose=false )

    if verbose 
        println( "~~~~~~~~~~~~~~~~~~~~~~~~~~" )
        println( "REDUCED MATRIX COMPUTATION" )
        println( "~~~~~~~~~~~~~~~~~~~~~~~~~~" )
        println()
        @show multiplets_atom 
        @show multiplets_operator
        println()
    end
    
    # irrep combinations in dictmat
    Gcombs = Set( (q_i[1:3],q_a[1:3],q_j[1:3]) 
                  for (q_i,q_a,q_j) in keys(dictmat) )

    # [(irrep,multiplicity)]
    irrmult_ij = get_irreps( Set(multiplets_atom) ; multiplicity=true )
    irrmult_a  = get_irreps( Set(multiplets_operator) ; multiplicity=true )

    # valid irrep combinations
    cg_o_combs = keys(cg_o_fullmatint)
    cg_s_combs = keys(cg_s_fullmatint)
    
    if verbose
        @show irrmult_ij 
        @show irrmult_a
        @show Gcombs
        println()
    end

    # reduced matrix with zeros
    #       < m_i || f^+_{m_a} || m_j > = 0
    irr2mult_ij = Dict( G=>R for (G,R) in irrmult_ij )
    irr2mult_a  = Dict( G=>R for (G,R) in irrmult_a  )
    redmat::Dict{NTuple{3,NTuple{3,Int64}},Array{ComplexF64,3}} = 
        Dict( (G1,G2,G3)=>zeros(ComplexF64,irr2mult_ij[G1],irr2mult_a[G2],irr2mult_ij[G3]) 
              for (G1,G2,G3) in Gcombs )

    # iterate over irreps (with multiplicities)
    @inbounds for (G_i,R_i) in irrmult_ij,
                  (G_a,R_a) in irrmult_a,
                  (G_j,R_j) in irrmult_ij

        # irrep quantum numbers
        (N_i,I_i,S_i) = G_i
        (N_j,I_j,S_j) = G_j
        (N_a,I_a,S_a) = G_a
        
        # early discard
        ((G_i,G_a,G_j) in Gcombs)     || continue
        N_i==(N_j+1)                  || continue 
        ((I_a,I_j,I_i) in cg_o_combs) || continue
        ((S_a,S_j,S_i) in cg_s_combs) || continue

        if verbose
            println( "IRREPS" )
            @show G_i,G_a,G_j
            println()
        end

        # clebsch-gordan matrices
        cgomat_aji = @view cg_o_fullmatint[I_a,I_j,I_i][:,:,:]
        cgsmat_aji = @view cg_s_fullmatint[S_a,S_j,S_i][:,:,:]

        # irrep dimensions
        (dI_a,dI_j,dI_i) = size(cgomat_aji)
        (dS_a,dS_j,dS_i) = size(cgsmat_aji)

        # iterate over irrep multiplets 
        @inbounds for r_i in 1:R_i,
                      r_a in 1:R_a,
                      r_j in 1:R_j 

            m_i = (G_i...,r_i)
            m_a = (G_a...,r_a)
            m_j = (G_j...,r_j)

            if verbose 
                println( "MATRIX ELEMENT MULTIPLETS" )
                @show m_i,m_a,m_j
                println()
            end

            # Store  < m_i || f^+_{m_a} || m_j > = 0
            # without having to search in array
            matel = zero(ComplexF64)

            # choose representative
            coinciding = Dict( (q_i,q_a,q_j)=>v
                              for ((q_i,q_a,q_j),v) in dictmat 
                              if ((q_i[1:3]...,q_i[6]),(q_a[1:3]...,q_a[6]),(q_j[1:3]...,q_j[6]))==(m_i,m_a,m_j) )
            length(coinciding)==0 && continue
            representative_key   = collect(keys(coinciding))[1]
            representative_value = coinciding[representative_key]
            q_i,q_a,q_j = representative_key
            (N_i,I_i,S_i,i_i,s_i,r_i) = q_i
            (N_a,I_a,S_a,i_a,s_a,r_a) = q_a
            (N_j,I_j,S_j,i_j,s_j,r_j) = q_j
            si_i = Int64((s_i+S_i)/2+1)
            si_a = Int64((s_a+S_a)/2+1)
            si_j = Int64((s_j+S_j)/2+1)

            if verbose 
                @show q_i,q_a,q_j
                @show representative_value 
                println()
            end

            # clebsch-gordan
            cgo = cgomat_aji[i_a,i_j,i_i] 
            cgs = cgsmat_aji[si_a,si_j,si_i]
            cg  = conj(cgo*cgs)

            # compute reduced matrix element
            matel = representative_value/cg

            # introduce value in final matrix
            redmat[(G_i,G_a,G_j)][r_i,r_a,r_j] = matel
            
        end
    end
    return redmat
end
function get_redmat( dictmat , 
                     multiplets_atom , 
                     multiplets_operator ,
                     cg_o_fullmatint , 
                     cg_s_fullmatint ;
                     verbose=false )

    if verbose 
        println( "~~~~~~~~~~~~~~~~~~~~~~~~~~" )
        println( "REDUCED MATRIX COMPUTATION" )
        println( "~~~~~~~~~~~~~~~~~~~~~~~~~~" )
        println()
        @show multiplets_atom 
        @show multiplets_operator
        println()
    end

    # [(irrep,multiplicity)]
    irrmult_ij = get_irreps( Set(multiplets_atom) ; multiplicity=true )
    irrmult_a  = get_irreps( Set(multiplets_operator) ; multiplicity=true )

    # valid irrep combinations
    cg_o_combs = keys(cg_o_fullmatint)
    cg_s_combs = keys(cg_s_fullmatint)
    
    # irrep combinations in dictmat
    Gcombs = Set( (q_i[1:3],q_a[1:3],q_j[1:3]) 
                  for (q_i,q_a,q_j) in keys(dictmat) )

    if verbose
        @show irrmult_ij 
        @show irrmult_a
        @show Gcombs
        println()
    end

    # reduced matrix with zeros
    #       < m_i || f^+_{m_a} || m_j > = 0
    irr2mult = Dict( G=>R for (G,R) in irrmult_ij )
    merge!( irr2mult , Dict( G=>R for (G,R) in irrmult_a ) )
    redmat::Dict{NTuple{3,NTuple{3,Int64}},Array{ComplexF64,3}} = 
        Dict( (G1,G2,G3)=>zeros(ComplexF64,irr2mult[G1],irr2mult[G2],irr2mult[G3]) 
              for (G1,G2,G3) in Gcombs )

    # iterate over irreps (with multiplicities)
    @inbounds for (G_i,R_i) in irrmult_ij,
                  (G_a,R_a) in irrmult_a,
                  (G_j,R_j) in irrmult_ij

        ((G_i,G_a,G_j) in Gcombs) || continue

        # irrep quantum numbers
        (N_i,I_i,S_i) = G_i
        (N_j,I_j,S_j) = G_j
        (N_a,I_a,S_a) = G_a
        
        # filter irrep combination
        N_i==(N_j+1)                  || continue 
        ((I_a,I_j,I_i) in cg_o_combs) || continue
        ((S_a,S_j,S_i) in cg_s_combs) || continue

        if verbose
            println( "IRREPS" )
            @show G_i,G_a,G_j
            println()
        end
        #verbose=false 
        #if G_i==(3,5,1) && G_j==(2,2,0) 
        #    verbose=true
        #end

        # clebsch-gordan matrices
        cgomat_aji = @view cg_o_fullmatint[I_a,I_j,I_i][:,:,:]
        cgsmat_aji = @view cg_s_fullmatint[S_a,S_j,S_i][:,:,:]

        # irrep dimensions
        (dI_a,dI_j,dI_i) = size(cgomat_aji)
        (dS_a,dS_j,dS_i) = size(cgsmat_aji)

        # iterate over irrep multiplets 
        @inbounds for r_i in 1:R_i,
                      r_a in 1:R_a,
                      r_j in 1:R_j 

            if verbose 
                println( "MATRIX ELEMENT MULTIPLETS" )
                m_i = (G_i...,r_i)
                m_a = (G_a...,r_a)
                m_j = (G_j...,r_j)
                @show m_i,m_a,m_j
                println()
            end

            # Store  < m_i || f^+_{m_a} || m_j > = 0
            # without having to search in array
            matel = zero(ComplexF64)

            # iterate over g_j,g_a,g_i,g'_i
            @inbounds for i_i  in 1:dI_i,
                          i_a  in 1:dI_a,
                          i_j  in 1:dI_j,
                          i_ip in 1:dI_i,
                          (si_i,s_i)   in enumerate((-S_i):2:S_i),
                          (si_a,s_a)   in enumerate((-S_a):2:S_a),
                          (si_j,s_j)   in enumerate((-S_j):2:S_j),
                          (si_ip,s_ip) in enumerate((-S_i):2:S_i)

                # full quantum numbers 
                q_i = (G_i...,i_i,s_i,r_i)
                q_a = (G_a...,i_a,s_a,r_a)
                q_j = (G_j...,i_j,s_j,r_j)

                # clebsch-gordan coefficient
                cg_o = cgomat_aji[i_a, i_j, i_ip]
                cg_s = cgsmat_aji[si_a,si_j,si_ip]
                cg   = cg_o*cg_s

                # dictmat coefficient
                dm   = get(dictmat,(q_i,q_a,q_j),zero(ComplexF64))

                # add to this matrix element
                matel += dm * cg

                if verbose
                    println( "matrix element contribution" )
                    @show q_i,q_a,q_j
                    @show cg, cg_o, cg_s 
                    @show dm
                    @show dm*cg
                    @show matel
                    println()
                end
            end
            
            # include dimension factor 
            matel /= (dI_i*dS_i)
            if verbose
                @show matel 
                println( "-----------------------------" )
            end

            # introduce value in final matrix
            redmat[(G_i,G_a,G_j)][r_i,r_a,r_j] = matel
            
        end
    end
    return redmat
end


function compute_pcgred_iaj_full(
                multiplets_block::Set{NTuple{4,Int64}}, 
                multiplets_shell::Set{NTuple{4,Int64}},
                multiplets_a_block::Vector{NTuple{4,Int64}}, 
                Csum_o_array::Array{ComplexF64,6} ,
                Csum_s_array::Array{ComplexF64,6} ,
                pcgred_block::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} },
                pcgred_shell::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} },
                cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
                cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
                combs_uvprima::Dict{ NTuple{2,NTuple{3,Int64}} , Vector{NTuple{6,NTuple{4,Int64}}} },
                combinations_uprima_new::Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} },
                irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} })

    println( "FULL IAJ MATRIX CALCULATION" )


    G2R_uv::Dict{NTuple{3,Int64},Int64} = Dict()
    @inbounds for (G_uv,mults) in combinations_uprima_new,
        (m_u,m_mu,m_i) in mults

        G_i,r_i = m_i[1:3],m_i[4]

        if G_i in keys(G2R_uv)
            if r_i>G2R_uv[G_i]
                G2R_uv[G_i] = r_i
            end
        else
            G2R_uv[G_i] = r_i
        end
    end
    G2R_a::Dict{NTuple{3,Int64},Int64} = 
        Dict( G=>R for (G,R) in get_irreps( Set(multiplets_a_block) ; multiplicity=true ) )

    pcgred_iaj_full::Dict{NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } = 
            Dict( (G_u,G_a,G_v)=>zeros(ComplexF64,R_u,R_a,R_v)
                            for (G_u,R_u) in G2R_uv for (G_v,R_v) in G2R_uv
                            for (G_a,R_a) in G2R_a )

    mm_i::Set{NTuple{4,Int64}} = Set( m_i for (G_uv,combs) in combinations_uprima_new 
                                          for (m_u,m_mu,m_i) in combs )
    GG_i::Set{NTuple{3,Int64}} = Set( m_i[1:3] for m_i in mm_i )
    G2rr_i::Dict{NTuple{3,Int64} , Vector{Int64}} = Dict( G_i=>sort([m_i[4] for m_i in mm_i if m_i[1:3]==G_i])
                   for G_i in GG_i ) 

   @inbounds for (G_i::NTuple{3,Int64},rr_i::Vector{Int64}) in G2rr_i,
                 (G_j::NTuple{3,Int64},rr_j::Vector{Int64}) in G2rr_i,
                 (G_a::NTuple{3,Int64},R_a::Int64) in G2R_a

        N_i,I_i,S_i    = G_i
        N_j,I_j,S_j    = G_j
        N_a,I_a,S_a    = G_a

        N_i==(N_j+1) || continue
        ((I_a,I_j,I_i) in keys(cg_o_fullmatint)) || continue
        ((S_a,S_j,S_i) in keys(cg_s_fullmatint)) || continue

        U_i::Matrix{ComplexF64} = irrEU[G_i][2]
        U_j::Matrix{ComplexF64} = irrEU[G_j][2]

        combs_uvprima_local::Vector{NTuple{6,NTuple{4,Int64}}} = combs_uvprima[(G_i,G_j)]

        zeromat::Array{ComplexF64,3} = zeros(ComplexF64,G2R_uv[G_i],G2R_a[G_a],G2R_uv[G_j])
        pcgred_gcomb::Array{ComplexF64,3} = zeromat

        @inbounds for r_j::Int64 in rr_j,
                      r_a::Int64 in 1:R_a,
                      r_i::Int64 in rr_i

            m_i::NTuple{4,Int64} = (G_i...,r_i)
            m_j::NTuple{4,Int64} = (G_j...,r_j)
            m_a::NTuple{4,Int64} = (G_a...,r_a)

            pcgred_gcomb[r_i,r_a,r_j] = 
                        compute_pcgred_iaj(
                            m_i , m_a , m_j ,
                            Csum_o_array ,
                            Csum_s_array ,
                            pcgred_block ,
                            combs_uvprima_local ,
                            U_i ,# = U_ut 
                            U_j ,# = U_vt
                            zeromat )

        end

        if !isapprox(sum(abs2.(pcgred_gcomb)),0.0)
            pcgred_iaj_full[(G_i,G_a,G_j)] = pcgred_gcomb
        end

    end

    return pcgred_iaj_full
end

# reduced matrix method
function matdiag_redmat( 
        multiplets_block::Set{NTuple{4,Int64}}, 
        multiplets_shell::Set{NTuple{4,Int64}},
        irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
        hop_symparams::Dict{ Int64 , Matrix{ComplexF64} },
        cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        Csum_o_array::Array{ComplexF64,6} ,
        Csum_s_array::Array{ComplexF64,6} ,
        Bsum_o_array::Array{ComplexF64,6} ,
        Bsum_s_array::Array{ComplexF64,6} ,
        pcgred_block::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} },
        pcgred_shell::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} },
        multiplets_a_block::Vector{NTuple{4,Int64}}, 
        multiplets_a_shell::Vector{NTuple{4,Int64}}, 
        combinations_uprima::Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} },
        oindex2dimensions::Vector{ Int64 };
        distributed=false,
        impinfo=false ,
        verbose=false ,
        precompute_iaj=true)
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
    #

    verbose && println( "CALCULATION OF <u||H||v> MATRIX\n" )

    # block-shell combination multiplets 
    multiplets_mui2u::Dict{NTuple{2,NTuple{4,Int64}},Vector{NTuple{4,Int64}}}  = 
            get_combination_multiplets( multiplets_block , 
                                        multiplets_shell , 
                                        Set(keys(cg_o_fullmatint)) ,
                                        Set(keys(cg_s_fullmatint)) ;
                                        verbose=false )
    
    # new block-shell combinations 
    GG_u = Set( m_u[1:3] for mm_u in values(multiplets_mui2u) 
                         for m_u in mm_u )
    combinations_uprima_new::Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} } = Dict( 
                G => NTuple{3,NTuple{4,Int64}}[ (m_u,m_mu,m_i) 
                                                for ((m_mu,m_i),mm_u) in multiplets_mui2u 
                                                for m_u in mm_u
                                                if m_u[1:3]==G 
                                              ] 
                for G in GG_u
    )
    combinations_uprima_new_vec = collect( (x[1],x[2]) for x in combinations_uprima_new )
    
    irrEU_new::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} } = Dict()

    combs_uvprima::Dict{ NTuple{2,NTuple{3,Int64}} , Vector{NTuple{6,NTuple{4,Int64}}} } = Dict()
    @inbounds for (G_u,mm_u) in combinations_uprima,
                  (G_v,mm_v) in combinations_uprima

        combs_uvprima[(G_u,G_v)] = [(m_u,m_mu,m_i,m_v,m_nu,m_j) 
                                   for (m_u,m_mu,m_i) in mm_u
                                   for (m_v,m_nu,m_j) in mm_v 
                                   if m_i==m_j]

    end

    pcgred_iaj_full::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } = Dict()
    if precompute_iaj
        pcgred_iaj_full= compute_pcgred_iaj_full( multiplets_block, 
                                                  multiplets_shell,
                                                  multiplets_a_block, 
                                                  Csum_o_array,
                                                  Csum_s_array,
                                                  pcgred_block,
                                                  pcgred_shell,
                                                  cg_o_fullmatint,
                                                  cg_s_fullmatint,
                                                  combs_uvprima,
                                                  combinations_uprima_new,
                                                  irrEU)
    end

    if distributed 
        @everywhere multiplets_mui2u    = $multiplets_mui2u
        @everywhere oindex2dimensions   = $oindex2dimensions::Vector{ Int64 }
        @everywhere hop_symparams       = $hop::Dict{ Int64 , Matrix{ComplexF64} }
        @everywhere cg_o_fullmatint     = $cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} }
        @everywhere cg_s_fullmatint     = $cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} }
        @everywhere Csum_o_array        = $Csum_o_array::Array{ComplexF64,6} 
        @everywhere Csum_s_array        = $Csum_s_array::Array{ComplexF64,6} 
        @everywhere Bsum_o_array        = $Bsum_o_array::Array{ComplexF64,6} 
        @everywhere Bsum_s_array        = $Bsum_s_array::Array{ComplexF64,6} 
        @everywhere pcg_block           = $pcgred_block::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} }
        @everywhere pcg_shell           = $pcgred_shell::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} }
        @everywhere multiplets_a_block  = $multiplets_a_block::Vector{NTuple{4,Int64}} 
        @everywhere multiplets_a_shell  = $multiplets_a_shell::Vector{NTuple{4,Int64}} 
        @everywhere combinations_uprima = $combinations_uprima::Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} }
        @everywhere irrEU               = $irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} }
        @everywhere combs_uvprima       = $combs_uvprima::Dict{ NTuple{2,NTuple{3,Int64}} , Vector{NTuple{6,NTuple{4,Int64}}} }
    end
    
    if distributed
        irrEU_new = @distributed (merge) for 
                        (G_uv,mults) in combinations_uprima_new_vec

            construct_diag_block(
                        G_uv ,
                        mults ,
                        irrEU ,
                        combs_uvprima ,
                        cg_o_fullmatint ,
                        cg_s_fullmatint ,
                        Csum_o_array ,
                        Csum_s_array ,
                        Bsum_o_array ,
                        Bsum_s_array ,
                        pcgred_block , 
                        pcgred_shell ,
                        multiplets_a_block ,
                        multiplets_a_shell ,
                        hop_symparams ,
                        oindex2dimensions ;
                        verbose=verbose )

        end
    else
        for (G_uv,mults) in combinations_uprima_new_vec 

            merge!( irrEU_new ,
                    construct_diag_block(
                        G_uv ,
                        mults ,
                        irrEU ,
                        combs_uvprima ,
                        cg_o_fullmatint ,
                        cg_s_fullmatint ,
                        Csum_o_array ,
                        Csum_s_array ,
                        Bsum_o_array ,
                        Bsum_s_array ,
                        pcgred_block , 
                        pcgred_shell ,
                        multiplets_a_block ,
                        multiplets_a_shell ,
                        hop_symparams ,
                        oindex2dimensions ;
                        verbose=verbose,
                        precompute_iaj=precompute_iaj,
                        pcgred_iaj_full=pcgred_iaj_full) )

        end
    end

    minE = minimum([e for (E,U) in values(irrEU_new) for e in E])
    irrEU_new = Dict( G=>(E.-minE,U) for (G,(E,U)) in irrEU_new )

    return ( irrEU_new , combinations_uprima_new )
end

function construct_diag_block(
            G_uv::NTuple{3,Int64} ,
            mults::Vector{NTuple{3,NTuple{4,Int64}}} ,
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            combs_uvprima::Dict{ NTuple{2,NTuple{3,Int64}} , Vector{NTuple{6,NTuple{4,Int64}}} },
            cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            Csum_o_array::Array{ComplexF64,6} ,
            Csum_s_array::Array{ComplexF64,6} ,
            Bsum_o_array::Array{ComplexF64,6} ,
            Bsum_s_array::Array{ComplexF64,6} ,
            pcgred_block::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} },
            pcgred_shell::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} },
            multiplets_a_block::Vector{NTuple{4,Int64}} ,
            multiplets_a_shell::Vector{NTuple{4,Int64}} ,
            hop_symparams::Dict{ Int64 , Matrix{ComplexF64} },
            oindex2dimensions::Vector{Int64} ;
            verbose=false ,
            precompute_iaj::Bool=false,
            pcgred_iaj_full::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } = Dict() )

        verbose && println( "COMPUTING AND DIAGONALIZING G_uv=$G_uv BLOCK...\n" )

        zeromat::Array{ComplexF64,3} = zeros(ComplexF64,20,length(multiplets_a_block),20)

        hblock::Array{ComplexF64,3} = zeros( ComplexF64 , length(mults) , length(mults) , 2 )

        @inbounds for (m_u::NTuple{4,Int64},m_mu::NTuple{4,Int64},m_i::NTuple{4,Int64}) in mults,
                      (m_v::NTuple{4,Int64},m_nu::NTuple{4,Int64},m_j::NTuple{4,Int64}) in mults

            if verbose 
                println( "Hamiltonian block" )
                @show m_u, m_mu, m_i 
                @show m_v, m_nu, m_j
                println()
            end

            # block and shell quantum numbers 
            (N_i, I_i, S_i, r_i ) = m_i
            (N_j, I_j, S_j, r_j ) = m_j
            G_i = (N_i, I_i, S_i)
            G_j = (N_j, I_j, S_j)
            (N_mu,I_mu,S_mu,r_mu) = m_mu
            (N_nu,I_nu,S_nu,r_nu) = m_nu
            r_u = m_u[4]
            r_v = m_v[4]

            # diagonal energy from previous step
            #   δ_{m_u,m_v} E^[n-1]_(m_i)
            # (m_u=m_v => m_i=m_j)
            hblock[r_u,r_v,1] = m_u==m_v ? irrEU[G_i][1][r_i] : 
                                           zero(ComplexF64)
            verbose && println( "δ_{m_u,m_v} E^[n-1]_(m_i) = $(hblock[r_u,r_v,1])\n" )

            # chosen partner: orbital 1, max z-spin
            p = (N,I,S,i,s) = ( G_uv... , 1 , G_uv[3] )
            verbose && println( "p_u = p_v = $((N,I,S,i,s))" )

            # U matrices (not transposed for faster loop)
            U_i::Matrix{ComplexF64} = irrEU[G_i][2]
            U_j::Matrix{ComplexF64} = irrEU[G_j][2]

            # hopping part
            verbose && println( "Occupations: (N_nu=$N_nu,N_mu=$N_mu) and (N_i=$N_i,N_j=$N_j)\n" )
            if (N_nu==(N_mu+1) && N_i==(N_j+1)) 

                verbose && println( "Hopping possible\n" )
                
                # (mu',i'=>u'), (nu',j'=>v') combinations 
                combs_uvprima_local = combs_uvprima[(G_i,G_j)]

                # hopping element of the hamiltonian matrix
                hblock[r_u,r_v,2] = compute_hopelement_redmat( 
                    m_i  , m_j , 
                    m_mu , m_nu ,
                    G_uv , 
                    cg_o_fullmatint ,
                    cg_s_fullmatint ,
                    Csum_o_array ,
                    Csum_s_array ,
                    Bsum_o_array ,
                    Bsum_s_array ,
                    pcgred_block , 
                    pcgred_shell ,
                    multiplets_a_block ,
                    multiplets_a_shell ,
                    hop_symparams ,
                    oindex2dimensions ,
                    U_i , U_j ,
                    combs_uvprima_local ,
                    zeromat ;
                    verbose=verbose ,
                    precompute_iaj=precompute_iaj,
                    pcgred_iaj_full=pcgred_iaj_full)
            end
        end
        
        # [E,hop] --> realE
        R::Int64 = length(mults)
        hblock_complete::Matrix{ComplexF64} = zeros(ComplexF64,R,R)
        @inbounds for r_u::Int64 in 1:R, 
                      r_v::Int64 in 1:R 
            hblock_complete[r_u,r_v] = sum(hblock[r_u,r_v,:]) + conj(hblock[r_v,r_u,2])
        end

        # diagonalize
        hblock_complete = 0.5*( hblock_complete + hblock_complete' )
        F = eigen( hblock_complete )
        (e,u) = ( real(F.values) , F.vectors )
        e = real.(e)

        if verbose 
            println( "RESULTS FOR G_uv=$G_uv BLOCK" )
            for m in mults 
                @show m 
            end
            @show e 
            @show u
            println()
        end

        return Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} }( G_uv => (e,u) )

end

# reduced pcgmat method
function compute_hopelement_redmat( 
        m_i::NTuple{4,Int64}  , m_j::NTuple{4,Int64} , 
        m_mu::NTuple{4,Int64} , m_nu::NTuple{4,Int64} , 
        G_uv::NTuple{3,Int64} ,
        cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        Csum_o_array::Array{ComplexF64,6} ,
        Csum_s_array::Array{ComplexF64,6} ,
        Bsum_o_array::Array{ComplexF64,6} ,
        Bsum_s_array::Array{ComplexF64,6} ,
        pcgred_block::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} },
        pcgred_shell::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} },
        multiplets_a_block::Vector{NTuple{4,Int64}}, 
        multiplets_a_shell::Vector{NTuple{4,Int64}}, 
        hop_symparams::Dict{ Int64 , Matrix{ComplexF64} } ,
        oindex2dimensions::Vector{ Int64 } ,
        U_i::Matrix{ComplexF64}, 
        U_j::Matrix{ComplexF64},
        combs_uvprima_local::Vector{NTuple{6,NTuple{4,Int64}}},
        zeromat::Array{ComplexF64,3} ;
        verbose::Bool=false ,
        precompute_iaj::Bool=false,
        pcgred_iaj_full::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } = Dict() )


    if verbose 
        println( "Computing hopping element for multiplets" )
        @show G_uv 
        @show m_i,m_mu 
        @show m_j,m_nu 
        println()
    end

    multiplets_a_combs::Vector{NTuple{2,NTuple{4,Int64}}} = [
        (m_a_block,m_a_shell) for m_a_block in multiplets_a_block 
                              for m_a_shell in multiplets_a_shell 
                              if m_a_block[2]==m_a_shell[2] 
    ]

    hopelement::ComplexF64 = zero(ComplexF64)

    # |u) and |v) irrep and partners 
    (N,I,S,i,s) = ( G_uv... , 1 , G_uv[3] )
    si = Int64((s+S)/2+1)

    # quantum numbers
    (N_i,I_i,S_i,r_i) = m_i
    (N_j,I_j,S_j,r_j) = m_j
    (N_mu,I_mu,S_mu,r_mu) = m_mu
    (N_nu,I_nu,S_nu,r_nu) = m_nu
    (N_uv,I_uv,S_uv) = G_uv
    G_i = (N_i,I_i,S_i)
    G_j = (N_j,I_j,S_j)
    G_mu = m_mu[1:3] 
    G_nu = m_nu[1:3]

    # sign factor 
    sign = (-1)^N_mu
    verbose && println( "Sign factor: $sign\n" )

    # hopping multiplet iteration 
    verbose && println( "Hopping multiplet iteration..." )
    @inbounds for (m_a_block,m_a_shell) in multiplets_a_combs

        verbose && (@show m_a_block; @show m_a_shell)

        # hopper quantum numbers 
        G_a = m_a_block[1:3] 
        (N_a,I_a,S_a) = G_a
        r_a_block = m_a_block[4]
        r_a_shell = m_a_shell[4]

        # early discard
        ((G_nu,G_a,G_mu) in keys(pcgred_shell))  || continue
        ((I_a,I_j,I_i) in keys(cg_o_fullmatint)) || continue
        ((S_a,S_j,S_i) in keys(cg_s_fullmatint)) || continue

        verbose && println( "symmetry discard passed\n" )

        # channel hopping parameter
        hoparam = hop_symparams[m_a_block[2]][r_a_block,r_a_shell]
        verbose && println( "hoparam = $hoparam\n" )
        hoparam==zero(hoparam) && continue

        # shell pcg 
        pcgred_nuamu = pcgred_shell[(G_nu,G_a,G_mu)][r_nu,r_a_shell,r_mu]
        verbose && println( "( m_nu || f^†_{m_a} || m_mu ) = $pcgred_nuamu\n" )
        pcgred_nuamu==zero(ComplexF64) && continue

        # clebsch gordan sum (B) 
        sidx = (S_mu,S_nu,S_i,S_j,S_uv,S_a).+1
        B = Bsum_o_array[I_mu,I_nu,I_i,I_j,I_uv,I_a]*
            Bsum_s_array[sidx...]
        verbose && println( "B sum = $B\n" )
        B==zero(B) && continue

        if precompute_iaj
            if ((G_i,G_a,G_j) in keys(pcgred_iaj_full))
                pcgred_iaj = pcgred_iaj_full[G_i,G_a,G_j][r_i,r_a_block,r_j]
            else
                pcgred_iaj = zero(pcgred_iaj)
            end
        else
            pcgred_iaj = compute_pcgred_iaj(
                            m_i , m_a_block , m_j ,
                            Csum_o_array ,
                            Csum_s_array ,
                            pcgred_block ,
                            combs_uvprima_local ,
                            U_i ,# = U_ut 
                            U_j ,# = U_vt
                            zeromat ;
                            verbose=verbose )
        end
        verbose && println( "( m_i || f^†_{m_a} || m_j ) = $pcgred_iaj\n" )

        hopelement += sign*
                      hoparam*
                      B*
                      pcgred_iaj*
                      conj(pcgred_nuamu)
        if verbose 
            c = sign*hoparam*B*pcgred_iaj*conj(pcgred_nuamu)
            println( "contribution to hopelement = $c\n" )
        end

    end

    verbose && println( "hopelement = $hopelement\n" )

    return hopelement 

end

#
# < m_i || f^†[n-1]_{m_a} || m_j >
#
function compute_pcgred_iaj(
            m_i::NTuple{4,Int64} ,
            m_a::NTuple{4,Int64} ,
            m_j::NTuple{4,Int64} ,
            Csum_o_array::Array{ComplexF64,6} ,
            Csum_s_array::Array{ComplexF64,6} ,
            pcgred_block::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} },
            combs_uvprima_local::Vector{NTuple{6,NTuple{4,Int64}}},
            U_ut::Matrix{ComplexF64} , 
            U_vt::Matrix{ComplexF64} ,
            zeromat::Array{ComplexF64,3} ;
            verbose=false )

    if verbose 
        println( "Computing ( i || f^†_a || j )" )
        @show m_i, m_a, m_j
        println()
    end

    pcgred_iaj::ComplexF64 = zero(ComplexF64)

    # redmat key quantum numbers
    G_ut::NTuple{3,Int64},r_ut::Int64 = m_i[1:3],m_i[4]
    G_vt::NTuple{3,Int64},r_vt::Int64 = m_j[1:3],m_j[4]
    G_a::NTuple{3,Int64}, r_a::Int64  = m_a[1:3],m_a[4]
    (N_a,I_a,S_a) = G_a 
    
    redmatel::ComplexF64 = zero(ComplexF64)

    # iterate over multiplets in G_ut=G_i and G_vt=G_j
    @inbounds for (m_utp::NTuple{4,Int64},m_mut::NTuple{4,Int64},m_it::NTuple{4,Int64},m_vtp::NTuple{4,Int64},m_nut::NTuple{4,Int64},m_jt::NTuple{4,Int64}) in combs_uvprima_local

        # early discard of unequal block parts
        #m_it==m_jt || continue
        m_ijt = m_it

        # multiplet quantum numbers 
        (N_ut,I_ut,S_ut,r_utp)    = m_utp 
        (N_vt,I_vt,S_vt,r_vtp)    = m_vtp
        (N_ijt,I_ijt,S_ijt,r_ijt) = m_ijt
        (N_mut,I_mut,S_mut,r_mut) = m_mut
        (N_nut,I_nut,S_nut,r_nut) = m_nut
        G_utp = m_utp[1:3] 
        G_vtp = m_vtp[1:3]
        G_ijt = m_ijt[1:3]
        (N_utp,I_utp,S_utp) = G_utp
        (N_vtp,I_vtp,S_vtp) = G_vtp
        (N_ijt,I_ijt,S_ijt) = G_ijt
        G_mut::NTuple{3,Int64},r_mut = m_mut[1:3],m_mut[4]
        G_nut::NTuple{3,Int64},r_nut = m_nut[1:3],m_nut[4]
        
        # U matrix cofficient
        uu::ComplexF64 = (conj(U_ut[r_utp,r_ut])*U_vt[r_vtp,r_vt])::ComplexF64
        #isapprox( uu , zero(uu) ) && continue
        
        # reduced matrix element
        #haskey( pcgred_block , (G_mut,G_a,G_nut) ) || continue
        #redmatel = pcgred_block[(G_mut,G_a,G_nut)][r_mut,r_a,r_nut]
        redmatel = get( pcgred_block , (G_mut,G_a,G_nut) , zeromat )[r_mut,r_a,r_nut]::ComplexF64

        # clebsch-gordan sum
        #C = get( Csum_o_dict , (I_utp,I_vtp,I_ijt,I_mut,I_nut,I_a) , zero(ComplexF64) )*
        #    get( Csum_s_dict , (S_utp,S_vtp,S_ijt,S_mut,S_nut,S_a) , zero(ComplexF64) )
        sidx = (S_utp,S_vtp,S_ijt,S_mut,S_nut,S_a).+1
        #os = Csum_o_array[I_utp,I_vtp,I_ijt,I_mut,I_nut,I_a]
        #ss = Csum_s_array[sidx...]
        C = Csum_o_array[I_utp,I_vtp,I_ijt,I_mut,I_nut,I_a]*
            Csum_s_array[sidx...]

        # put all together
        pcgred_iaj += uu * C * redmatel

    end

    verbose && (@show pcgred_iaj)

    return pcgred_iaj::ComplexF64

end


# CG sum related to the decompositions 
# u->mu,i  and  v->nu,j 
function compute_CG_Bsum(
        oindex2dimensions::Vector{Int64} ,
        cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        IS_mu::NTuple{2,Int64},
        IS_nu::NTuple{2,Int64},
        IS_i::NTuple{2,Int64},
        IS_j::NTuple{2,Int64},
        IS_uv::NTuple{2,Int64},
        IS_a::NTuple{2,Int64} )

    # irrep quantum numbers
    (I_mu,S_mu) = IS_mu
    (I_nu,S_nu) = IS_nu
    (I_i,S_i)   = IS_i
    (I_j,S_j)   = IS_j
    (I_uv,S_uv) = IS_uv
    (I_a,S_a)   = IS_a

    # selected partner 
    i_uv = 1 
    si_uv = S_uv+1

    # check irrep combination possibility
    ((I_a,I_j,I_i)   in keys(cg_o_fullmatint)) || (return zero(ComplexF64))
    ((S_a,S_j,S_i)   in keys(cg_s_fullmatint)) || (return zero(ComplexF64))
    ((I_a,I_mu,I_nu) in keys(cg_o_fullmatint)) || (return zero(ComplexF64))
    ((S_a,S_mu,S_nu) in keys(cg_s_fullmatint)) || (return zero(ComplexF64))

    # local matrices 
    @inbounds begin 

        cgomat_muiu = @view cg_o_fullmatint[(I_mu,I_i,I_uv)][:,:,i_uv]
        cgsmat_muiu = @view cg_s_fullmatint[(S_mu,S_i,S_uv)][:,:,si_uv]

        cgomat_nujv = @view cg_o_fullmatint[(I_nu,I_j,I_uv)][:,:,i_uv] 
        cgsmat_nujv = @view cg_s_fullmatint[(S_nu,S_j,S_uv)][:,:,si_uv]

        cgomat_aji = @view cg_o_fullmatint[(I_a,I_j,I_i)][:,:,:]
        cgsmat_aji = @view cg_s_fullmatint[(S_a,S_j,S_i)][:,:,:]

        cgomat_amunu = @view cg_o_fullmatint[(I_a,I_mu,I_nu)][:,:,:]
        cgsmat_amunu = @view cg_s_fullmatint[(S_a,S_mu,S_nu)][:,:,:]

    end

    # irrep dimensions 
    dI_mu = oindex2dimensions[I_mu]
    dI_nu = oindex2dimensions[I_nu]
    dI_i  = oindex2dimensions[I_i]
    dI_j  = oindex2dimensions[I_j]
    dI_a  = oindex2dimensions[I_a]
    dS_mu = S_mu + 1
    dS_nu = S_nu + 1
    dS_i  = S_i + 1
    dS_j  = S_j + 1
    dS_a  = S_a + 1

    # orbital sum 
    sum_o::ComplexF64 = zero(ComplexF64)
    @inbounds for i_i  in 1:dI_i,
                  i_j  in 1:dI_j,
                  i_nu in 1:dI_nu,
                  i_mu in 1:dI_mu,
                  i_a  in 1:dI_a

        sum_o += conj(cgomat_muiu[i_mu,i_i])*
                 cgomat_nujv[i_nu,i_j]*
                 conj(cgomat_amunu[i_a,i_mu,i_nu])*
                 cgomat_aji[i_a,i_j,i_i]

    end
    # spin sum 
    sum_s::ComplexF64 = zero(ComplexF64)
    @inbounds for si_i  in 1:dS_i,
                  si_j  in 1:dS_j,
                  si_nu in 1:dS_nu,
                  si_mu in 1:dS_mu,
                  si_a  in 1:dS_a

        sum_s += conj(cgsmat_muiu[si_mu,si_i])*
                 cgsmat_nujv[si_nu,si_j]*
                 conj(cgsmat_amunu[si_a,si_mu,si_nu])*
                 cgsmat_aji[si_a,si_j,si_i]

    end

    return (sum_o*sum_s)::ComplexF64

end
function compute_CG_Bsum_orbital(
        oindex2dimensions::Vector{Int64} ,
        cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        I_mu::Int64,
        I_nu::Int64,
        I_i::Int64,
        I_j::Int64,
        I_uv::Int64,
        I_a::Int64 )

    # selected partner 
    i_uv = 1 

    # check irrep combination possibility
    ((I_a,I_j,I_i)   in keys(cg_o_fullmatint)) || (return zero(ComplexF64))
    ((I_a,I_mu,I_nu) in keys(cg_o_fullmatint)) || (return zero(ComplexF64))

    # local matrices 
    @inbounds begin 
        cgomat_muiu  = @view cg_o_fullmatint[(I_mu,I_i,I_uv)][:,:,i_uv]
        cgomat_nujv  = @view cg_o_fullmatint[(I_nu,I_j,I_uv)][:,:,i_uv] 
        cgomat_aji   = @view cg_o_fullmatint[(I_a,I_j,I_i)][:,:,:]
        cgomat_amunu = @view cg_o_fullmatint[(I_a,I_mu,I_nu)][:,:,:]
    end

    # irrep dimensions 
    dI_mu = oindex2dimensions[I_mu]
    dI_nu = oindex2dimensions[I_nu]
    dI_i  = oindex2dimensions[I_i]
    dI_j  = oindex2dimensions[I_j]
    dI_a  = oindex2dimensions[I_a]

    # orbital sum 
    sum_o::ComplexF64 = zero(ComplexF64)
    @inbounds for i_i  in 1:dI_i,
                  i_j  in 1:dI_j,
                  i_nu in 1:dI_nu,
                  i_mu in 1:dI_mu,
                  i_a  in 1:dI_a

        sum_o += conj(cgomat_muiu[i_mu,i_i])*
                 cgomat_nujv[i_nu,i_j]*
                 cgomat_amunu[i_a,i_mu,i_nu]*
                 conj(cgomat_aji[i_a,i_j,i_i])

    end

    return sum_o::ComplexF64

end
function compute_CG_Bsum_spin(
        cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        S_mu::Int64, 
        S_nu::Int64,
        S_i::Int64,
        S_j::Int64,
        S_uv::Int64,
        S_a::Int64)

    # selected partner 
    si_uv = S_uv+1

    # check irrep combination possibility
    ((S_a,S_j,S_i)   in keys(cg_s_fullmatint)) || (return zero(ComplexF64))
    ((S_a,S_mu,S_nu) in keys(cg_s_fullmatint)) || (return zero(ComplexF64))

    # local matrices 
    @inbounds begin 
        cgsmat_muiu  = @view cg_s_fullmatint[(S_mu,S_i,S_uv)][:,:,si_uv]
        cgsmat_nujv  = @view cg_s_fullmatint[(S_nu,S_j,S_uv)][:,:,si_uv]
        cgsmat_aji   = @view cg_s_fullmatint[(S_a,S_j,S_i)][:,:,:]
        cgsmat_amunu = @view cg_s_fullmatint[(S_a,S_mu,S_nu)][:,:,:]
    end

    # irrep dimensions 
    dS_mu = S_mu + 1
    dS_nu = S_nu + 1
    dS_i  = S_i + 1
    dS_j  = S_j + 1
    dS_a  = S_a + 1

    # spin sum 
    sum_s::ComplexF64 = zero(ComplexF64)
    @inbounds for si_i  in 1:dS_i,
                  si_j  in 1:dS_j,
                  si_nu in 1:dS_nu,
                  si_mu in 1:dS_mu,
                  si_a  in 1:dS_a

        sum_s += conj(cgsmat_muiu[si_mu,si_i])*
                 cgsmat_nujv[si_nu,si_j]*
                 cgsmat_amunu[si_a,si_mu,si_nu]*
                 conj(cgsmat_aji[si_a,si_j,si_i])

    end

    return sum_s::ComplexF64

end

#
# CG sum related to the decomposition 
#
#    utp->mut,it and vtp->nut,jt,
#
# where we have the condition that i=j
#
function compute_CG_Csum( 
            oindex2dimensions::Vector{Int64} ,
            cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            IS_utp::NTuple{2,Int64} ,
            IS_vtp::NTuple{2,Int64} ,
            IS_ijt::NTuple{2,Int64} ,
            IS_mut::NTuple{2,Int64} ,
            IS_nut::NTuple{2,Int64} ,
            IS_a::NTuple{2,Int64} )
    
    # * we drop the "t" (tilde) and "p" (prima)
    # * for an easier reading

    # irrep quantum numbers
    (I_u,S_u)   = IS_utp
    (I_v,S_v)   = IS_vtp
    (I_ij,S_ij) = IS_ijt
    (I_mu,S_mu) = IS_mut
    (I_nu,S_nu) = IS_nut
    (I_a,S_a)   = IS_a

    # irrep dimensions 
    dI_u::Int64  = oindex2dimensions[I_u] 
    dI_v::Int64  = oindex2dimensions[I_v]
    dI_ij::Int64 = oindex2dimensions[I_ij] 
    dI_mu::Int64 = oindex2dimensions[I_mu] 
    dI_nu::Int64 = oindex2dimensions[I_nu]
    dI_a::Int64  = oindex2dimensions[I_a]
    dS_u::Int64  = S_u+1
    dS_v::Int64  = S_v+1
    dS_ij::Int64 = S_ij+1 
    dS_mu::Int64 = S_mu+1
    dS_nu::Int64 = S_nu+1
    dS_a::Int64  = S_a+1

    # clebsch-gordan matrices
    @inbounds begin 

        cgomat_muiu  = @view cg_o_fullmatint[(I_mu,I_ij,I_u)][:,:,:]
        cgsmat_muiu  = @view cg_s_fullmatint[(S_mu,S_ij,S_u)][:,:,:]

        cgomat_nujv  = @view cg_o_fullmatint[(I_nu,I_ij,I_v)][:,:,:]
        cgsmat_nujv  = @view cg_s_fullmatint[(S_nu,S_ij,S_v)][:,:,:]

        cgomat_avu   = @view cg_o_fullmatint[(I_a,I_v,I_u)][:,:,:]
        cgsmat_avu   = @view cg_s_fullmatint[(S_a,S_v,S_u)][:,:,:]

        cgomat_anumu = @view cg_o_fullmatint[(I_a,I_nu,I_mu)][:,:,:]
        cgsmat_anumu = @view cg_s_fullmatint[(S_a,S_nu,S_mu)][:,:,:]

    end

    # orbital sum 
    sum_o::ComplexF64 = zero(ComplexF64)
    @inbounds for i_u  in 1:dI_u,
                  i_up in 1:dI_u,
                  i_v  in 1:dI_v,
                  i_ij in 1:dI_ij,
                  i_mu in 1:dI_mu,
                  i_nu in 1:dI_nu,
                  i_a  in 1:dI_a

        sum_o += conj(cgomat_avu[i_a,i_v,i_up])*
                 conj(cgomat_muiu[i_mu,i_ij,i_u])*
                 cgomat_nujv[i_nu,i_ij,i_v]*
                 cgomat_anumu[i_a,i_nu,i_mu]

    end
    # spin sum 
    sum_s::ComplexF64 = zero(ComplexF64)
    @inbounds for si_u  in 1:dS_u,
                  si_up in 1:dS_u,
                  si_v  in 1:dS_v,
                  si_ij in 1:dS_ij,
                  si_mu in 1:dS_mu,
                  si_nu in 1:dS_nu,
                  si_a  in 1:dS_a

        sum_s += conj(cgsmat_avu[si_a,si_v,si_up])*
                 conj(cgsmat_muiu[si_mu,si_ij,si_u])*
                 cgsmat_nujv[si_nu,si_ij,si_v]*
                 cgsmat_anumu[si_a,si_nu,si_mu]

    end

    return (sum_o*sum_s/(dI_u*dS_u))::ComplexF64

end
function compute_CG_Csum_orbital( 
            oindex2dimensions::Vector{Int64} ,
            cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            I_utp::Int64,
            I_vtp::Int64,
            I_ijt::Int64,
            I_mut::Int64,
            I_nut::Int64,
            I_a::Int64)
    
    # irrep quantum numbers
    I_u  = I_utp
    I_v  = I_vtp
    I_ij = I_ijt
    I_mu = I_mut
    I_nu = I_nut

    # irrep dimensions 
    dI_u::Int64  = oindex2dimensions[I_u] 
    dI_v::Int64  = oindex2dimensions[I_v]
    dI_ij::Int64 = oindex2dimensions[I_ij] 
    dI_mu::Int64 = oindex2dimensions[I_mu] 
    dI_nu::Int64 = oindex2dimensions[I_nu]
    dI_a::Int64  = oindex2dimensions[I_a]

    # clebsch-gordan matrices
    @inbounds begin 
        cgomat_muiu  = @view cg_o_fullmatint[(I_mu,I_ij,I_u)][:,:,:]
        cgomat_nujv  = @view cg_o_fullmatint[(I_nu,I_ij,I_v)][:,:,:]
        cgomat_avu   = @view cg_o_fullmatint[(I_a,I_v,I_u)][:,:,:]
        cgomat_anumu = @view cg_o_fullmatint[(I_a,I_nu,I_mu)][:,:,:]
    end

    # orbital sum 
    sum_o::ComplexF64 = zero(ComplexF64)
    #@inbounds for i_u  in 1:dI_u,
    #              i_up in 1:dI_u,
    #              i_v  in 1:dI_v,
    #              i_ij in 1:dI_ij,
    #              i_mu in 1:dI_mu,
    #              i_nu in 1:dI_nu,
    #              i_a  in 1:dI_a

    #    sum_o += cgomat_avu[i_a,i_v,i_up]*
    #             conj(cgomat_muiu[i_mu,i_ij,i_u])*
    #             cgomat_nujv[i_nu,i_ij,i_v]*
    #             conj(cgomat_anumu[i_a,i_nu,i_mu])

    #end
    @inbounds for i_u  in 1:dI_u,
                  i_v  in 1:dI_v,
                  i_ij in 1:dI_ij,
                  i_mu in 1:dI_mu,
                  i_nu in 1:dI_nu,
                  i_a  in 1:dI_a

        sum_o += cgomat_avu[i_a,i_v,i_u]*
                 conj(cgomat_muiu[i_mu,i_ij,i_u])*
                 cgomat_nujv[i_nu,i_ij,i_v]*
                 conj(cgomat_anumu[i_a,i_nu,i_mu])

    end
    
    return (sum_o/dI_u)::ComplexF64

end
function compute_CG_Csum_spin( 
            cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            S_utp::Int64,  
            S_vtp::Int64,
            S_ijt::Int64,
            S_mut::Int64,
            S_nut::Int64,
            S_a::Int64)
    
    # * we drop the "t" (tilde) and "p" (prima)
    # * for an easier reading

    # irrep quantum numbers
    S_u  = S_utp
    S_v  = S_vtp
    S_ij = S_ijt
    S_mu = S_mut
    S_nu = S_nut

    # irrep dimensions 
    dS_u::Int64  = S_u+1
    dS_v::Int64  = S_v+1
    dS_ij::Int64 = S_ij+1 
    dS_mu::Int64 = S_mu+1
    dS_nu::Int64 = S_nu+1
    dS_a::Int64  = S_a+1

    # clebsch-gordan matrices
    @inbounds begin 
        cgsmat_muiu  = @view cg_s_fullmatint[(S_mu,S_ij,S_u)][:,:,:]
        cgsmat_nujv  = @view cg_s_fullmatint[(S_nu,S_ij,S_v)][:,:,:]
        cgsmat_avu   = @view cg_s_fullmatint[(S_a,S_v,S_u)][:,:,:]
        cgsmat_anumu = @view cg_s_fullmatint[(S_a,S_nu,S_mu)][:,:,:]
    end

    # spin sum 
    sum_s::ComplexF64 = zero(ComplexF64)
    #@inbounds for si_u  in 1:dS_u,
    #              si_up in 1:dS_u,
    #              si_v  in 1:dS_v,
    #              si_ij in 1:dS_ij,
    #              si_mu in 1:dS_mu,
    #              si_nu in 1:dS_nu,
    #              si_a  in 1:dS_a

    #    sum_s += cgsmat_avu[si_a,si_v,si_up]*
    #             conj(cgsmat_muiu[si_mu,si_ij,si_u])*
    #             cgsmat_nujv[si_nu,si_ij,si_v]*
    #             conj(cgsmat_anumu[si_a,si_nu,si_mu])

    #end
    @inbounds for si_u  in 1:dS_u,
                  si_v  in 1:dS_v,
                  si_ij in 1:dS_ij,
                  si_mu in 1:dS_mu,
                  si_nu in 1:dS_nu,
                  si_a  in 1:dS_a

        sum_s += cgsmat_avu[si_a,si_v,si_u]*
                 conj(cgsmat_muiu[si_mu,si_ij,si_u])*
                 cgsmat_nujv[si_nu,si_ij,si_v]*
                 conj(cgsmat_anumu[si_a,si_nu,si_mu])

    end

    return (sum_s/dS_u)::ComplexF64

end

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# PRECOMPUTED CG SUM DICTIONARIES #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# B-type sums
function compute_CG_Bsumdict_o(
                oindex2dimensions::Vector{Int64},
                cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
                Mo_tot::Int64 ,
                II_a::Vector{Int64} )

    Bsum_o_dict = Dict{ NTuple{6,Int64} , ComplexF64 }()

    for I_mu in 1:Mo_tot,
        I_nu in 1:Mo_tot,
        I_i  in 1:Mo_tot,
        I_j  in 1:Mo_tot,
        I_uv in 1:Mo_tot,
        I_a  in II_a 

        haskey( cg_o_fullmatint , (I_mu,I_i,I_uv) ) || continue
        haskey( cg_o_fullmatint , (I_nu,I_j,I_uv) ) || continue

        B = compute_CG_Bsum_orbital(
                oindex2dimensions ,
                cg_o_fullmatint ,
                I_mu ,
                I_nu ,
                I_i ,
                I_j ,
                I_uv ,
                I_a )
        isapprox( B , zero(B) ) || (Bsum_o_dict[I_mu,I_nu,I_i,I_j,I_uv,I_a]=B)

    end

    return Bsum_o_dict

end
function compute_CG_Bsumdict_s(
                cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
                Ms_tot::Int64 ,
                Ms_shell::Int64 )

    Bsum_s_dict = Dict{ NTuple{6,Int64} , ComplexF64 }()

    S_a::Int64 = 1

    for S_mu in 0:Ms_shell,
        S_nu in 0:Ms_shell,
        S_i  in 0:Ms_tot,
        S_j  in 0:Ms_tot,
        S_uv in 0:Ms_tot

        haskey( cg_s_fullmatint , (S_mu,S_i,S_uv) ) || continue
        haskey( cg_s_fullmatint , (S_nu,S_j,S_uv) ) || continue

        B = compute_CG_Bsum_spin(
                cg_s_fullmatint ,
                S_mu ,
                S_nu ,
                S_i ,
                S_j ,
                S_uv ,
                S_a )
        isapprox( B , zero(B) ) || (Bsum_s_dict[S_mu,S_nu,S_i,S_j,S_uv,S_a]=B)

    end

    return Bsum_s_dict

end

# C-type sums
function compute_CG_Csumdict_o(
                oindex2dimensions::Vector{Int64} ,
                cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
                Mo_tot::Int64 ,
                II_a::Vector{Int64} )

    Csum_o_dict = Dict{ NTuple{6,Int64} , ComplexF64 }()

    for I_utp::Int64 in 1:Mo_tot, 
        I_vtp::Int64 in 1:Mo_tot, 
        I_ijt::Int64 in 1:Mo_tot, 
        I_mut::Int64 in 1:Mo_tot, 
        I_nut::Int64 in 1:Mo_tot, 
        I_a::Int64 in II_a
    
        ((I_mut,I_ijt,I_utp) in keys(cg_o_fullmatint)) || continue
        ((I_nut,I_ijt,I_vtp) in keys(cg_o_fullmatint)) || continue
        ((I_a,I_vtp,I_utp)   in keys(cg_o_fullmatint)) || continue
        ((I_a,I_nut,I_mut)   in keys(cg_o_fullmatint)) || continue
    
        c = compute_CG_Csum_orbital(
                    oindex2dimensions ,
                    cg_o_fullmatint ,
                    I_utp ,
                    I_vtp ,
                    I_ijt ,
                    I_mut ,
                    I_nut ,
                    I_a )::ComplexF64
        abs(c)≈0.0 || (Csum_o_dict[(I_utp,I_vtp,I_ijt,I_mut,I_nut,I_a)]=c)
         
    end

    return Csum_o_dict 

end
function compute_CG_Csumdict_s(
            cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            Ms_tot::Int64 ,
            Ms_shell::Int64 )
    
    S_a::Int64 = 1
    Csum_s_dict = Dict{ NTuple{6,Int64} , ComplexF64 }()

    for S_utp::Int64 in 0:Ms_tot, 
        S_vtp::Int64 in 0:Ms_tot, 
        S_ijt::Int64 in 0:Ms_tot, 
        S_mut::Int64 in 0:Ms_shell, 
        S_nut::Int64 in 0:Ms_shell
    
        ((S_mut,S_ijt,S_utp) in keys(cg_s_fullmatint)) || continue
        ((S_nut,S_ijt,S_vtp) in keys(cg_s_fullmatint)) || continue
        ((S_a,S_vtp,S_utp)   in keys(cg_s_fullmatint)) || continue
        ((S_a,S_nut,S_mut)   in keys(cg_s_fullmatint)) || continue
    
        c = compute_CG_Csum_spin(
                    cg_s_fullmatint ,
                    S_utp ,
                    S_vtp ,
                    S_ijt ,
                    S_mut ,
                    S_nut ,
                    S_a )::ComplexF64
        abs(c)≈0.0 || (Csum_s_dict[(S_utp,S_vtp,S_ijt,S_mut,S_nut,S_a)]=c)
         
    end

    return Csum_s_dict

end

function compute_CG_Ksum( 
            oindex2dimensions::Vector{Int64} ,
            cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            G_u::NTuple{3,Int64},
            G_v::NTuple{3,Int64},
            G_a::NTuple{3,Int64},
            G_munu::NTuple{3,Int64} ,
            G_i::NTuple{3,Int64},
            G_j::NTuple{3,Int64} )

    (_,I_u,S_u)      = G_u 
    (_,I_v,S_v)      = G_v 
    (_,I_a,S_a)      = G_a 
    (N_munu,I_munu,S_munu) = G_munu 
    (_,I_i,S_i)      = G_i 
    (_,I_j,S_j)      = G_j 

    Ko::ComplexF64 = compute_Ksum_orbital( 
            oindex2dimensions ,
            cg_o_fullmatint ,
            I_u ,
            I_v ,
            I_a ,
            I_munu ,
            I_i ,
            I_j )
    Ks::ComplexF64 = compute_Ksum_spin( 
            cg_s_fullmatint ,
            S_u ,
            S_v ,
            S_a ,
            S_munu ,
            S_i ,
            S_j )

    K::ComplexF64 = (-1)^N_munu * Ko * Ks

    return K
    
end
function compute_Ksum_orbital( 
            oindex2dimensions::Vector{Int64} ,
            cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            I_u::Int64,
            I_v::Int64,
            I_a::Int64,
            I_munu::Int64 ,
            I_i::Int64,
            I_j::Int64 )

    dI_u::Int64    = oindex2dimensions[I_u]
    dI_v::Int64    = oindex2dimensions[I_v]
    dI_a::Int64    = oindex2dimensions[I_a]
    dI_munu::Int64 = oindex2dimensions[I_munu]
    dI_i::Int64    = oindex2dimensions[I_i]
    dI_j::Int64    = oindex2dimensions[I_j]

    haskey( cg_o_fullmatint , (I_a,I_v,I_u) )    || (return zero(ComplexF64))
    haskey( cg_o_fullmatint , (I_a,I_j,I_i) )    || (return zero(ComplexF64))
    haskey( cg_o_fullmatint , (I_munu,I_i,I_u) ) || (return zero(ComplexF64))
    haskey( cg_o_fullmatint , (I_munu,I_j,I_v) ) || (return zero(ComplexF64))

    cgomat_muiu = @view cg_o_fullmatint[(I_munu,I_i,I_u)][:,:,:] 
    cgomat_nujv = @view cg_o_fullmatint[(I_munu,I_j,I_v)][:,:,:]
    cgomat_avu  = @view cg_o_fullmatint[(I_a,I_v,I_u)][:,:,:]
    cgomat_aji  = @view cg_o_fullmatint[(I_a,I_j,I_i)][:,:,:]

    osum::ComplexF64 = zero(ComplexF64)
    #for i_u    in 1:dI_u,
    #    i_ut   in 1:dI_u,
    #    i_v    in 1:dI_v,
    #    i_a    in 1:dI_a,
    #    i_munu in 1:dI_munu,
    #    i_i    in 1:dI_i,
    #    i_j    in 1:dI_j 

    #    osum += cgomat_avu[i_a,i_v,i_ut]*
    #            conj(cgomat_muiu[i_munu,i_i,i_u])*
    #            cgomat_nujv[i_munu,i_j,i_v]*
    #            conj(cgomat_aji[i_a,i_j,i_i])
    #    
    #end
    @inbounds for i_u    in 1:dI_u,
                  i_v    in 1:dI_v,
                  i_a    in 1:dI_a,
                  i_munu in 1:dI_munu,
                  i_i    in 1:dI_i,
                  i_j    in 1:dI_j 

        osum += cgomat_avu[i_a,i_v,i_u]*
                conj(cgomat_muiu[i_munu,i_i,i_u])*
                cgomat_nujv[i_munu,i_j,i_v]*
                conj(cgomat_aji[i_a,i_j,i_i])
        
    end

    return (osum/dI_u)::ComplexF64
    
end
function compute_Ksum_spin( 
            cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            S_u::Int64,
            S_v::Int64,
            S_a::Int64,
            S_munu::Int64 ,
            S_i::Int64,
            S_j::Int64 )

    dS_u::Int64    = S_u+1
    dS_v::Int64    = S_v+1
    dS_a::Int64    = S_a+1
    dS_munu::Int64 = S_munu+1
    dS_i::Int64    = S_i+1
    dS_j::Int64    = S_j+1

    haskey( cg_s_fullmatint , (S_a,S_v,S_u) )    || (return zero(ComplexF64))
    haskey( cg_s_fullmatint , (S_a,S_j,S_i) )    || (return zero(ComplexF64))
    haskey( cg_s_fullmatint , (S_munu,S_i,S_u) ) || (return zero(ComplexF64))
    haskey( cg_s_fullmatint , (S_munu,S_j,S_v) ) || (return zero(ComplexF64))

    cgsmat_muiu = @view cg_s_fullmatint[(S_munu,S_i,S_u)][:,:,:]
    cgsmat_nujv = @view cg_s_fullmatint[(S_munu,S_j,S_v)][:,:,:]
    cgsmat_avu  = @view cg_s_fullmatint[(S_a,S_v,S_u)][:,:,:]
    cgsmat_aji  = @view cg_s_fullmatint[(S_a,S_j,S_i)][:,:,:]

    ssum::ComplexF64 = zero(ComplexF64)
    #for si_u    in 1:dS_u,
    #    si_ut   in 1:dS_u,
    #    si_v    in 1:dS_v,
    #    si_a    in 1:dS_a,
    #    si_munu in 1:dS_munu,
    #    si_i    in 1:dS_i,
    #    si_j    in 1:dS_j 

    #    ssum += cgsmat_avu[si_a,si_v,si_ut]*
    #            conj(cgsmat_muiu[si_munu,si_i,si_u])*
    #            cgsmat_nujv[si_munu,si_j,si_v]*
    #            conj(cgsmat_aji[si_a,si_j,si_i])
    #    
    #end
    @inbounds for si_u    in 1:dS_u,
                  si_v    in 1:dS_v,
                  si_a    in 1:dS_a,
                  si_munu in 1:dS_munu,
                  si_i    in 1:dS_i,
                  si_j    in 1:dS_j 

        ssum += cgsmat_avu[si_a,si_v,si_u]*
                conj(cgsmat_muiu[si_munu,si_i,si_u])*
                cgsmat_nujv[si_munu,si_j,si_v]*
                conj(cgsmat_aji[si_a,si_j,si_i])
        
    end

    return (ssum/dS_u)::ComplexF64
    
end
function compute_Kdict_orbital( 
            oindex2dimensions::Vector{Int64},
            cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            Mo_tot::Int64 ,
            II_a::Vector{Int64} )

    Kdict_orbital::Dict{ NTuple{6,Int64} , ComplexF64 } = Dict()

    for I_u in 1:Mo_tot,
        I_v in 1:Mo_tot,
        I_i in 1:Mo_tot,
        I_j in 1:Mo_tot,
        I_munu in 1:Mo_tot,
        I_a in II_a

        Kdict_orbital[(I_u,I_v,I_a,I_munu,I_i,I_j)] = 
            compute_Ksum_orbital(
                    oindex2dimensions ,
                    cg_o_fullmatint ,
                    I_u ,
                    I_v ,
                    I_a ,
                    I_munu ,
                    I_i , 
                    I_j )

    end

    return Kdict_orbital

end
function compute_Kdict_spin( 
            cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            Ms_tot::Int64 )

    S_a = 1

    Kdict_spin::Dict{ NTuple{6,Int64} , ComplexF64 } = Dict()

    for S_u    in 0:Ms_tot,
        S_v    in 0:Ms_tot,
        S_i    in 0:Ms_tot,
        S_j    in 0:Ms_tot,
        S_munu in 0:Ms_tot

        Kdict_spin[(S_u,S_v,S_a,S_munu,S_i,S_j)] = 
            compute_Ksum_spin(
                    cg_s_fullmatint ,
                    S_u ,
                    S_v ,
                    S_a ,
                    S_munu ,
                    S_i , 
                    S_j )

    end

    return Kdict_spin

end
function Kdict2array_orbital( Kdict_orbital::Dict{ NTuple{6,Int64} , ComplexF64 } )

    Mo_u    = maximum( k[1] for k in keys(Kdict_orbital) )
    Mo_v    = maximum( k[2] for k in keys(Kdict_orbital) )
    Mo_a    = maximum( k[3] for k in keys(Kdict_orbital) )
    Mo_munu = maximum( k[4] for k in keys(Kdict_orbital) )
    Mo_i    = maximum( k[5] for k in keys(Kdict_orbital) )
    Mo_j    = maximum( k[6] for k in keys(Kdict_orbital) )

    Karray_orbital = zeros( ComplexF64 , Mo_u , Mo_v , Mo_a , Mo_munu , Mo_i , Mo_j )

    for (Ikey,coeff) in Kdict_orbital 

        Karray_orbital[Ikey...] = coeff

    end

    return Karray_orbital
    
end
function Kdict2array_spin( Kdict_spin::Dict{ NTuple{6,Int64} , ComplexF64 } )

    Ms_u    = maximum( k[1] for k in keys(Kdict_spin) )+1 
    Ms_v    = maximum( k[2] for k in keys(Kdict_spin) )+1
    Ms_a    = maximum( k[3] for k in keys(Kdict_spin) )+1
    Ms_munu = maximum( k[4] for k in keys(Kdict_spin) )+1
    Ms_i    = maximum( k[5] for k in keys(Kdict_spin) )+1
    Ms_j    = maximum( k[6] for k in keys(Kdict_spin) )+1

    Karray_spin = zeros( ComplexF64 , Ms_u , Ms_v , Ms_a , Ms_munu , Ms_i , Ms_j )

    for (SSSSSS,coeff) in Kdict_spin

        Skey = SSSSSS.+1
        Karray_spin[Skey...] = coeff

    end

    return Karray_spin
    
end
function compute_Ksum_arrays(
            oindex2dimensions::Vector{Int64},
            cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            Mo_tot::Int64 ,
            II_a::Vector{Int64} ,
            Ms_tot::Int64 ,
            Ms_shell::Int64 )

    Kdict_orbital = compute_Kdict_orbital(
                        oindex2dimensions ,
                        cg_o_fullmatint ,
                        Mo_tot ,
                        II_a )
    Kdict_spin    = compute_Kdict_spin( 
                        cg_s_fullmatint ,
                        Ms_tot )

    Karray_orbital = Kdict2array_orbital( Kdict_orbital )
    Karray_spin    = Kdict2array_spin( Kdict_spin )

    return (Karray_orbital,Karray_spin)

end






