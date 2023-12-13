using LinearAlgebra
using DelimitedFiles

# ============ #
# FILE WRITING #
# ============ #

# file header 
const spectralheader = "# 1 energy | 2 spectral function\n"

# Naming convention
function spectral_filename( 
            label::String ;
            z::Float64=0.0 ,
            zavg::Bool=false ,
            orbital::Int64=0 ,
            tail::String="" )

    z_string = zavg ? "_zavg" : "_z$(z)"
    o_string = orbital==0 ? "" : "_o$(orbital)"
    return "spectral/spectral_$(label)$(z_string)$(o_string)$(tail).dat"

end

# Writing to file
function write_spectral_function(
            filename::String ,
            spectral_data::Matrix{Float64} ;
            header::String=spectralheader ,
            orbitalresolved_header::String="" )
    
    open( filename , write=true ) do f 
        if orbitalresolved_header!==""
            write( f , orbitalresolved_header )
        end
        write( f , spectralheader )
        writedlm( f , spectral_data )
    end

end

# Reading from file
function read_spectral_data( filename::String )
    return readdlm( filename , comments=true , comment_char='#' )
end

###############

function is_in_interval( omega , emin , emax ) 
    return (omega>=emin && omega<=emax)
end

# initial M matrix
function get_M0( pcg , qq_a )
    M::Dict{ NTuple{6,Int64} , Dict{NTuple{2,NTuple{6,Int64}},ComplexF64} } = Dict()
    for q_a in qq_a 
        M[q_a] = Dict( (k[3],k[1])=>conj(v) for (k,v) in pcg if k[2]==q_a ) 
    end
    return M 
end

# update M (dict)
function get_new_M( 
            M_old::Dict{ NTuple{6,Int64}, Dict{Tuple{NTuple{6,Int64},NTuple{6,Int64 }}, ComplexF64}},
            qq_a::Vector{NTuple{6,Int64}},
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            oindex2dimensions::Vector{ Int64 },
            combinations_uprima::Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} },
            cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} };
            verbose=false )
    
    M_new::Dict{ NTuple{6,Int64}, Dict{Tuple{NTuple{6,Int64},NTuple{6,Int64 }}, ComplexF64}} = Dict( q_a=>Dict() for q_a in qq_a )

    # iterate over annihilation operators 
    # and irreps of all irreps
    #@inbounds for q_a::NTuple{6,Int64} in qq_a,
    #             (G_u::NTuple{3,Int64},(E_u,U_u)::Tuple{Vector{Float64},Matrix{ComplexF64}}) in irrEU,
    #             (G_v::NTuple{3,Int64},(E_v,U_v)::Tuple{Vector{Float64},Matrix{ComplexF64}}) in irrEU
    @inbounds for q_a in qq_a,
                 (G_u,(E_u,U_u)) in irrEU,
                 (G_v,(E_v,U_v)) in irrEU

        if verbose 
            println( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" )
            println( "Irrep combination: $G_u & $G_v" )
            println( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" )
            println()
        end

        # irrep quantum numbers
        (N_u,I_u,S_u) = G_u
        (N_v,I_v,S_v) = G_v

        N_u==(N_v-1) || continue

        # within each irrep, iterate over 
        # multiplets and partners
        @inbounds for r_u::Int64 in 1:length(E_u),
                      r_v::Int64 in 1:length(E_v),
                      i_u::Int64 in 1:oindex2dimensions[I_u]::Int64,
                      i_v::Int64 in 1:oindex2dimensions[I_v]::Int64,
                      (si_u::Int64,s_u::Int64) in enumerate((-S_u):2:S_u),
                      (si_v::Int64,s_v::Int64) in enumerate((-S_v):2:S_v)

            # pair quantum numbers
            q_u = (G_u...,i_u,s_u,r_u)
            q_v = (G_v...,i_v,s_v,r_v)

            if verbose 
                println( "Q pair: $q_u & $q_v" )
                println( "..............................." )
                println()
            end

            # M(u,v)
            M_loc = zero(ComplexF64)

            # within irreps G_u and G_v, iterate over
            # multiplets m_u' and m_v' in order to 
            # undo the diagonalization
            combs_Gup = combinations_uprima[G_u]::Vector{NTuple{3,NTuple{4,Int64}}} 
            combs_Gvp = combinations_uprima[G_v]::Vector{NTuple{3,NTuple{4,Int64}}} 
            #@inbounds for (m_up,m_mu,m_i)::NTuple{3,NTuple{4,Int64}} in combs_Gup,
            #              (m_vp,m_nu,m_j)::NTuple{3,NTuple{4,Int64}} in combs_Gvp
            @inbounds for (m_up,m_mu,m_i) in combs_Gup,
                          (m_vp,m_nu,m_j) in combs_Gvp

                # outer shell states must coincide,
                # therefore so do multiplets
                m_mu==m_nu || continue
                
                if verbose 
                    println( "multiplet decomposition:" )
                    @show m_up, m_mu, m_i
                    @show m_vp, m_nu, m_j
                    println()
                end

                # previous-step quantum numbers
                (N_mu,I_mu,S_mu,r_mu) = m_mu 
                (N_i, I_i ,S_i ,r_i ) = m_i 
                (N_j, I_j ,S_j ,r_j ) = m_j 
                
                # permutation factor 
                fac = (-1)^N_mu

                #((N_i==N_j+1) || (N_i==N_j-1)) || continue 
                (N_i==N_j-1) || continue 
                
                # clebsch-gordan matrices
                @inbounds begin
                    cgomat_muiu = @view cg_o_fullmatint[I_mu,I_i,I_u][:,:,i_u]
                    cgomat_mujv = @view cg_o_fullmatint[I_mu,I_j,I_v][:,:,i_v]
                    cgsmat_muiu = @view cg_s_fullmatint[S_mu,S_i,S_u][:,:,si_u]
                    cgsmat_mujv = @view cg_s_fullmatint[S_mu,S_j,S_v][:,:,si_v]
                end

                # outer multiplicity for u' and v'
                r_up = m_up[end] 
                r_vp = m_vp[end]

                # diagonalization matrix U term 
                uterm = conj(U_u[r_up,r_u])*U_v[r_vp,r_v]
                isapprox(abs2(uterm),0.0) && continue
                if verbose
                    println( "uterm test passed" )
                    println()
                end

                # iterate over partners
                @inbounds for i_i::Int64  in 1:oindex2dimensions[I_i ]::Int64,
                              i_j::Int64  in 1:oindex2dimensions[I_j ]::Int64,
                              i_mu::Int64 in 1:oindex2dimensions[I_mu]::Int64,
                              (si_i::Int64,s_i::Int64)    in enumerate((-S_i):2:S_i),
                              (si_j::Int64,s_j::Int64)    in enumerate((-S_j):2:S_j),
                              (si_mu::Int64,s_mu::Int64)  in enumerate((-S_mu):2:S_mu)

                    # i and j quantum numbers 
                    q_i = (N_i,I_i,S_i,i_i,s_i,r_i)
                    q_j = (N_j,I_j,S_j,i_j,s_j,r_j)

                    if verbose
                        println( "q_i = $q_i" )
                        println( "q_a = $q_a" ) 
                        println( "q_j = $q_j" )
                        println()
                    end

                    # old M coefficient
                    Mold_coeff = get( M_old[q_a] ,
                                      (q_i,q_j) , 
                                      zero(ComplexF64) )::ComplexF64
                    Mold_coeff==zero(Mold_coeff) && continue
                    if verbose
                        println( "mold coeff test passed" )
                    end

                    # clebsch-gordan coefficients 
                    cg_muiu = cgomat_muiu[i_mu,i_i]*cgsmat_muiu[si_mu,si_i]
                    cg_mujv = cgomat_mujv[i_mu,i_j]*cgsmat_mujv[si_mu,si_j]
                    cg_total = conj(cg_muiu)*cg_mujv
                    cg_total==zero(cg_total) && continue
                    if verbose
                        println( "cg test passed" )
                    end

                    # compound coefficient 
                    M_loc += fac * uterm * cg_total * Mold_coeff

                    if verbose 
                        contrib = uterm * cg_total * Mold_coeff
                        println( "contribution: $contrib" )
                        println()
                    end
                end
            end
            isapprox(abs2(M_loc),0.0) || (M_new[q_a][q_u,q_v]=M_loc)
        end
    end

    return M_new 

end

# get A as
#   { q_a => { m => (E[m],[coeff-,coeff+]) } }
function get_A( M::Dict{ NTuple{6,Int64} , Dict{NTuple{2,NTuple{6,Int64}},ComplexF64} },
                irrEU::Dict{ Tuple{Int64,Int64,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} } )

    A::Dict{ NTuple{6,Int64}, Dict{ NTuple{4,Int64}, Tuple{Float64,Vector{Float64}}}} = Dict( q_a=>Dict((G...,r)=>(E[r],[0.0,0.0]) for (G,(E,U)) in irrEU for r=1:length(E)) for q_a in qq_a )

    @inbounds for q_a in keys(M),
                  ((q_u,q_v),matel) in M[q_a]

        G_u = q_u[1:3]
        G_v = q_v[1:3]

        m_u = (G_u...,q_u[6])
        m_v = (G_v...,q_v[6])
        
        r_u = q_u[6]
        r_v = q_v[6]

        E_u = irrEU[G_u][1][r_u]
        E_v = irrEU[G_v][1][r_v]

        E_v==0.0 && (A[q_a][m_u][2][1] += abs2(matel))
        E_u==0.0 && (A[q_a][m_v][2][2] += abs2(matel))

    end

    for (q_a,v) in A, 
        (m,(e,coeffs)) in v 

        sum(coeffs)==0.0 && pop!( A[q_a] , m )

    end

    return A
end

# pcgmat method
function NRG_mat_spectral( 
            iterations::Int64, 
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
            mm_i ,
            M::Dict{ NTuple{6,Int64}, Dict{Tuple{NTuple{6,Int64},NTuple{6,Int64 }}, ComplexF64}} ,
            AA::Vector{ Dict{NTuple{6,Int64},Dict{NTuple{4,Int64},Tuple{Float64,Vector{Float64}}}} } ;
            spinarray=true ,
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

    println( "M size in NRG: $(length(collect(keys(M[qq_a[1]]))))" )

    temperatures   = []
    magnetizations = []
    energies       = []
    numbers        = []
    partitions     = []
    entropies      = []
    impspins       = []
    impnums        = []
    impmults       = []
    partitions0    = []

    performance = []
    maxes = []
    eigenvals5::Matrix{Float64} = zeros(Float64,iterations-1,5)

    ξ = compute_xi_vector( iterations , z , L ; discretization=discretization )

    for n=2:iterations

        println( "ITERATION n=$n" )

        # cut off and block block multiplets
        print( "Applying cutoff to obtain block multiplets... " )
        @time (multiplets_block, discarded) = 
                    cut_off!( irrEU ; 
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
        if distributed
            # spinarray must be true
            ppp = @timed (irrEU,combinations_uprima) = 
                        matdiag_distributed_mat( 
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
            global ppp = @timed (irrEU,combinations_uprima) = 
                        matdiag_serial_mat( 
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
        push!( impmults , m_imp )

        # spectral thermo 
        G0 = [G for (G,(E,U)) in irrEU if E[1]==0][1]
        I0,S0 = G0[2:3]
        D0s = S0+1
        D0o = oindex2dimensions[I0]
        part0 = D0o*D0s
        push!( partitions0 , part0 )

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

        println( "M size before: $(length(collect(keys(M[qq_a[1]]))))" )
        M = get_new_M( M , 
                       qq_a , 
                       irrEU , 
                       oindex2dimensions , 
                       combinations_uprima ,
                       cg_o_fullmatint , 
                       cg_s_fullmatint )
        println( "M size after: $(length(collect(keys(M[qq_a[1]]))))" )
        push!( AA , get_A(M,irrEU) )

    end

    maxe_avg = sum(maxes)/length(maxes)
    @show maxe_avg

    return ( t=temperatures , 
             m=magnetizations , 
             e=energies , 
             p=partitions ,
             n=numbers , 
             entr=entropies,
             perf=performance ,
             impmults=impmults ,
             AA=AA ,
             maxes=maxes ,
             partitions0 )
end

# pcgmat method with reduced M
function NRG_mat_spectral_redM( 
            iterations::Int64, 
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
            mm_i ,
            Mred,
            AA;
            spinarray=true ,
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
    partitions0    = []

    performance = []
    maxes = []
    eigenvals5::Matrix{Float64} = zeros(Float64,iterations-1,5)

    multiplets_a = collect(Set( (q[1:3]...,q[6]) for q in qq_a ))

    ξ = compute_xi_vector( iterations , z , L ; discretization=discretization )

    for n=2:iterations

        println( "ITERATION n=$n" )

        # cut off and block block multiplets
        print( "Applying cutoff to obtain block multiplets... " )
        @time (multiplets_block, discarded) = 
                    cut_off!( irrEU ; 
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
        if distributed
            # spinarray must be true
            ppp = @timed (irrEU,combinations_uprima) = 
                        matdiag_distributed_mat( 
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
            global ppp = @timed (irrEU,combinations_uprima) = 
                        matdiag_serial_mat( 
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

        Mred, AA = update_redmat_AA(
                    Mred ,
                    irrEU ,
                    combinations_uprima ,
                    multiplets_a ,
                    cg_o_fullmatint ,
                    cg_s_fullmatint ,
                    AA )

    end

    maxe_avg = sum(maxes)/length(maxes)
    @show maxe_avg

    return ( t=temperatures , 
             m=magnetizations , 
             e=energies , 
             p=partitions ,
             n=numbers , 
             entr=entropies,
             perf=performance ,
             impmults=impmults ,
             AA=AA ,
             maxes=maxes )
end

# broadened deltas
function P( w_pm_E::N , 
            eta::N ; 
            distribution="gaussian" , 
            omega=0.0 , 
            E=0.0 ) where {N<:Real}

    if distribution=="gaussian" 
        return P_Gaussian( w_pm_E , eta ) 
    elseif distribution=="loggaussian" 
        return P_LogGaussian( w_pm_E , eta , omega , E )
    end
end
function P_Gaussian( w_pm_E::N , eta::N ) where {N<:Real}
    return exp(-(w_pm_E/eta)^2) / (eta*sqrt(pi))
end
function P_LogGaussian( w_pm_E::N , eta::N , omega::N , E::N ) where {N<:Real}
    iszero(E) && (return 0.0)
    return exp(-eta^2/4.0)/(eta*E*sqrt(pi)) * exp(-(log(abs(omega)/E)/eta)^2)
end


# given reduced M, computes A
function redM2A( 
        Mred::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } ,
        multiplets_a::Vector{NTuple{4,Int64}} ,
        cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
        partition::Int64 ;
        verbose=false ,
        density_matrix::Dict{IntIrrep,Matrix{Float64}}=Dict{IntIrrep,Matrix{Float64}}() ,
        multiplets_kept::Set{IntMultiplet}=Set{IntMultiplet}() )

    if verbose 
        println( "~~~~~~~~~~~~~~~~~~" )
        println( "COMPUTING A MATRIX" )
        println( "~~~~~~~~~~~~~~~~~~" )
        println()
    end

    use_density_matrix = length(density_matrix)!==0

    irreps_kept = Set( m[1:3] for m in multiplets_kept )

    # m => ( E_m , [ weight- , weight+ ] )
    A = Dict( (G...,r)=>[E[r],[0.0,0.0]]
              for (G,(E,U)) in irrEU 
              for r in 1:length(E) )

    # ground multiplet 
    GG0::Set{NTuple{3,Int64}} = Set( G for (G,(E,U)) in irrEU if iszero(E[1]) )

    if verbose 
        println( "ground irreps: $GG0")
        println()
    end

    # iterate over irreducible matrices
    if verbose 
        println( "iterating over reduced matrix..." )
        println()
    end
    if !use_density_matrix
        for ((G_u,G_a,G_v),redmat) in Mred 

            psector = ((G_u in GG0) && (G_u!==G_v))
            nsector = ((G_v in GG0) && (G_v!==G_u))
            (psector || nsector) || continue

            if verbose 
                Gcomb = (G_u,G_a,G_v)
                sector = psector ? "positive" : "negative"
                println( "$Gcomb ==> $sector sector" )
                println()
            end

            # irrep quantum numbers
            (N_u,I_u,S_u) = G_u 
            (N_a,I_a,S_a) = G_a
            (N_v,I_v,S_v) = G_v

            ((I_a,I_v,I_u) in keys(cg_o_fullmatint)) || continue
            ((S_a,S_v,S_u) in keys(cg_s_fullmatint)) || continue

            # clebsch-gordan contribution and partition function
            @views begin
                cgomat = cg_o_fullmatint[(I_a,I_v,I_u)]
                cgsmat = cg_s_fullmatint[(S_a,S_v,S_u)]
            end
            dim_Ia,dim_Iv,dim_Iu = size(cgomat)
            dim_Sa,dim_Sv,dim_Su = size(cgsmat)
            cg = sum(abs2.(cgomat))*
                 sum(abs2.(cgsmat))/
                 (dim_Ia*dim_Sa)


            # iterate over multiplets in irrep combination
            for rrr in CartesianIndices(redmat)

                # outer multiplicities
                (r_u,r_a,r_v) = Tuple(rrr)

                # multiplets
                m_u = (G_u...,r_u)
                m_a = (G_a...,r_a)
                m_v = (G_v...,r_v)

                # coefficient from reduced matrix element
                redmatel = abs2(redmat[rrr])

                # total contribution 
                w = cg*redmatel

                # insert in A
                if (nsector && r_v==1)

                    A[m_u][2] .+= [w/partition,0.0]

                end
                if (psector && r_u==1)

                    A[m_v][2] .+= [0.0,w/partition]

                end
            end
        end
    elseif use_density_matrix
        for ((G_u,G_a,G_v),redmat) in Mred

            # irrep quantum numbers
            (_,I_u,S_u) = G_u 
            (_,I_a,S_a) = G_a
            (_,I_v,S_v) = G_v

            ((I_a,I_v,I_u) in keys(cg_o_fullmatint)) || continue
            ((S_a,S_v,S_u) in keys(cg_s_fullmatint)) || continue
            Gu_is_ground = haskey(density_matrix,G_u)
            Gv_is_ground = haskey(density_matrix,G_v)
            (Gu_is_ground || Gv_is_ground) || continue

            Gu_is_ground && (@views dm_u=density_matrix[G_u])
            Gv_is_ground && (@views dm_v=density_matrix[G_v])

            Ru_kept = Gu_is_ground ? size(dm_u,1) : 0
            Rv_kept = Gv_is_ground ? size(dm_v,1) : 0

            @show Gu_is_ground, Ru_kept
            @show Gv_is_ground, Rv_kept
            @show size(redmat)
            println()

            # clebsch-gordan contribution and partition function
            @views begin
                cgomat = cg_o_fullmatint[(I_a,I_v,I_u)]
                cgsmat = cg_s_fullmatint[(S_a,S_v,S_u)]
            end
            dim_Ia,dim_Iv,dim_Iu = size(cgomat)
            dim_Sa,dim_Sv,dim_Su = size(cgsmat)
            cg = (dim_Iu*dim_Su)/(dim_Ia*dim_Sa)

            # iterate over multiplets in irrep combination
            for r1_u in axes(redmat,1),
                r1_v in axes(redmat,3),
                r2_u in axes(redmat,1),
                r2_v in axes(redmat,3),
                r_a  in axes(redmat,2)

                # multiplets
                m_u = (G_u...,r1_u)
                m_v = (G_v...,r1_v)

                # insert in A
                # G_v as the kept "ground" multiplet, G_u discarded
                if Gv_is_ground && (r1_u==r2_u) && r1_u>Ru_kept && r1_v<=Rv_kept && r2_v<=Rv_kept

                    # reduced excitation matrix element product
                    redmatel = conj(redmat[r1_u,r_a,r1_v])*redmat[r1_u,r_a,r2_v]
                    w = cg*redmatel
                    @show redmatel
                    @show length(irrEU[G_u][1]),r1_u
                    @show length(irrEU[G_v][1]),r1_v,r2_v

                    A[m_u][1] = irrEU[G_u][1][r1_u] - 0.5*(irrEU[G_v][1][r1_v]+irrEU[G_v][1][r2_v])
                    A[m_u][2] .+= [w*dm_v[r1_v,r2_v],0.0]

                end
                # G_u as the kept "ground" multiplet, G_v discarded
                if Gu_is_ground && (r1_v==r2_v) && r1_v>Rv_kept  && r1_u<=Ru_kept && r2_u<=Ru_kept

                    # reduced excitation matrix element product
                    redmatel = redmat[r1_u,r_a,r1_v]*conj(redmat[r2_u,r_a,r1_v])
                    w = cg*redmatel

                    A[m_v][1] = irrEU[G_v][1][r1_v] - 0.5*(irrEU[G_u][1][r1_u]+irrEU[G_u][1][r2_u])
                    A[m_v][2] .+= [0.0,w*dm_u[r1_u,r2_u]]

                end
            end
        end
    end

    return A
end
function redM2A_selfenergy( 
        Mred::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } ,
        Mred_se::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } ,
        multiplets_a::Vector{NTuple{4,Int64}} ,
        cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
        partition::Int64 ;
        verbose=false )

    if verbose 
        println( "~~~~~~~~~~~~~~~~~~" )
        println( "COMPUTING A MATRIX" )
        println( "~~~~~~~~~~~~~~~~~~" )
        println()
    end

    # m => ( E_m , [ weight- , weight+ ] )
    A::Dict{ IntMultiplet , Tuple{Float64,Vector{Float64}} } = 
        Dict( (G...,r)=>(E[r],[0.0,0.0]) 
              for (G,(E,U)) in irrEU 
              for r in 1:length(E) )

    # ground multiplet 
    #G0::NTuple{3,Int64} = first(collect( G for (G,(E,U)) in irrEU if isapprox( E[1] , zero(E[1]) ) ))
    GG0::Set{NTuple{3,Int64}} = Set( G for (G,(E,U)) in irrEU if iszero(E[1]) )
    println( "A matrix" )
    @show GG0

    if verbose 
        println( "ground irreps: $GG0")
        println()
    end

    # iterate over irreducible matrices
    if verbose 
        println( "iterating over reduced matrix..." )
        println()
    end
    for ((G_u,G_a,G_v),redmat) in Mred 

        haskey( Mred_se , (G_u,G_a,G_v) ) || continue
        redmat_se = Mred_se[G_u,G_a,G_v]

        # sector: negative=false, positive=true
        #if ((G_u in GG0) && !(G_v in GG0))
        #    psector = true
        #elseif ((G_v in GG0) && !(G_u in GG0))
        #    psector = false 
        #else
        #    continue
        #end
        psector = ((G_u in GG0) && (G_u!==G_v))
        nsector = ((G_v in GG0) && (G_v!==G_u))
        (psector || nsector) || continue

        if verbose 
            Gcomb = (G_u,G_a,G_v)
            sector = psector ? "positive" : "negative"
            println( "$Gcomb ==> $sector sector" )
            println()
        end

        # irrep quantum numbers
        (N_u,I_u,S_u) = G_u 
        (N_a,I_a,S_a) = G_a
        (N_v,I_v,S_v) = G_v

        ((I_a,I_v,I_u) in keys(cg_o_fullmatint)) || continue
        ((S_a,S_v,S_u) in keys(cg_s_fullmatint)) || continue

        # clebsch-gordan contribution and partition function
        @views begin
            cgomat = cg_o_fullmatint[(I_a,I_v,I_u)]
            cgsmat = cg_s_fullmatint[(S_a,S_v,S_u)]
        end
        dim_Ia,dim_Iv,dim_Iu = size(cgomat)
        dim_Sa,dim_Sv,dim_Su = size(cgsmat)
        cg = sum(abs2.(cgomat))*
             sum(abs2.(cgsmat))/
             (dim_Ia*dim_Sa)


        # iterate over multiplets in irrep combination
        for rrr in CartesianIndices(redmat)

            # outer multiplicities
            (r_u,r_a,r_v) = Tuple(rrr)

            # exclude non-contributing elements
            #((!psector && r_v==1) || ( psector && r_u==1)) || continue

            # multiplets
            m_u = (G_u...,r_u)
            m_a = (G_a...,r_a)
            m_v = (G_v...,r_v)

            # coefficient from reduced matrix element
            redmatel = conj(redmat_se[rrr])*redmat[rrr]

            # total contribution 
            w = cg*redmatel

            # insert in A
            if (nsector && r_v==1)

                A[m_u][2] .+= [w/partition,0.0]

            end
            if (psector && r_u==1)

                A[m_v][2] .+= [0.0,w/partition]

            end
        end
    end

    return A
end

# Compute orbital resolved A
function redM2A_orbitalresolved( 
        Mred::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } ,
        multiplets_a::Vector{NTuple{4,Int64}} ,
        cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
        irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
        partition::Int64 ;
        verbose=false )

    if verbose 
        println( "~~~~~~~~~~~~~~~~~~" )
        println( "COMPUTING A MATRIX" )
        println( "~~~~~~~~~~~~~~~~~~" )
        println()
    end

    Norbitals = size(collect(values(Mred))[1])[2]

    # orbital
    A::Dict{ NTuple{2,IntMultiplet} , Tuple{Float64,Vector{Float64}} } = Dict( 
        ((G...,r),m_a)=>(E[r],[0.0,0.0])
        for m_a in multiplets_a
        for (G,(E,U)) in irrEU
        for r in 1:length(E) 
    )

    # ground multiplet 
    #G0::NTuple{3,Int64} = first(collect( G for (G,(E,U)) in irrEU if isapprox( E[1] , zero(E[1]) ) ))
    GG0::Set{NTuple{3,Int64}} = Set( G for (G,(E,U)) in irrEU if isapprox( E[1] , zero(E[1]) ) )

    if verbose 
        println( "ground irreps: $GG0")
        println()
    end

    # iterate over irreducible matrices
    if verbose 
        println( "iterating over reduced matrix..." )
        println()
    end
    for ((G_u,G_a,G_v),redmat) in Mred 

        # sector: negative=false, positive=true
        if ((G_u in GG0) && !(G_v in GG0))
            psector = true
        elseif ((G_v in GG0) && !(G_u in GG0))
            psector = false 
        else
            continue
        end

        if verbose 
            Gcomb = (G_u,G_a,G_v)
            sector = psector ? "positive" : "negative"
            println( "$Gcomb ==> $sector sector" )
            println()
        end

        # irrep quantum numbers
        (N_u,I_u,S_u) = G_u 
        (N_a,I_a,S_a) = G_a
        (N_v,I_v,S_v) = G_v

        ((I_a,I_v,I_u) in keys(cg_o_fullmatint)) || continue
        ((S_a,S_v,S_u) in keys(cg_s_fullmatint)) || continue
        
        # clebsch-gordan contribution 
        @views begin
            cgomat = cg_o_fullmatint[(I_a,I_v,I_u)]
            cgsmat = cg_s_fullmatint[(S_a,S_v,S_u)]
        end
        d_Ia,_,d_Iu = size(cgomat)
        d_Sa,_,d_Su = size(cgsmat)
        cg = d_Iu*d_Su*d_Iu*d_Su*d_Ia*d_Sa

        # iterate over multiplets in irrep combination
        for rrr in CartesianIndices(redmat)

            # outer multiplicities
            (r_u,r_a,r_v) = Tuple(rrr)

            # exclude non-contributing elements
            ((!psector && r_v==1) || ( psector && r_u==1)) || continue

            # multiplets
            m_u = (G_u...,r_u)
            m_a = (G_a...,r_a)
            m_v = (G_v...,r_v)

            # coefficient from reduced matrix element
            redmatel = abs2(redmat[rrr])

            # total contribution 
            w = cg*redmatel

            # insert in A
            if !psector

                A[m_u,m_a][2] .+= [w/partition,0.0]

            elseif psector

                A[m_v,m_a][2] .+= [0.0,w/partition]

            end
        end
    end

    return A
end

function get_new_blockredmat( 
            ijredmat::Dict{NTuple{3,NTuple{3,Int64}},Array{ComplexF64,3}} ,
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            combinations_uprima::Dict{NTuple{3,Int64}, Vector{NTuple{3,NTuple{4,Int64}}}},
            multiplets_a::Vector{NTuple{4,Int64}} ,
            cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            oindex2dimensions::Vector{Int64})

    uvredmat::Dict{NTuple{3,NTuple{3,Int64}},Array{ComplexF64,3}} = Dict()

    # hopping irreps and multiplicities
    irrmult_a = get_irreps(Set(multiplets_a);multiplicity=true)
    
    # iterate over irrep combinations G_u and G_v
    # in the new step
    @inbounds for (G_u::NTuple{3,Int64},ucombs::Vector{NTuple{3,NTuple{4,Int64}}}) in combinations_uprima,
                  (G_v::NTuple{3,Int64},vcombs::Vector{NTuple{3,NTuple{4,Int64}}}) in combinations_uprima

        # irrep quantum numbers
        (N_u,I_u,S_u) = G_u 
        (N_v,I_v,S_v) = G_v 

        # early discard
        N_u==(N_v+1) || continue
        
        # diagonalization matrices 
        U_u = conj.(irrEU[G_u][2])::Matrix{ComplexF64}  
        U_v=        irrEU[G_v][2]::Matrix{ComplexF64} 

        # irrep multiplicities
        R_u = length(ucombs)
        R_v = length(vcombs)

        # dimension factor 
        dI_u  = oindex2dimensions[I_u]::Int64
        dS_u  = (S_u+1)::Int64
        dim_u = dI_u*dS_u
        
        # allocating clebsch-gordan matrices
        cgo_contracted::Matrix{ComplexF64} = zeros( ComplexF64 , dI_u , dI_u )
        cgs_contracted::Matrix{ComplexF64} = zeros( ComplexF64 , dS_u , dS_u )

        # iterate over hopping irreps G_a
        @inbounds for (G_a::NTuple{3,Int64},R_a::Int64) in irrmult_a

            # irrep quantum numbers
            (N_a,I_a,S_a) = G_a

            # initiate new matrix entry
            thisredmat = zeros(ComplexF64,R_u,R_a,R_v)
            transformedmat = zeros(ComplexF64,R_u,R_a,R_v)

            # iterate over primed combinations, to be 
            # combined later by u matrices.
            # here:
            #   u = u'
            #   v = v'
            @inbounds for (m_u::NTuple{4,Int64},m_mu::NTuple{4,Int64},m_i::NTuple{4,Int64}) in ucombs,
                          (m_v::NTuple{4,Int64},m_nu::NTuple{4,Int64},m_j::NTuple{4,Int64}) in vcombs

                # for block operators, same shell states
                m_mu==m_nu || continue

                # multiplet quantum numbers
                r_u = m_u[4]
                r_v = m_v[4]
                (N_mu,I_mu,S_mu,r_mu) = m_mu
                (N_i, I_i, S_i, r_i ) = m_i 
                (N_j, I_j, S_j, r_j ) = m_j
                G_mu = (N_mu,I_mu,S_mu)
                G_i = (N_i,I_i,S_i)
                G_j = (N_j,I_j,S_j)

                ((G_i,G_a,G_j) in keys(ijredmat)) || continue
                ((I_a,I_v,I_u) in keys(cg_o_fullmatint)) || continue
                ((S_a,S_v,S_u) in keys(cg_s_fullmatint)) || continue
                ((I_a,I_j,I_i) in keys(cg_o_fullmatint)) || continue
                ((S_a,S_j,S_i) in keys(cg_s_fullmatint)) || continue

                # clebsch-gordan matrices
                #cgomat_muiu = @view cg_o_fullmatint[(I_mu,I_i,I_u)][:,:,:]
                #cgsmat_muiu = @view cg_s_fullmatint[(S_mu,S_i,S_u)][:,:,:]
                #cgomat_mujv = @view cg_o_fullmatint[(I_mu,I_j,I_v)][:,:,:]
                #cgsmat_mujv = @view cg_s_fullmatint[(S_mu,S_j,S_v)][:,:,:]
                #cgomat_avu  = @view cg_o_fullmatint[(I_a, I_v,I_u)][:,:,:]
                #cgsmat_avu  = @view cg_s_fullmatint[(S_a, S_v,S_u)][:,:,:]
                #cgomat_aji  = @view cg_o_fullmatint[(I_a, I_j,I_i)][:,:,:]
                #cgsmat_aji  = @view cg_s_fullmatint[(S_a, S_j,S_i)][:,:,:]
                cgomat_muiu::Array{ComplexF64,3} = cg_o_fullmatint[(I_mu,I_i,I_u)][:,:,:]
                cgsmat_muiu::Array{ComplexF64,3} = cg_s_fullmatint[(S_mu,S_i,S_u)][:,:,:]
                cgomat_mujv::Array{ComplexF64,3} = cg_o_fullmatint[(I_mu,I_j,I_v)][:,:,:]
                cgsmat_mujv::Array{ComplexF64,3} = cg_s_fullmatint[(S_mu,S_j,S_v)][:,:,:]
                cgomat_avu::Array{ComplexF64,3}  = cg_o_fullmatint[(I_a, I_v,I_u)][:,:,:]
                cgsmat_avu::Array{ComplexF64,3}  = cg_s_fullmatint[(S_a, S_v,S_u)][:,:,:]
                cgomat_aji::Array{ComplexF64,3}  = cg_o_fullmatint[(I_a, I_j,I_i)][:,:,:]
                cgsmat_aji::Array{ComplexF64,3}  = cg_s_fullmatint[(S_a, S_j,S_i)][:,:,:]


                # clebsch-gordan sum coefficient
                #@tensor cgo_contracted[i_ut,i_u] = 
                #            conj.(cgomat_avu)[i_a,i_v,i_ut]*
                #            conj.(cgomat_muiu)[i_mu,i_i,i_u]*
                #            cgomat_mujv[i_mu,i_j,i_v]*
                #            cgomat_aji[i_a,i_j,i_i]
                #@tensor cgs_contracted[s_ut,s_u] = 
                #            conj.(cgsmat_avu)[s_a,s_v,s_ut]*
                #            conj.(cgsmat_muiu)[s_mu,s_i,s_u]*
                #            cgsmat_mujv[s_mu,s_j,s_v]*
                #            cgsmat_aji[s_a,s_j,s_i]
                #cgsum = sum(cgo_contracted)*sum(cgs_contracted)::ComplexF64
                cgsum = compute_cgsum(
                            cgomat_avu , 
                            cgomat_muiu , 
                            cgomat_mujv ,   
                            cgomat_aji ,  
                            cgsmat_avu ,   
                            cgsmat_muiu ,   
                            cgsmat_mujv ,   
                            cgsmat_aji ,
                            cgo_contracted ,
                            cgs_contracted )

                # sign factor
                sign = ComplexF64((-1)^N_mu)

                # reduced matrix sector
                ijsector = @view ijredmat[(G_i,G_a,G_j)][r_i,:,r_j]

                # iterate over hopping multiplets m_a
                @inbounds for r_a::Int64 in 1:R_a

                    # reduced matrix element
                    redmatel = ijsector[r_a]::ComplexF64

                    thisredmat[r_u,r_a,r_v] = (sign/dim_u)*
                                              redmatel*
                                              cgsum::ComplexF64

                end
            end

            transformedmat::Array{ComplexF64,3} = zeros(ComplexF64,R_u,R_a,R_v)
            compute_transformedmat!( transformedmat , U_u , U_v , thisredmat )

            # insert result in final dict
            uvredmat[(G_u,G_a,G_v)] = transformedmat
        end
    end

    return uvredmat

end

@inline function compute_cgsum( 
    cgomat_avu::Array{ComplexF64,3} , 
    cgomat_muiu::Array{ComplexF64,3} , 
    cgomat_mujv::Array{ComplexF64,3} ,   
    cgomat_aji::Array{ComplexF64,3} ,  
    cgsmat_avu::Array{ComplexF64,3} ,   
    cgsmat_muiu::Array{ComplexF64,3} ,   
    cgsmat_mujv::Array{ComplexF64,3} ,   
    cgsmat_aji::Array{ComplexF64,3} ,
    cgo_contracted::Matrix{ComplexF64} ,
    cgs_contracted::Matrix{ComplexF64} )

    dI_a, dI_v,dI_u = size(cgomat_avu) 
    dI_mu,dI_i,_    = size(cgomat_muiu)
    _,    dI_j,dI_v = size(cgomat_mujv)
    dS_a, dS_v,dS_u = size(cgsmat_avu) 
    dS_mu,dS_i,_    = size(cgsmat_muiu)
    _,    dS_j,dS_v = size(cgsmat_mujv)

    cgo_contracted .= zero(ComplexF64)
    cgs_contracted .= zero(ComplexF64)

    @inbounds for i_u::Int64  in 1:dI_u,
                  i_ut::Int64 in 1:dI_u,
                  i_v::Int64  in 1:dI_v,
                  i_i::Int64  in 1:dI_i,
                  i_j::Int64  in 1:dI_j,
                  i_mu::Int64 in 1:dI_mu,
                  i_a::Int64  in 1:dI_a

        cgo_contracted[i_ut,i_u] += 
                conj(cgomat_avu[i_a,i_v,i_ut])*
                conj(cgomat_muiu[i_mu,i_i,i_u])*
                cgomat_mujv[i_mu,i_j,i_v]*
                cgomat_aji[i_a,i_j,i_i]::ComplexF64

    end
    @inbounds for s_u::Int64  in 1:dS_u,
                  s_ut::Int64 in 1:dS_u,
                  s_v::Int64  in 1:dS_v,
                  s_i::Int64  in 1:dS_i,
                  s_j::Int64  in 1:dS_j,
                  s_mu::Int64 in 1:dS_mu,
                  s_a::Int64  in 1:dS_a

        cgs_contracted[s_ut,s_u] += 
                conj(cgsmat_avu[s_a,s_v,s_ut])*
                conj(cgsmat_muiu[s_mu,s_i,s_u])*
                cgsmat_mujv[s_mu,s_j,s_v]*
                cgsmat_aji[s_a,s_j,s_i]::ComplexF64

    end

    return sum(cgo_contracted)*sum(cgs_contracted)::ComplexF64

end

function compute_transformedmat_mul!(
            redmat::Array{ComplexF64,3},
            U_u::Matrix{ComplexF64},
            U_v::Matrix{ComplexF64} )

    U_udagger = U_u'
    tmp = zeros( ComplexF64 , size(U_u,1) , size(U_v,1) )
    @views for j in 1:size(redmat,2)
        mul!( tmp , redmat[:,j,:] , U_v )
        mul!( redmat[:,j,:] , U_udagger , tmp )
    end

end
function compute_transformedmat!(
            transformedmat::Array{ComplexF64,3},
            U_u::Matrix{ComplexF64},
            U_v::Matrix{ComplexF64},
            thisredmat::Array{ComplexF64,3} )

    (R_u::Int64,R_a::Int64,R_v::Int64) = size(transformedmat)

    c::ComplexF64 = zero(ComplexF64)

    this::Matrix{ComplexF64}        = zeros(ComplexF64,R_u,R_v)
    transformed::Matrix{ComplexF64} = zeros(ComplexF64,R_u,R_v)
    @inbounds for r_a in 1:R_a 

        this        .= thisredmat[:,r_a,:]
        transformed .= zero(ComplexF64)

        @inbounds for r_v  in 1:R_v,
                      r_u  in 1:R_u

            c = zero(c)

            @inbounds for r_vp in 1:R_v,
                          r_up in 1:R_u

                c += conj(U_u[r_up,r_u]) * U_v[r_vp,r_v] * this[r_up,r_vp]

            end

            transformed[r_u,r_v] = c

        end

        transformedmat[:,r_a,:] .= transformed
    end
end

# Wrapper for the various spectral calculation methods
function compute_spectral_function(
            AA ,
            L ,
            iterations ,
            first_asymptotic_hopping_amplitude ;
            spectral_broadening::Float64=1.0 ,
            broadening_distribution::String="gaussian",
            method::String="sakai1989" ,
            label::String="" ,
            z::Float64=0.0 ,
            K_factor::Float64=2.0 ,
            orbitalresolved::Bool=false ,
            triple_excitation::Bool=false )


    # currently working
    if method=="sakai1989"
        compute_spectral_function_Sakai1989(
                AA ,
                L ,
                iterations ,
                first_asymptotic_hopping_amplitude ,
                spectral_broadening ,
                label ,
                z ,
                K_factor=K_factor ,
                orbitalresolved=orbitalresolved ,
                broadening_distribution=broadening_distribution ,
                triple_excitation=triple_excitation )
    # probably remove?
    elseif method=="fullrange"
        return compute_spectral_function_FullRange(
                AA ,
                L ,
                iterations ,
                first_asymptotic_hopping_amplitude ,
                spectral_broadening 
        )
    elseif method=="frota1986"
        return compute_spectral_function_Frota1986(
                AA ,
                L ,
                iterations ,
                first_asymptotic_hopping_amplitude
        )
    end
end

# Compute spectral functions for energies in the
# range between omega_N and maximum( maximum_energy , L )
function compute_spectral_function_FullRange( 
            AA ,
            L ,
            iterations ,
            alpha ,
            spectral_broadening ;
            z::Float64=0.0 )
    
    # combine energies in linear and logarithmic scales 
    omegas_log = Float64[ sign*L^(-(x-2)/2.0) for sign in [-1.0,1.0] 
                   for x in 2:0.1:(iterations-2)]
    omegas_lin = Float64[ sign*x for sign in [-1.0,1.0] 
                                   for x in 0.001:0.001:1]
    omegas::Vector{Float64} = sort(vcat( omegas_log , omegas_lin ))

    # spectral function
    spectral = Float64[0.0 for o in omegas]

    for (i,oo) in enumerate(omegas)

        o = -oo

        # odd iterations only in order to minimize oscillations
        for N in 3:2:iterations 

            # local A dictionary
            A = AA[N+1]

            # energy scale for this step
            omegaN  = Float64( alpha * L^(-(N-2)/2.0) )#* sqrt(L) )

            # maximum energy achieved in this step
            maximum_step_energy = maximum(collect( v[1] for v in values(A) ))

            # set reliable energy range 
            widthfac = maximum((L,maximum_step_energy))
            emin    = omegaN
            emax    = omegaN*widthfac

            # broadening factor for the step
            eta     = spectral_broadening*omegaN
            
            # negative interval 
            if is_in_interval( o , -emax , -emin )
                for (m,(e,coeffs)) in A
                    Δ = o+e*omegaN
                    spectral[i] += coeffs[1] *
                                   P(Δ,eta)
                end
            end
    
            #positive interval 
            if is_in_interval( o , emin , emax )
                for (m,(e,coeffs)) in A
                    Δ = (o-e*omegaN)
                    spectral[i] += coeffs[2]*
                                   P(Δ,eta)
                end
            end
        end
    end

    return [omegas spectral]

end

# Spectral function as computed by Sakai, Shimizu, Kasuya (1989).
# - Compute function only at one energy for each step, then
#   interpolate the results.
# - Not suitable for z averaging.
using Interpolations
function linear_interpolate_spectral_function(
            omegas_old::Vector{Float64},
            spectral_old::Vector{Float64},
            omegas_new::Vector{Float64} )

    # define interpolating function
    interpolator = linear_interpolation( omegas_old , spectral_old , extrapolation_bc=Line() )

    # interpolate to new values
    spectral_new = Float64[ interpolator(x) for x in omegas_new ]

    return spectral_new 
end
function spline_interpolate_spectral_function( 
        xx::Vector{Float64} , 
        yy::ST ;
        step_reductor::Int64=100 ,
        orbitalresolved::Bool=false ) where {ST<:Union{ Vector{Float64} , Dict{IntMultiplet,Vector{Float64}} }}

    # integer range for interpolation 
    number_of_energies = orbitalresolved ? length(yy[collect(keys(yy))[1]]) : length(yy)
    range_input = 1:number_of_energies
    range_output = 1:(1.0/(step_reductor-1)):number_of_energies

    # interpolate x and y
    interpolator_x = cubic_spline_interpolation( range_input , xx )
    xx_dense = map( interpolator_x , range_output )
    if !orbitalresolved
        interpolator_y = cubic_spline_interpolation( range_input , yy )
        yy_dense = map( interpolator_y , range_output )
    elseif orbitalresolved
        yy_dense = Dict()
        for (multiplet,spectral_function) in yy
            interpolator_y = cubic_spline_interpolation( range_input , spectral_function )
            push!(
                yy_dense ,
                multiplet => map( interpolator_y , range_output )
            )
        end
    end

    return xx_dense,yy_dense

end
function compute_spectral_function_Sakai1989( 
            AA ,
            L ,
            iterations ,
            first_asymptotic_hopping_amplitude::Float64 ,
            spectral_broadening::Float64 ,
            label::String ,
            z::Float64 ;
            broadening_distribution::String="gaussian",
            K_factor::Float64=2.0 ,
            orbitalresolved::Bool=false ,
            width_two::Bool=false ,
            triple_excitation::Bool=false )

    # even and odd iteration ranges
    even_iterator = 2:2:iterations
    odd_iterator  = 1:2:iterations
    # the chosen energies are 
    #
    #   e = K_factor * omega_N
    #
    omegas_odd = sort([ 
        ( K_factor * sign * first_asymptotic_hopping_amplitude * L^(-(N-2)/2.0) )
        for sign in [-1.0,1.0]
        for N in odd_iterator
    ])
    omegas_even = sort([ 
        ( K_factor * sign * first_asymptotic_hopping_amplitude * L^(-(N-2)/2.0) )
        for sign in [-1.0,1.0]
        for N in even_iterator
    ])
    scalings_odd = [ 
        ( first_asymptotic_hopping_amplitude * L^(-(N-2)/2.0) )
        for N in odd_iterator
    ]
    scalings_even = [ 
        ( first_asymptotic_hopping_amplitude * L^(-(N-2)/2.0) )
        for N in even_iterator
    ]

    # spectral function
    #
    # standard
    spectral_odd  = [0.0 for _ in omegas_odd]
    spectral_even = [0.0 for _ in omegas_even]
    # orbital-resolved
    if orbitalresolved
        excitation_multiplets = Set( multiplet for (_,multiplet) in keys(AA[1]) )
        spectral_odd  = Dict{IntMultiplet,Vector{Float64}}( 
            multiplet=>copy(spectral_odd) for multiplet in excitation_multiplets
        )
        spectral_even  = Dict{IntMultiplet,Vector{Float64}}( 
            multiplet=>copy(spectral_even) for multiplet in excitation_multiplets
        )
    end

    omega_positive_rescaled = +K_factor
    omega_negative_rescaled = -K_factor
    eta_rescaled = spectral_broadening
    for (i,N) in enumerate(odd_iterator)

        # energy at which to compute the spectral function
        omega_positive = -omegas_odd[i]
        omega_negative =  omegas_odd[i]

        # energy scale for step N
        omegaN = omega_positive/K_factor
        omegaN = scalings_odd[i]

        # local A dictionary
        A = AA[N+1]

        # rescaled broadening
        eta = spectral_broadening*omegaN

        if !orbitalresolved
            for (m,(eigenenergy,coeffs)) in A

                # positive energy range
                Delta_positive = omega_positive - eigenenergy*omegaN
                Delta_positive_rescaled = omega_positive_rescaled - eigenenergy
                #contribution = coeffs[1]*P(Delta_positive,eta)
                contribution = coeffs[1]*P(Delta_positive_rescaled,eta_rescaled)/omegaN
                isnan(contribution) && error( "nan 1: $contribution")
                if broadening_distribution=="loggaussian" && !iszero(coeffs[1])
                    contribution = coeffs[1]*P(Delta_positive_rescaled,eta_rescaled;distribution="loggaussian",E=eigenenergy,omega=omega_positive_rescaled)/omegaN
                    isnan(contribution) && error( "nan 2: $contribution")
                end
                spectral_odd[end-(i-1)] += contribution

                # negative energy range
                Delta_negative = omega_negative + eigenenergy*omegaN
                Delta_negative_rescaled = omega_negative_rescaled + eigenenergy
                #contribution = coeffs[2]*P(Delta_negative,eta)
                contribution = coeffs[2]*P(Delta_negative_rescaled,eta_rescaled)/omegaN
                isnan(contribution) && error( "nan 3: $contribution")
                if broadening_distribution=="loggaussian" && !iszero(coeffs[2])
                    contribution = coeffs[2]*P(Delta_negative_rescaled,eta_rescaled;distribution="loggaussian",E=eigenenergy,omega=omega_negative_rescaled)/omegaN
                    isnan(contribution) && error( "nan 4: $contribution")
                end
                spectral_odd[i] += contribution

            end
        elseif orbitalresolved
            for ((m,excitation_multiplet),(eigenenergy,coeffs)) in A

                # positive energy range
                Delta_positive = omega_positive - eigenenergy*omegaN
                Delta_positive_rescaled = omega_positive_rescaled - eigenenergy
                #contribution = coeffs[1]*P(Delta_positive,eta)
                contribution = coeffs[1]*P(Delta_positive_rescaled,eta_rescaled)/omegaN
                if broadening_distribution=="loggaussian" && !iszero(coeffs[1])
                    contribution = coeffs[1]*P(Delta_positive_rescaled,eta_rescaled;distribution="loggaussian",E=eigenenergy,omega=omega_positive_rescaled)/omegaN
                end
                spectral_odd[excitation_multiplet][end-(i-1)] += contribution

                # negative energy range
                Delta_negative = omega_negative + eigenenergy*omegaN
                Delta_negative_rescaled = omega_negative_rescaled + eigenenergy
                #contribution = coeffs[2]*P(Delta_negative,eta)
                contribution = coeffs[2]*P(Delta_negative_rescaled,eta_rescaled)/omegaN
                if broadening_distribution=="loggaussian" && !iszero(coeffs[2])
                    contribution = coeffs[2]*P(Delta_negative_rescaled,eta_rescaled;distribution="loggaussian",E=eigenenergy,omega=omega_negative_rescaled)/omegaN
                end
                spectral_odd[excitation_multiplet][i] += contribution

            end

        end

    end
    for (i,N) in enumerate(even_iterator)

        # energy at which to compute the spectral function
        omega_positive = -omegas_even[i]
        omega_negative =  omegas_even[i]

        # energy scale for step N
        omegaN = omega_positive/K_factor
        omegaN = scalings_even[i]

        # local A dictionary
        A = AA[N+1]

        # rescaled broadening
        eta = spectral_broadening*omegaN

        if !orbitalresolved
            for (_,(eigenenergy,coeffs)) in A

                # positive energy range
                Delta_positive = omega_positive - eigenenergy*omegaN
                Delta_positive_rescaled = omega_positive_rescaled - eigenenergy
                #contribution = coeffs[1]*P(Delta_positive,eta)
                contribution = coeffs[1]*P(Delta_positive_rescaled,eta_rescaled)/omegaN
                isnan(contribution) && error( "nan 5: $contribution")
                if broadening_distribution=="loggaussian" && !iszero(coeffs[1])
                    contribution = coeffs[1]*P(Delta_positive_rescaled,eta_rescaled;distribution="loggaussian",E=eigenenergy,omega=omega_positive_rescaled)/omegaN
                    isnan(contribution) && error( "nan 6: $contribution")
                end
                spectral_even[end-(i-1)] += contribution

                # negative energy range
                Delta_negative = omega_negative + eigenenergy*omegaN
                Delta_negative_rescaled = omega_negative_rescaled + eigenenergy
                #contribution = coeffs[2]*P(Delta_negative,eta)
                contribution = coeffs[2]*P(Delta_negative_rescaled,eta_rescaled)/omegaN
                isnan(contribution) && error( "nan 7: $contribution")
                if broadening_distribution=="loggaussian" && !iszero(coeffs[2])
                    contribution = coeffs[2]*P(Delta_negative_rescaled,eta_rescaled;distribution="loggaussian",E=eigenenergy,omega=omega_negative_rescaled)/omegaN
                    isnan(contribution) && error( "nan 8: $contribution $eigenenergy")
                end
                spectral_even[i] += contribution

            end
        elseif orbitalresolved
            for ((_,excitation_multiplet),(eigenenergy,coeffs)) in A

                # positive energy range
                Delta_positive = omega_positive - eigenenergy*omegaN
                Delta_positive_rescaled = omega_positive_rescaled - eigenenergy
                #contribution = coeffs[1]*P(Delta_positive,eta)
                contribution = coeffs[1]*P(Delta_positive_rescaled,eta_rescaled)/omegaN
                if broadening_distribution=="loggaussian" && !iszero(coeffs[1])
                    contribution = coeffs[1]*P(Delta_positive_rescaled,eta_rescaled;distribution="loggaussian",E=eigenenergy,omega=omega_positive_rescaled)/omegaN
                end
                spectral_even[excitation_multiplet][end-(i-1)] += contribution

                # negative energy range
                Delta_negative = omega_negative + eigenenergy*omegaN
                Delta_negative_rescaled = omega_negative_rescaled + eigenenergy
                #contribution = coeffs[2]*P(Delta_negative,eta)
                contribution = coeffs[2]*P(Delta_negative_rescaled,eta_rescaled)/omegaN
                if broadening_distribution=="loggaussian" && !iszero(coeffs[2])
                    contribution = coeffs[2]*P(Delta_negative_rescaled,eta_rescaled;distribution="loggaussian",E=eigenenergy,omega=omega_negative_rescaled)/omegaN
                end
                spectral_even[excitation_multiplet][i] += contribution

            end
        end

    end

    # average even and odd calculations with linear interpolation
    #
    # collect energy (omega) values
    omegas_evenodd = sort(vcat(omegas_even,omegas_odd))
    if width_two
        filter!( x->(abs(x)<=1.0) , omegas_evenodd )
        omegas_evenodd[1]!==-1.0  && insert!( omegas_evenodd , 1 , -1.0 )
        omegas_evenodd[end]!==1.0 && push!( omegas_evenodd , 1.0 )
    else
        max_energy = minimum((maximum(abs.(omegas_even)),maximum(abs.(omegas_odd))))
        min_energy = maximum((minimum(abs.(omegas_even)),minimum(abs.(omegas_odd))))
        filter!( x->abs(x)<=max_energy , omegas_evenodd )
        filter!( x->abs(x)>=min_energy , omegas_evenodd )
        omegas_evenodd[1]>-max_energy && insert!( omegas_evenodd , 1 , -max_energy )
        omegas_evenodd[end]<max_energy && push!( omegas_evenodd , max_energy )
    end
    # interpolate even and odd spectral functions to new omegas
    if !orbitalresolved
        spectral_even_interpolated = linear_interpolate_spectral_function( 
            omegas_even  , 
            spectral_even , 
            omegas_evenodd 
        )
        spectral_odd_interpolated  = linear_interpolate_spectral_function( 
            omegas_odd , 
            spectral_odd , 
            omegas_evenodd 
        )
    elseif orbitalresolved
        spectral_even_interpolated = Dict(
            excitation_multiplet => linear_interpolate_spectral_function( 
                omegas_even , 
                spectral_even_orbital , 
                omegas_evenodd 
            )
            for (excitation_multiplet,spectral_even_orbital) in spectral_even
        )
        spectral_odd_interpolated = Dict(
            excitation_multiplet => linear_interpolate_spectral_function( 
                omegas_odd  , 
                spectral_odd_orbital , 
                omegas_evenodd 
            )
            for (excitation_multiplet,spectral_odd_orbital) in spectral_odd
        )
    end
    # compute average 
    if !orbitalresolved
        spectral_evenodd = 0.5*( spectral_even_interpolated + spectral_odd_interpolated )
    elseif orbitalresolved
        spectral_evenodd = Dict(
            excitation_multiplet => 0.5*( spectral_even_interpolated[excitation_multiplet] + 
                                          spectral_odd_interpolated[excitation_multiplet]   )
            for excitation_multiplet in excitation_multiplets
        )
    end

    # compute spline-interpolated smooth spectral function
    omegas_evenodd_spline,spectral_evenodd_spline = spline_interpolate_spectral_function( 
        omegas_evenodd ,
        spectral_evenodd ,
        orbitalresolved=orbitalresolved
   )

    # write even, odd, average and spline-interpolated spectral data
    if !orbitalresolved
        if triple_excitation
            write_spectral_function( 
                spectral_filename(label,z=z,tail="_triple_even") ,
                [omegas_even spectral_even]
            )
            write_spectral_function(
                spectral_filename(label,z=z,tail="_triple_odd") ,
                [omegas_odd spectral_odd]
            )
            write_spectral_function(
                spectral_filename(label,z=z,tail="_triple") ,
                [omegas_evenodd spectral_evenodd]
            )
            write_spectral_function(
                spectral_filename(label,z=z,tail="_triple_splined") ,
                [omegas_evenodd_spline spectral_evenodd_spline]
            )
        else
            write_spectral_function( 
                spectral_filename(label,z=z,tail="_even") ,
                [omegas_even spectral_even]
            )
            write_spectral_function(
                spectral_filename(label,z=z,tail="_odd") ,
                [omegas_odd spectral_odd]
            )
            write_spectral_function(
                spectral_filename(label,z=z) ,
                [omegas_evenodd spectral_evenodd]
            )
            write_spectral_function(
                spectral_filename(label,z=z,tail="_splined") ,
                [omegas_evenodd_spline spectral_evenodd_spline]
            )
        end
    elseif orbitalresolved
        for (orbital_idx,multiplet) in enumerate(collect(excitation_multiplets))
            orbital_header = "# Excitation multiplet: $multiplet , index=$orbital_idx\n"
            write_spectral_function( 
                spectral_filename(label,z=z,tail="_even",orbital=orbital_idx) ,
                [omegas_even spectral_even[multiplet]] ,
                orbitalresolved_header=orbital_header
            )
            write_spectral_function(
                spectral_filename(label,z=z,tail="_odd",orbital=orbital_idx) ,
                [omegas_odd spectral_odd[multiplet]] ,
                orbitalresolved_header=orbital_header
            )
            write_spectral_function(
                spectral_filename(label,z=z,orbital=orbital_idx) ,
                [omegas_evenodd spectral_evenodd[multiplet]] ,
                orbitalresolved_header=orbital_header
            )
            write_spectral_function(
                spectral_filename(label,z=z,tail="_splined",orbital=orbital_idx) ,
                [omegas_evenodd_spline spectral_evenodd_spline[multiplet]] ,
                orbitalresolved_header=orbital_header
            )
        end
    end

end

# Discrete spectral function in Frota & Oliveira 1986
function compute_spectral_function_Frota1986( 
            AA ,
            L ,
            iterations ,
            first_asymptotic_hopping_amplitude )

    spectral_vector = Matrix{Float64}[]
    for N in 3:2:iterations 

        # spectral weights
        A = AA[N+1]

        # energy scale
        omegaN  = Float64( first_asymptotic_hopping_amplitude * L^(-(N-2)/2.0) )
        
        # introduce spectral weights
        for (m,(e,coeffs)) in A

            # energy
            omega = e * omegaN

            # weights in the positive and negative regions
            isapprox( coeffs[2] , 0.0 ) || push!( spectral_vector , [ omega coeffs[2]] )
            isapprox( coeffs[1] , 0.0 ) || push!( spectral_vector , [-omega coeffs[1]] )
        end

    end

    # sort vector by energies
    sort!( spectral_vector , by=x->x[1,1] )

    # concatenate to obtain matrix
    spectral_matrix = vcat(spectral_vector...)

    return spectral_matrix

end

function compute_spectral_function_orbitalresolved(
             AA ,
             L ,
             iterations ,
             alpha ,
             spectral_broadening ;
             widthfac::Float64=2.0 )

     omegas_log = [ sign*L^(-(x-2)/2.0) for sign in [-1.0,1.0]
                                    for x in 2:0.1:(iterations-10)]
     omegas_lin = [ sign*x for sign in [-1,1] for x in 0.001:0.001:1]
     omegas = vcat( omegas_log , omegas_lin )
     sort!( omegas )

     multiplets_a = Set([ k[2] for k in keys(AA[1]) ])
     spectral = Dict( m_a=> [0.0 for o in omegas] for m_a in multiplets_a )

     for (i,oo) in enumerate(omegas)

         o = -oo

         for N in 3:2:iterations

             omegaN = Float64(alpha * L^(-(N-2)/2.0) )
             emin    = omegaN/sqrt(L)
             emax    = omegaN*sqrt(L) #widthfac
             local A = AA[N+1]
             eta     = spectral_broadening*omegaN

             # negative interval
             if is_in_interval( o , -emax , -emin )
                 for ((m,m_a),(e,coeffs)) in A
                     D = o+e*omegaN
                     spectral[m_a][i] += coeffs[1] * P(D,eta)
                 end
             end

             #positive interval
             if is_in_interval( o , emin , emax )
                 for ((m,m_a),(e,coeffs)) in A
                     D = (o-e*omegaN)
                     spectral[m_a][i] += coeffs[2]* P(D,eta)
                 end
             end
         end
     end

     return Dict( m_a => [omegas spectralo.*(2/integ(omegas,spectralo))] for (m_a,spectralo) in spectral )

 end

function integ( X , Y ) 

    integrala = 0.0

    for (i,x) in enumerate(X[1:(end-1)])

        (a,b) = (x,X[i+1]) 

        Δ = (b-a) 
        V = (Y[i]+Y[i+1])/2.0

        integrala += Δ*V

    end

    return integrala 

end

function setup_redmat_AA( 
            matdict ,
            multiplets_atom ,
            multiplets_operator ,
            cg_o_fullmatint ,
            cg_s_fullmatint ,
            irrEU ,
            oindex2dimensions ;
            verbose=false )

    # reduced matrix
    Mred = get_redmat3( matdict ,
                        multiplets_atom ,
                        multiplets_operator ,
                        cg_o_fullmatint ,
                        cg_s_fullmatint ;
                        verbose=verbose)
    if verbose
        println( "M reduced" )
        print_dict( Mred )
        println()
    end

    # spectral thermo 
    GG0 = [G for (G,(E,e)) in irrEU for e in E if iszero(e)]
    part0 = sum(map( G->oindex2dimensions[G[2]]*(G[3]+1) , GG0 ))
    if verbose 
        @show part0
        println()
    end

    # normalize irrEU just in case
    irrEU = normalize_irrEU( irrEU )

    # spectral info
    A = redM2A( Mred,
                collect(multiplets_operator),
                cg_o_fullmatint,
                cg_s_fullmatint,
                irrEU,
                part0;
                verbose=verbose)
    if verbose
        @show A 
        println()
    end
    AA = [A]

    return (Mred,AA) 

end
function setup_selfenergy( 
            Mred ,
            matdict ,
            multiplets_atom ,
            multiplets_operator ,
            cg_o_fullmatint ,
            cg_s_fullmatint ,
            irrEU ,
            oindex2dimensions ;
            verbose=false )

    # reduced matrix
    Mred_se = get_redmat3( matdict ,
                           multiplets_atom ,
                           multiplets_operator ,
                           cg_o_fullmatint ,
                           cg_s_fullmatint ;
                           verbose=verbose)

    if verbose
        println( "M reduced" )
        print_dict( Mred )
        println()
    end

    # spectral thermo 
    GG0 = [G for (G,(E,e)) in irrEU for e in E if iszero(e)]
    part0 = sum(map( G->oindex2dimensions[G[2]]*(G[3]+1) , GG0 ))
    if verbose 
        @show part0
        println()
    end

    # normalize irrEU just in case
    irrEU = normalize_irrEU( irrEU )

    # spectral info
    A_se = redM2A_selfenergy( 
                Mred,
                Mred_se ,
                collect(multiplets_operator),
                cg_o_fullmatint,
                cg_s_fullmatint,
                irrEU,
                part0;
                verbose=verbose)
    print_dict(A_se)

    if verbose
        @show A_se 
        println()
    end
    AA_se = [A_se]

    return (Mred_se,AA_se) 

end
function setup_redmat_AA_orbitalresolved(
            matdict ,
            multiplets_atom ,
            multiplets_operator ,
            cg_o_fullmatint ,
            cg_s_fullmatint ,
            irrEU ,
            oindex2dimensions ;
            verbose=false )

    # reduced matrix
    Mred = get_redmat3( matdict ,
                        multiplets_atom ,
                        multiplets_operator ,
                        cg_o_fullmatint ,
                        cg_s_fullmatint ;
                        verbose=verbose)
    if verbose
        println( "M reduced" )
        print_dict( Mred )
        println()
    end

    # spectral thermo
    G0 = [G for (G,(E,U)) in irrEU if E[1]==0][1]
    I0,S0 = G0[2:3]
    D0s = S0+1
    D0o = oindex2dimensions[I0]
    part0 = D0o*D0s
    if verbose
        @show part0
        println()
    end

    # normalize irrEU just in case
    irrEU = normalize_irrEU( irrEU )

    # spectral info
    A = redM2A_orbitalresolved( 
                Mred,
                collect(multiplets_operator),
                cg_o_fullmatint,
                cg_s_fullmatint,
                irrEU,
                part0;
                verbose=verbose
    )
    if verbose
        @show A
        println()
    end
    AA = [A]

    return (Mred,AA)

end

function get_partition0( irrEU , 
                         oindex2dimensions ;
                         verbose=false )

    G0 = [G for (G,(E,U)) in irrEU if E[1]==0][1]
    I0,S0 = G0[2:3]
    D0s = S0+1
    D0o = oindex2dimensions[I0]
    part0 = D0o*D0s
    if verbose 
        @show part0
        println()
    end

    return part0 
end

function update_redmat_AA(
            Mred ,
            irrEU ,
            combinations_uprima ,
            multiplets_a ,
            cg_o_fullmatint ,
            cg_s_fullmatint ,
            AA ,
            oindex2dimensions ;
            verbose=false )

    Mred = get_new_blockredmat( 
                Mred , 
                irrEU ,
                combinations_uprima ,
                collect(multiplets_a) ,
                cg_o_fullmatint ,
                cg_s_fullmatint ,
                oindex2dimensions )

    # spectral thermo 
    G0 = [G for (G,(E,U)) in irrEU if E[1]==0][1]
    I0,S0 = G0[2:3]
    D0s = S0+1
    D0o = oindex2dimensions[I0]
    part0 = D0o*D0s

    push!(AA,
         redM2A(Mred,
                collect(multiplets_a),
                cg_o_fullmatint,
                cg_s_fullmatint,
                irrEU,
                part0;
                verbose=verbose) 
         )
    return ( Mred , AA )

end
function update_redmat_AA_CGsummethod(
            Mred::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } ,
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            combinations_uprima::Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} },
            multiplets_a::Vector{NTuple{4,Int64}} ,
            cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            Karray_orbital::Array{ComplexF64,6} , 
            Karray_spin::Array{ComplexF64,6} , 
            AA ,
            oindex2dimensions::Vector{Int64} ;
            verbose=false ,
            density_matrix::Dict{IntIrrep,Matrix{Float64}}=Dict{IntIrrep,Matrix{Float64}}() ,
            multiplets_kept::Set{IntMultiplet}=Set{IntMultiplet}() )

    Mred = update_blockredmat( 
                multiplets_a,
                Karray_orbital ,
                Karray_spin ,
                Mred,
                combinations_uprima ,
                irrEU )

    if verbose 
        println( "M" )
        print_dict( Mred ) 
    end

    # spectral thermo 
    G0 = [G for (G,(E,U)) in irrEU if E[1]==0][1]
    I0,S0 = G0[2:3]
    D0s = S0+1
    D0o = oindex2dimensions[I0]
    part0 = D0o*D0s

    push!(
        AA,
        redM2A( 
            Mred,
            multiplets_a,
            cg_o_fullmatint,
            cg_s_fullmatint,
            irrEU,
            part0;
            verbose=false,
            density_matrix=density_matrix,
            multiplets_kept=multiplets_kept
        ) 
    )

    if verbose 
        println( "AA" )
        print_dict( AA )
    end

    return ( Mred , AA )

end
function update_selfenergy_CGsummethod(
            Mred::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } ,
            Mred_se::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } ,
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            combinations_uprima::Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} },
            multiplets_a::Vector{NTuple{4,Int64}} ,
            cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            Karray_orbital::Array{ComplexF64,6} , 
            Karray_spin::Array{ComplexF64,6} , 
            AA_se::Vector{Dict{NTuple{4,Int64},Tuple{Float64,Vector{Float64}}}} ,
            oindex2dimensions::Vector{Int64} ;
            verbose=false )

    Mred_se = update_blockredmat( 
                multiplets_a,
                Karray_orbital ,
                Karray_spin ,
                Mred_se,
                combinations_uprima ,
                irrEU )

    if verbose 
        println( "M_se" )
        print_dict( Mred_se ) 
    end

    # spectral thermo 
    G0 = [G for (G,(E,U)) in irrEU if E[1]==0][1]
    I0,S0 = G0[2:3]
    D0s = S0+1
    D0o = oindex2dimensions[I0]
    part0 = D0o*D0s

    push!(
        AA_se,
        redM2A_selfenergy( 
            Mred,
            Mred_se,
            multiplets_a,
            cg_o_fullmatint,
            cg_s_fullmatint,
            irrEU,
            part0;
            verbose=false
        ) 
    )

    if verbose 
        println( "AA_se" )
        print_dict( AA_se )
    end

    return ( Mred_se , AA_se )

end
function update_redmat_AA_CGsummethod_orbitalresolved(
            Mred::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } ,
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            combinations_uprima::Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} },
            multiplets_a::Vector{NTuple{4,Int64}} ,
            cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            Karray_orbital::Array{ComplexF64,6} , 
            Karray_spin::Array{ComplexF64,6} , 
            AA::Vector{Dict{NTuple{2,IntMultiplet},Tuple{Float64,Vector{Float64}}}} ,
            oindex2dimensions::Vector{Int64} ;
            verbose=false )

    Mred = update_blockredmat( 
                multiplets_a,
                Karray_orbital ,
                Karray_spin ,
                Mred,
                combinations_uprima ,
                irrEU )
    #Mred = get_new_blockredmat_CGsummethod3( 
    #            Mred,
    #            collect(multiplets_a),
    #            combinations_uprima ,
    #            irrEU ,
    #            Karray_orbital ,
    #            Karray_spin ,
    #            cg_o_fullmatint ,
    #            cg_s_fullmatint )


    if verbose 
        println( "M" )
        print_dict( Mred ) 
    end

    # spectral thermo 
    G0 = [G for (G,(E,U)) in irrEU if E[1]==0][1]
    I0,S0 = G0[2:3]
    D0s = S0+1
    D0o = oindex2dimensions[I0]
    part0 = D0o*D0s

    push!(
        AA,
        redM2A_orbitalresolved(
                Mred,
                multiplets_a,
                cg_o_fullmatint,
                cg_s_fullmatint,
                irrEU,
                part0;
                verbose=false) 
    )

    if verbose 
        println( "AA" )
        print_dict( AA )
    end

    return ( Mred , AA )

end


function get_new_blockredmat_CGsummethod( 
                iaj::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } ,
                multiplets_a::Vector{NTuple{4,Int64}} ,
                combinations_uprima::Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} },
                irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
                Karray_orbital::Array{ComplexF64,6} ,
                Karray_spin::Array{ComplexF64,6} ,
                cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
                cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} })

    combs_uvprima::Dict{ NTuple{2,NTuple{3,Int64}} , Vector{NTuple{6,NTuple{4,Int64}}} } = Dict()
    mm_uv::Vector{NTuple{6,NTuple{4,Int64}}} = []
    @inbounds for (G_u,mm_u) in combinations_uprima,
                  (G_v,mm_v) in combinations_uprima

        mm_uv = [(m_u,m_mu,m_i,m_v,m_nu,m_j) 
                 for (m_u,m_mu,m_i) in mm_u
                 for (m_v,m_nu,m_j) in mm_v 
                 if (m_mu==m_nu && m_i[1]==(m_j[1]+1))]

        combs_uvprima[(G_u,G_v)] = mm_uv 

    end

    iajsizes = Set( size(m) for m in values(iaj) )
    max_i = maximum([ s[1] for s in iajsizes ])
    max_a = maximum([ s[2] for s in iajsizes ])
    max_j = maximum([ s[3] for s in iajsizes ])
    zeromat = zeros( ComplexF64 , max_i , max_a , max_j )

    GG_a::Vector{NTuple{3,Int64}} = collect(Set( m[1:3] for m in multiplets_a ))
    G2R_a::Dict{ NTuple{3,Int64} , Int64 } = Dict( G=>length(Set(m for m in multiplets_a if m[1:3]==G)) 
                                                   for G in GG_a )
    G2R_uv::Dict{ NTuple{3,Int64} , Int64 } = Dict( 
                G=>length(ucombs) 
                for (G,ucombs) in combinations_uprima )
    GG_uv::Vector{NTuple{3,Int64}} = collect(keys(combinations_uprima))


    #uav::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } = Dict(
    #        (G_u,G_a,G_v)=>zeros(ComplexF64,R_u,R_a,R_v)
    #        for (G_u,R_u) in G2R_uv for (G_v,R_v) in G2R_uv for (G_a,R_a) in G2R_a
    #       )
    uav::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } = Dict()

    # iterate over irreps of u, a, v
    @inbounds for (G_u::NTuple{3,Int64},R_u::Int64) in G2R_uv,
                  (G_v::NTuple{3,Int64},R_v::Int64) in G2R_uv,
                  (G_a::NTuple{3,Int64},R_a::Int64) in G2R_a

        # irrep quantum numbers
        (N_u::Int64,I_u::Int64,S_u::Int64) = G_u
        (N_a::Int64,I_a::Int64,S_a::Int64) = G_a
        (N_v::Int64,I_v::Int64,S_v::Int64) = G_v
        N_u==(N_v+1) || continue
        haskey( cg_o_fullmatint , (I_a,I_v,I_u) ) || continue 
        haskey( cg_s_fullmatint , (S_a,S_v,S_u) ) || continue
        
        # diagonalization matrices
        U_u::Matrix{ComplexF64} = irrEU[G_u][2]
        U_v::Matrix{ComplexF64} = irrEU[G_v][2]

        # G_u, G_a, G_v block
        uav_block::Array{ComplexF64,3} = zeros(ComplexF64,R_u,R_a,R_v)

        # iterate over block matrix elements
        @inbounds for r_u::Int64 in 1:R_u,
                      r_v::Int64 in 1:R_v,
                      r_a::Int64 in 1:R_a
            
            # matrix element for r_u, r_a, r_v
            uav_matel::ComplexF64 = zero(ComplexF64)

            # iterate over multiplet combinations in G_u,G_v where m_mu=m_nu
            @inbounds for (m_up::NTuple{4,Int64},
                           m_munu::NTuple{4,Int64},
                           m_i::NTuple{4,Int64},
                           m_vp::NTuple{4,Int64},
                           m_nu::NTuple{4,Int64},
                           m_j::NTuple{4,Int64}) in combs_uvprima[(G_u,G_v)]

                # quantum numbers
                r_up = m_up[4]
                r_vp = m_vp[4]
                G_i::NTuple{3,Int64}, r_i = m_i[1:3],m_i[4]
                G_j::NTuple{3,Int64}, r_j = m_j[1:3],m_j[4]
                G_munu,r_munu = m_munu[1:3],m_munu[4]
                (N_i,   I_i,   S_i)    = G_i
                (N_munu,I_munu,S_munu) = G_munu
                (N_j,   I_j,   S_j)    = G_j

                # reduced matrix element from previous step
                #haskey( iaj , (G_i,G_a,G_j) ) || continue
                iaj_matel::ComplexF64 = get( iaj , (G_i,G_a,G_j) , zeromat )[r_i,r_a,r_j]
                iaj_matel==zero(iaj_matel) && continue
                # diagonalization matrix element
                uu::ComplexF64 = conj(U_u[r_up,r_u])*U_v[r_vp,r_v]
                # sign factor
                sign::ComplexF64 = ComplexF64((-1)^N_munu)
                # clebsch-gordan sum
                K::ComplexF64 = Karray_orbital[I_u,I_v,I_a,I_munu,I_i,I_j]*
                                Karray_spin[((S_u,S_v,S_a,S_munu,S_i,S_j).+1)...]

                # insert everything 
                uav_matel += uu*sign*iaj_matel*K

            end

            # assign matrix element to block
            uav_block[r_u,r_a,r_v] = uav_matel

        end

        #isapprox( sum(abs.(uav_block)) , zero(ComplexF64) ) && continue
        uav[G_u,G_a,G_v] = uav_block

    end

    return uav

end
# iterate over <i||a||j> instead of <u||a||v>
function get_new_blockredmat_CGsummethod2( 
                iaj::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } ,
                multiplets_a::Vector{NTuple{4,Int64}} ,
                combinations_uprima::Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} },
                irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
                Karray_orbital::Array{ComplexF64,6} ,
                Karray_spin::Array{ComplexF64,6} ,
                cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
                cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} })

    GG_a::Vector{NTuple{3,Int64}} = collect(Set( m[1:3] for m in multiplets_a ))
    G2R_a::Dict{ NTuple{3,Int64} , Int64 } = Dict( G=>length(Set(m for m in multiplets_a if m[1:3]==G)) 
                                                   for G in GG_a )
    G2R_uv::Dict{ NTuple{3,Int64} , Int64 } = Dict( 
                G=>length(ucombs) 
                for (G,ucombs) in combinations_uprima )
    GG_uv::Vector{NTuple{3,Int64}} = collect(keys(combinations_uprima))

    uav::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } = Dict(
            (G_u,G_a,G_v)=>zeros(ComplexF64,R_u,R_a,R_v)
            for (G_u,R_u) in G2R_uv for (G_v,R_v) in G2R_uv for (G_a,R_a) in G2R_a
           )
    uav_diag::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } = Dict(
            (G_u,G_a,G_v)=>zeros(ComplexF64,R_u,R_a,R_v)
            for (G_u,R_u) in G2R_uv for (G_v,R_v) in G2R_uv for (G_a,R_a) in G2R_a
           )

    combinations_i::Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} } = Dict()
    for (G_u,combs) in combinations_uprima,
        (m_u,m_mu,m_i) in combs

        G_i = m_i[1:3]
        if haskey( combinations_i , G_i ) 
            push!( combinations_i[G_i] , (m_u,m_mu,m_i) )
        else 
            combinations_i[G_i] = [(m_u,m_mu,m_i)]
        end

    end

    #uav::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } = Dict()
    @inbounds for ((G_i,G_a,G_j),iaj_mat) in iaj 

        # irrep quantum numbers
        (N_i,I_i,S_i) = G_i
        (N_a,I_a,S_a) = G_a
        (N_j,I_j,S_j) = G_j

        # discard cut-off multiplets
        (haskey(combinations_i,G_i) && haskey(combinations_i,G_j)) || continue
        
        # a dimension 
        R_a = size(iaj_mat)[2]

        # iterate over u' and v' combinations
        # generated from i and j
        @inbounds for (m_u,m_mu,m_i) in combinations_i[G_i],
                      (m_v,m_nu,m_j) in combinations_i[G_j]

            # quantum numbers
            (N_u, I_u, S_u, r_u ) = m_u
            (N_v, I_v, S_v, r_v ) = m_v
            (N_mu,I_mu,S_mu,r_mu) = m_mu
            G_u = (N_u,I_u,S_u)
            G_v = (N_v,I_v,S_v)
            r_i = m_i[4]
            r_j = m_j[4]
    
            # early discard
            m_mu==m_nu || continue
            
            # clebsch-gordan K sum
            K::ComplexF64 = Karray_orbital[I_u,I_v,I_a,I_mu,I_i,I_j]*
                            Karray_spin[((S_u,S_v,S_a,S_mu,S_i,S_j).+1)...]
            K==zero(K) && continue

            # sign factor
            sign::ComplexF64 = ComplexF64((-1)^N_mu)
            
            uav[(G_u,G_a,G_v)][r_u,:,r_v] = sign*K*iaj_mat[r_i,:,r_j]

        end

    end
    
    
    count = 0
    @inbounds for ((G_u,G_a,G_v),diagmat) in uav_diag 

        @show size(diagmat)
        count += 1


        U_u = irrEU[G_u][2]
        U_v = irrEU[G_v][2]

        compute_transformedmat!( mat , U_u , U_v , uav[G_u,G_a,G_v] )

        nondiagmat = uav[(G_u,G_a,G_v)]
        #@einsum diagmat[r_u,r_a,r_v] = conj(U_u[r_up,r_u])*U_v[r_vp,r_v]*nondiagmat[r_up,r_a,r_vp]

        #@inbounds for r_u in 1:R_u,
        #              r_a in 1:R_a,
        #              r_v in 1:R_v 

        #    matel::ComplexF64 = zero(ComplexF64)

        #    @inbounds for r_up in 1:R_u,
        #                  r_vp in 1:R_v

        #        matel += conj(U_u[r_up,r_u])*U_v[r_vp,r_v]*uavmat[r_up,r_a,r_vp]

        #    end

        #    mat[r_u,r_a,r_v] = matel

        #end

    end
    @show count
    println()

    return uav_diag

end
# same method as get_new_blockredmat() 
# but with K arrays instead of computing cgsum
function get_new_blockredmat_CGsummethod3( 
            ijredmat::Dict{NTuple{3,NTuple{3,Int64}},Array{ComplexF64,3}} ,
            multiplets_a::Vector{NTuple{4,Int64}} ,
            combinations_uprima::Dict{NTuple{3,Int64}, Vector{NTuple{3,NTuple{4,Int64}}}},
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            Karray_orbital::Array{ComplexF64,6} ,
            Karray_spin::Array{ComplexF64,6},
            cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} })

    uvredmat::Dict{NTuple{3,NTuple{3,Int64}},Array{ComplexF64,3}} = Dict()

    # zeromat
    iajsizes = Set( size(m) for m in values(ijredmat) )
    length(iajsizes)==0 && (return uvredmat)
    max_i = maximum([ s[1] for s in iajsizes ])
    max_a = maximum([ s[2] for s in iajsizes ])
    max_j = maximum([ s[3] for s in iajsizes ])
    zeromat = zeros( ComplexF64 , max_i , max_a , max_j )

    # hopping irreps and multiplicities
    irrmult_a = get_irreps(Set(multiplets_a);multiplicity=true)
    
    # iterate over irrep combinations G_u and G_v
    # in the new step
    for (G_u::NTuple{3,Int64},ucombs::Vector{NTuple{3,NTuple{4,Int64}}}) in combinations_uprima,
        (G_v::NTuple{3,Int64},vcombs::Vector{NTuple{3,NTuple{4,Int64}}}) in combinations_uprima

        # irrep quantum numbers
        (N_u,I_u,S_u) = G_u 
        (N_v,I_v,S_v) = G_v 

        # early discard
        N_u==(N_v+1) || continue
        
        # diagonalization matrices 
        U_u = irrEU[G_u][2]::Matrix{ComplexF64}  
        U_v=  irrEU[G_v][2]::Matrix{ComplexF64} 

        # irrep multiplicities
        R_u = length(ucombs)
        R_v = length(vcombs)

        # iterate over hopping irreps G_a
        for (G_a::NTuple{3,Int64},R_a::Int64) in irrmult_a

            # irrep quantum numbers
            (N_a,I_a,S_a) = G_a

            # early discard
            haskey( cg_o_fullmatint , (I_a,I_v,I_u) ) || continue
            haskey( cg_s_fullmatint , (S_a,S_v,S_u) ) || continue

            # initiate new matrix entry
            thisredmat     = zeros(ComplexF64,R_u,R_a,R_v)
            transformedmat = zeros(ComplexF64,R_u,R_a,R_v)

            # iterate over primed combinations, to be 
            # combined later by u matrices.
            # here:
            #   u = u'
            #   v = v'
            for (m_u::NTuple{4,Int64},m_mu::NTuple{4,Int64},m_i::NTuple{4,Int64}) in ucombs,
                (m_v::NTuple{4,Int64},m_nu::NTuple{4,Int64},m_j::NTuple{4,Int64}) in vcombs

                # for block matrix elements, same shell states
                m_mu==m_nu || continue

                # multiplet quantum numbers
                r_u = m_u[4]
                r_v = m_v[4]
                (N_mu,I_mu,S_mu,r_mu) = m_mu
                (N_i, I_i, S_i, r_i ) = m_i 
                (N_j, I_j, S_j, r_j ) = m_j
                G_mu = (N_mu,I_mu,S_mu)
                G_i = (N_i,I_i,S_i)
                G_j = (N_j,I_j,S_j)

                ((G_i,G_a,G_j) in keys(ijredmat)) || continue

                # clebsch-gordan sum coefficient
                @inbounds Ko = Karray_orbital[I_u,I_v,I_a,I_mu,I_i,I_j]
                @inbounds Ks = Karray_spin[((S_u,S_v,S_a,S_mu,S_i,S_j).+1)...]
                K::ComplexF64 = Ko * Ks
                K==zero(K) && continue
                #isnan(K) && error( "NaN in K" )

                # sign factor
                sign = ComplexF64((-1)^N_mu)

                # reduced matrix sector
                ijsector = @view ijredmat[(G_i,G_a,G_j)][r_i,:,r_j]
                #any(isnan.(ijsector)) && error( "NaN in ijsector" )
                @inbounds @. thisredmat[r_u,:,r_v] = sign*K*ijsector[:]
            end

            if any(isnan.(thisredmat)) 
                error( "NaN in impurity excitation matrix. This is probably because the calculation requires a value of the optional parameter max_spin2 higher than the one given as input (2S=10 is the maximum by default). The maximum value of 2S required is printed in the output." )
            end

            transformedmat::Array{ComplexF64,3} = zeros(ComplexF64,R_u,R_a,R_v)
            #compute_transformedmat!( transformedmat , U_u , U_v , thisredmat )
            compute_transformedmat_mul!( thisredmat , U_u , U_v )

            # insert result in final dict
            uvredmat[(G_u,G_a,G_v)] = thisredmat
        end
    end

    return uvredmat

end

function print_A( A ; orbitalresolved::Bool=false )

    if orbitalresolved
        excitation_multiplets = Set( m_a for (m,m_a) in keys(A) )

        @printf "%-20s %-18s %-15s %-10s\n" "impurity excitation" "excited multiplet" "energy" "amplitudes"
        for excitation_multiplet in excitation_multiplets
            for ((m,m_a),(eigenenergy,coefficients)) in A
                m_a==excitation_multiplet || continue
                @printf "%-20s %-18s %-15.3f %-10.3e\n" m_a m eigenenergy coefficients[1]
                @printf "%-20s %-18s %-15s %-10.3e\n" "" "" "" coefficients[2]
            end
        end
    end

end

# optimized version of M updating
function get_new_blockredmat_old( 
            ijredmat::Dict{NTuple{3,NTuple{3,Int64}},Array{ComplexF64,3}} ,
            multiplets_a::Vector{NTuple{4,Int64}} ,
            combinations_uprima::Dict{NTuple{3,Int64}, Vector{NTuple{3,NTuple{4,Int64}}}},
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            Karray_orbital::Array{ComplexF64,6} ,
            Karray_spin::Array{ComplexF64,6},
            cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} },
            cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} })

    # initialize result matrix
    uvredmat::Dict{NTuple{3,NTuple{3,Int64}},Array{ComplexF64,3}} = Dict()

    # hopping irreps and multiplicities
    irrmult_a = get_irreps(Set(multiplets_a);multiplicity=true)
    
    # iterate over irrep combinations G_u and G_v
    # in the new step
    for (G_u::NTuple{3,Int64},ucombs::Vector{NTuple{3,NTuple{4,Int64}}}) in combinations_uprima,
        (G_v::NTuple{3,Int64},vcombs::Vector{NTuple{3,NTuple{4,Int64}}}) in combinations_uprima

        # irrep quantum numbers
        (N_u,I_u,S_u) = G_u 
        (N_v,I_v,S_v) = G_v 

        # early discard
        N_u==(N_v+1) || continue
        
        # diagonalization matrices 
        U_u = irrEU[G_u][2]::Matrix{ComplexF64}  
        U_v=  irrEU[G_v][2]::Matrix{ComplexF64} 

        # irrep multiplicities
        R_u = length(ucombs)
        R_v = length(vcombs)

        # iterate over hopping irreps G_a
        for (G_a::NTuple{3,Int64},R_a::Int64) in irrmult_a

            # irrep quantum numbers
            (N_a,I_a,S_a) = G_a

            # early discard
            haskey( cg_o_fullmatint , (I_a,I_v,I_u) ) || continue
            haskey( cg_s_fullmatint , (S_a,S_v,S_u) ) || continue

            # initiate new matrix entry
            thisredmat     = zeros(ComplexF64,R_u,R_a,R_v)
            transformedmat = zeros(ComplexF64,R_u,R_a,R_v)

            # iterate over primed combinations, to be 
            # combined later by u matrices.
            # here:
            #   u = u'
            #   v = v'
            for (m_u::NTuple{4,Int64},m_mu::NTuple{4,Int64},m_i::NTuple{4,Int64}) in ucombs,
                (m_v::NTuple{4,Int64},m_nu::NTuple{4,Int64},m_j::NTuple{4,Int64}) in vcombs

                # for block matrix elements, same shell states
                m_mu==m_nu || continue

                # multiplet quantum numbers
                r_u = m_u[4]
                r_v = m_v[4]
                (N_mu,I_mu,S_mu,r_mu) = m_mu
                (N_i, I_i, S_i, r_i ) = m_i 
                (N_j, I_j, S_j, r_j ) = m_j
                G_mu = (N_mu,I_mu,S_mu)
                G_i = (N_i,I_i,S_i)
                G_j = (N_j,I_j,S_j)

                ((G_i,G_a,G_j) in keys(ijredmat)) || continue

                # clebsch-gordan sum coefficient
                @inbounds Ko = Karray_orbital[I_u,I_v,I_a,I_mu,I_i,I_j]
                @inbounds Ks = Karray_spin[((S_u,S_v,S_a,S_mu,S_i,S_j).+1)...]
                K::ComplexF64 = Ko * Ks
                K==zero(K) && continue
                #isnan(K) && error( "NaN in K" )

                # sign factor
                sign = ComplexF64((-1)^N_mu)

                # reduced matrix sector
                ijsector = @view ijredmat[(G_i,G_a,G_j)][r_i,:,r_j]
                #any(isnan.(ijsector)) && error( "NaN in ijsector" )
                @inbounds @. thisredmat[r_u,:,r_v] = sign*K*ijsector[:]
            end

            if any(isnan.(thisredmat)) 
                error( "NaN in impurity excitation matrix. This is probably because the calculation requires a value of the optional parameter max_spin2 higher than the one given as input (2S=10 is the maximum by default). The maximum value of 2S required is printed in the output." )
            end

            transformedmat::Array{ComplexF64,3} = zeros(ComplexF64,R_u,R_a,R_v)
            #compute_transformedmat!( transformedmat , U_u , U_v , thisredmat )
            compute_transformedmat_mul!( thisredmat , U_u , U_v )

            # insert result in final dict
            uvredmat[(G_u,G_a,G_v)] = thisredmat
        end
    end

    return uvredmat

end

function update_blockredmat(
                multiplets_a_block::Vector{NTuple{4,Int64}}, 
                Ksum_o_array::Array{ComplexF64,6} ,
                Ksum_s_array::Array{ComplexF64,6} ,
                redmat_iaj::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} },
                combinations_uprima::Dict{ IntIrrep , Vector{NTuple{3,IntMultiplet}}} ,
                irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } )

    # G => multiplicity of G
    #
    # i (n) / u (n-1)
    G2R_uv::Dict{IntIrrep,Int64} = Dict( G_u=>size(U_u,1) for (G_u,(_,U_u)) in irrEU )
    # a
    G2R_a::Dict{NTuple{3,Int64},Int64} = Dict( 
        G=>R for (G,R) in get_irreps( Set(multiplets_a_block) ; multiplicity=true ) 
    )
    
    # irrep -> decomposition -> multiplet
    Gu2GiGmu2uimumults::Dict{IntIrrep,Dict{NTuple{2,IntIrrep},Vector{NTuple{3,Int64}}}} = Dict(
        G_u => Dict(
            (G_i,G_mu) => [(m_u[4],m_i[4],m_mu[4]) for (m_u,m_mu,m_i) in multiplet_combinations_Gu if (m_i[1:3]==G_i && m_mu[1:3]==G_mu)]
            for (G_i,G_mu) in Set( (m_i[1:3],m_mu[1:3]) for (_,m_mu,m_i) in combinations_uprima[G_u])
        )
        for (G_u,multiplet_combinations_Gu) in combinations_uprima
    )

    # full-sized matrices for avoiding allocations
    R_uv_max = maximum(values(G2R_uv))
    R_a_max = maximum(values(G2R_a))
    tmp_full = zeros(ComplexF64,R_uv_max,R_uv_max)
    uav_matrix_full = zeros(ComplexF64,R_uv_max,R_a_max,R_uv_max)

    # initialize result dictionary
    redmat_uav::Dict{NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} } = Dict()

    # G_u, G_v iteration
    for (G_u,GiGmu2uimumults) in Gu2GiGmu2uimumults,
        (G_v,GjGnu2vjnumults) in Gu2GiGmu2uimumults

        # irrep quantum numbers 
        N_u,I_u,S_u = G_u
        N_v,I_v,S_v = G_v

        # multiplicities
        R_u= G2R_uv[G_u]
        R_v= G2R_uv[G_v]

        # transformation matrices
        @views begin
            U_u::Matrix{ComplexF64} = irrEU[G_u][2]
            U_v::Matrix{ComplexF64} = irrEU[G_v][2]
        end
        # temporary matrix for transformation
        @views tmp = tmp_full[1:R_u,1:R_v]
        tmp .= zero(ComplexF64)

        # G_a iteration
        for (G_a,R_a) in G2R_a

            # irrep quantum numbers
            N_a,I_a,S_a = G_a

            # early discard 
            N_u==(N_v+N_a) || continue

            # < G_u || f^\dagger_{G_a} || G_v >
            @views uav_matrix = uav_matrix_full[1:R_u,1:R_a,1:R_v]
            uav_matrix .= zero(ComplexF64)

            # G_i,G_mu,G_j,G_nu iteration
            for ((G_i,G_mu),uimumults) in GiGmu2uimumults,
                ((G_j,G_nu),vjnumults) in GjGnu2vjnumults

                # irrep quantum numbers 
                N_i,I_i,S_i = G_i
                N_j,I_j,S_j = G_j
                N_munu,I_munu,S_munu = G_mu

                # permutation factor
                sign = (-1)^N_munu

                # early discard
                G_mu==G_nu     || continue
                N_i==(N_j+N_a) || continue

                # clebsch-gordan sum
                @inbounds K = Ksum_o_array[I_u,I_v,I_a,I_munu,I_i,I_j]*
                              Ksum_s_array[((S_u,S_v,S_a,S_munu,S_i,S_j).+1)...]
                K==zero(K) && continue

                # total prefactor
                sign_times_K = sign * K

                # < G_mu || f^\dagger_{G_a} || G_j >
                haskey( redmat_iaj , (G_i,G_a,G_j) ) || continue
                @views iaj_matrix = redmat_iaj[(G_i,G_a,G_j)]

                # u,i,mu , v,j,nu  multiplet iteration
                for (r_u,r_i,r_mu) in uimumults,
                    (r_v,r_j,r_nu) in vjnumults

                    r_mu==r_nu || continue

                    @inbounds for r_a in 1:R_a 
                        uav_matrix[r_u,r_a,r_v] += sign_times_K * iaj_matrix[r_i,r_a,r_j]
                    end

                end # u,i,mu , v,j,nu , a multiplet iteration

            end # G_i,G_mu,G_j,G_nu iteration

            # final discard
            is_matrix_zero(uav_matrix) && continue

            # transform matrix
            @views for r_a in 1:G2R_a[G_a]
                mul!( tmp , uav_matrix[:,r_a,:] , U_v )
                mul!( uav_matrix[:,r_a,:] , U_u' , tmp )
            end

            # redmat_uav
            push!( redmat_uav , (G_u,G_a,G_v) => uav_matrix[:,:,:] )
            #@inbounds redmat_uav[(G_u,G_a,G_v)] = uav_matrix_full[1:R_u,1:R_a,1:R_v]
            #@inbounds redmat_uav[(G_u,G_a,G_v)] = uav_matrix[:,:,:]

        end # G_a iteration
    end # G_u, G_v iteration

    return redmat_uav
end

function is_matrix_zero( matrix::SubArray{ComplexF64,3} )
    s::Float64 = zero(Float64)
    @inbounds for i in eachindex(matrix)
        s += abs2(matrix[i])
        s > zero(s) && return false
    end
    return true
end

# ------ #
# DM-NRG #
# ------ #
function backwards_dm_reduced(
        dm::Dict{IntIrrep,Matrix{Float64}} ,
        diagonalizers::Dict{IntIrrep,Matrix{ComplexF64}} ,
        combinations_uprima::Dict{ IntIrrep , Vector{NTuple{3,IntMultiplet}} } ,
        oindex2dimensions::Vector{Int64} ,
        multiplets_kept::Set{IntMultiplet} )

    combinations_uprima_modified = Dict(
        G_u => comb 
        for (G_u,comb) in combinations_uprima
        if haskey(dm,G_u)
    )

    GG_i = Set( 
        m_i[1:3] 
        for (G_u,combs_u) in combinations_uprima
        for (_,_,m_i) in combs_u
    )
    G2mm_i = Dict( 
        G_i=>Set(
            m_i 
            for (G_u,combs_u) in combinations_uprima for 
            (_,_,m_i) in combs_u
            if m_i[1:3]==G_i
        )
        for G_i in GG_i
    )
    G2multiplicity_i = Dict( G_i=>length(mm_i) for (G_i,mm_i) in G2mm_i )

    dm_back::Dict{IntIrrep,Matrix{Float64}} = Dict(
        G_i=>zeros(Float64,multiplicity_G_i,multiplicity_G_i)
        for (G_i,multiplicity_G_i) in G2multiplicity_i
    )

    for (G_u,combinations_G_u) in combinations_uprima_modified

        # diagonalization matrices and DM from next step
        du = dm[G_u]
        u = diagonalizers[G_u]

        _,I_u,S_u = G_u

        for (m_u,m_mu,m_i) in combinations_G_u,
            (m_v,m_nu,m_j) in combinations_G_u

            G_ij = m_i[1:3]

            # trace over mu and symmetry of the DM
            m_mu==m_nu     || continue
            G_ij==m_j[1:3] || continue

            # quantum numbers
            r_u = m_u[4]
            r_v = m_v[4]
            _,I_ij,S_ij,r_i = m_i
            r_j = m_j[4]

            dim_i = oindex2dimensions[I_ij]*(S_ij+1)
            dim_u = oindex2dimensions[I_u]*(S_u+1)
            cgsum = dim_u/dim_i

            contrib = dot(u[r_u,1:size(du,1)],du,conj.(u[r_v,1:size(du,1)]))

            dm_back[G_ij][r_i,r_j] += real(cgsum*contrib)

        end

    end
    println()

    dm_back = Dict( G=>m for (G,m) in dm_back if !isapprox(norm(m),0.0;atol=1e-6) )

    return dm_back
end
