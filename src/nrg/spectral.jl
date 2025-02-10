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
            tail::String="" ,
            spectraldir::String="spectral"
        )

    z_string = zavg ? "_zavg" : "_z$(z)"
    o_string = orbital==0 ? "" : "_o$(orbital)"
    return "$(spectraldir)/spectral_$(label)$(z_string)$(o_string)$(tail).dat"

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
        write( f , header )
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
@inline function C( spectral_energy::Float64 ,
                    excitation_energy::Float64 ,
                    broadening_factor::Float64 ,
                    distribution::String ;
                    tq::Float64=0.0)

    if distribution=="gaussian"

        return exp(-((spectral_energy-excitation_energy)/broadening_factor)^2) / (broadening_factor*sqrt(pi))

    elseif distribution=="loggaussian" && !iszero(excitation_energy)

        return exp(-broadening_factor^2/4.0) /
               (broadening_factor*abs(excitation_energy)*sqrt(pi)) *
               exp(-(log(abs(spectral_energy/excitation_energy))/broadening_factor)^2)

    elseif distribution=="lorentzian"

        return broadening_factor/(2*pi*((spectral_energy-excitation_energy)^2+broadening_factor^2))

    elseif distribution=="interpolated"

        return (broadening_factor*sqrt(pi))^-1 *
               exp(-(interpolator(spectral_energy,tq)-interpolator(excitation_energy,tq))^2/broadening_factor^2) *
               interpolatorder(excitation_energy,tq)

    end
    error( "Could not compute peak convolution." )
end
@inline function interpolator(ω::Float64,tq::Float64)
    frac = ω/tq
    return 0.5 * tanh(frac) * log(frac^2+exp(tq))
end
@inline function interpolatorder(ω::Float64,tq::Float64)
    frac = ω/tq 
    e = exp(tq)
    return (2*tq)^-1 * sech(frac)^2 * log((frac)^2+e) +
           frac/tq * tanh(frac) * (frac^2+e)^-1
end
@inline function P( w_pm_E::Float64 , 
                    eta::Float64 ; 
                    distribution::String="gaussian" , 
                    omega::Float64=0.0 , 
                    E::Float64=0.0 )

    if distribution=="gaussian" 
        return P_Gaussian( w_pm_E , eta ) 
    elseif distribution=="loggaussian" 
        return P_LogGaussian( w_pm_E , eta , omega , E )
    elseif distribution=="lorentzian"
        return P_Lorentzian( w_pm_E , eta )
    end
end
@inline function P_Gaussian( w_pm_E::Float64 , eta::Float64 )
    return exp(-(w_pm_E/eta)^2) / (eta*sqrt(pi))
end
@inline function P_LogGaussian( w_pm_E::Float64 , eta::Float64 , omega::Float64 , E::Float64 )
    iszero(E) && (return 0.0)
    return exp(-eta^2/4.0)/(eta*abs(E)*sqrt(pi)) * exp(-(log(abs(omega/E))/eta)^2)
end
@inline P_Lorentzian( w_pm_E::Float64 , eta::Float64 ) = eta/(2*pi*(w_pm_E^2+eta^2))


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
function is_matrix_zero( matrix::SubArray{ComplexF64,4} )
    s::Float64 = zero(Float64)
    @inbounds for i in eachindex(matrix)
        s += abs2(matrix[i])
        s > zero(s) && return false
    end
    return true
end
function is_matrix_zero( matrix::Array{ComplexF64,4} )
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
function backwards_dm_reduced_T0(
            dm::Dict{IntIrrep,Matrix{Float64}} ,
            diagonalizers::Dict{IntIrrep,Matrix{ComplexF64}} ,
            combinations_uprima::Dict{ IntIrrep , Vector{NTuple{3,IntMultiplet}} } ,
            oindex2dimensions::Vector{Int64} )

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
            for (G_u,combs_u) in combinations_uprima 
            for (_,_,m_i) in combs_u
            if m_i[1:3]==G_i
        )
        for G_i in GG_i
    )
    G2multiplicity_i = Dict( G_i=>length(mm_i) for (G_i,mm_i) in G2mm_i )

    # density matrix for the kept states |i> entering |u> as |μ>⊗|i>->|u>
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

            @views G_ij = m_i[1:3]

            # trace over mu and symmetry of the DM
            m_mu==m_nu     || continue
            @views G_ij==m_j[1:3] || continue

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

    dm_back = Dict( G=>m for (G,m) in dm_back if !iszero(norm(m)) )

    return dm_back
end
function backwards_dm_reduced_Tnonzero(
        dm::Dict{IntIrrep,Matrix{Float64}} ,
        diagonalizers::Dict{IntIrrep,Matrix{ComplexF64}} ,
        combinations_uprima::Dict{ IntIrrep , Vector{NTuple{3,IntMultiplet}} } ,
        oindex2dimensions::Vector{Int64} ,
        multiplets_kept::Set{IntMultiplet} ,
        multiplets_disc::Set{IntMultiplet} ,
        T::Float64 ,
        Z::Float64 ,
        shell_dimension::Int64 ,
        iterations::Int64 ,
        iteration::Int64 ,
        iteration_scale::Float64 ,
        spectrum::Dict{IntIrrep,Vector{Float64}} )

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

    GG_i_disc = Set( m_i[1:3] for m_i in multiplets_disc )
    G2rr_i_disc = Dict(
        G=>[m[end] for m in filter(m->m[1:3]==G,multiplets_disc)]
        for G in GG_i_disc
    )

    mm_i_all = union(multiplets_kept,multiplets_disc)
    GG_i_all = Set( m[1:3] for m in mm_i_all )
    G2mm_i_all = Dict( G=>filter(m->m[1:3]==G,mm_i_all) for G in GG_i_all )
    G2R_i_all = Dict( G=>length(mm) for (G,mm) in G2mm_i_all )

    dm_back::Dict{IntIrrep,Matrix{Float64}} = Dict(
        G_i=>zeros(Float64,R_i,R_i)
        for (G_i,R_i) in G2R_i_all
    )

    # contribution from next steps
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

    # diagonal contribution from discarded states
    for ((N,I,S),rr_i) in G2rr_i_disc

        @views dm_b = dm_back[N,I,S]
        @views multiplet_spectrum = spectrum[N,I,S]

        for r_i in rr_i

            dm_b[r_i,r_i] = shell_dimension^(iterations-iteration)*
                            exp(-multiplet_spectrum[r_i]*iteration_scale/T)/
                            Z

        end
    end

    return dm_back
end

# ----------------------------- #
# General correlation functions #
# ----------------------------- #

# peaks for correlation function 
#       < [ A , B ] >,
# with A given as A^† for multiplets 
# A and B transforming according to the
# same irreducible representations
function compute_correlation_peaks( 
        A::Dict{ IntTripleG , Array{ComplexF64,3} } ,
        B::Dict{ IntTripleG , Array{ComplexF64,3} } ,
        oindex2dimensions::Vector{Int64} ,
        irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
        z::Float64 ,
        iteration::Int64 ;
        correlation_type::String="spectral" ,
        T::Float64=0.0 ,
        iteration_scale::Float64=1.0 ,
        use_density_matrix::Bool=false )

    if (correlation_type=="spectral" && A==B)
        compute_spectral_peaks(A,oindex2dimensions,irrEU,z,iteration;T=T,iteration_scale=iteration_scale,use_density_matrix=use_density_matrix)
    end

end
function compute_spectral_peaks(
        A::Dict{ IntTripleG , Array{ComplexF64,3} } ,
        oindex2dimensions::Vector{Int64} ,
        irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
        z::Float64 ,
        iteration::Int64;
        T::Float64=0.0 ,
        iteration_scale::Float64=1.0 ,
        use_density_matrix::Bool=false )

    if !use_density_matrix

        if iszero(T)
            compute_spectral_peaks_T0(A,oindex2dimensions,irrEU,z,iteration)
        else
            compute_spectral_peaks_Tnonzero(A,oindex2dimensions,irrEU,z,iteration,T,iteration_scale)
        end

    elseif use_density_matrix

        compute_spectral_peaks_density_matrix_T0(A,oindex2dimensions,irrEU,z,iteration)

    end
end
function compute_spectral_peaks_T0(
        A::Dict{ IntTripleG , Array{ComplexF64,3} } ,
        oindex2dimensions::Vector{Int64} ,
        irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
        z::Float64 ,
        iteration::Int64 )

    # multiplets of the A operator
    multiplets_a::Set{IntMultiplet} = Set( (G_a...,r_a) 
                                            for ((_,G_a,_),redmat) in A
                                            for r_a in 1:size(redmat,2) )

    # peak information
    peaks::Dict{ IntMultiplet , Vector{Vector{Float64}} } = Dict( 
        m_a=>Vector{Float64}[] for m_a in multiplets_a
    )

    # ground multiplets and partition function
    G2R_0::Dict{IntIrrep,Int64} = Dict( 
        G => length(filter(iszero,E))
        for (G,(E,U)) in irrEU 
    )
    partition = sum( R*oindex2dimensions[I]*(S+1) for ((_,I,S),R) in G2R_0 )

    # iterate over irreducible matrices
    for ((G_u,G_a,G_v),redmat) in A 

        # irrep quantum numbers
        (_,I_u,S_u) = G_u 
        (_,I_a,S_a) = G_a

        @views energies_Gu = irrEU[G_u][1]
        @views energies_Gv = irrEU[G_v][1]

        # clebsch-gordan contribution 
        dim_u = oindex2dimensions[I_u]*(S_u+1)
        dim_a = oindex2dimensions[I_a]*(S_a+1)
        cg = dim_u/dim_a

        # iterate over multiplets in irrep combination
        for rrr in CartesianIndices(redmat)

            # outer multiplicities
            (r_u,r_a,r_v) = Tuple(rrr)

            # multiplets
            m_a = (G_a...,r_a)

            # coefficient from reduced matrix element
            redmatel = abs2(redmat[rrr])

            e_u = energies_Gu[r_u]
            e_v = energies_Gv[r_v]

            # total contribution 
            w = cg*redmatel

            # particle (positive) contribution: m_v ground, m_u not ground
            if (iszero(e_v) && !iszero(e_u)) || (iszero(e_u) && !iszero(e_v))

                push!( peaks[m_a] , [ e_u-e_v , w/partition ] )

            end
        end
    end

    # write to file
    z_dir = "spectral/peaks/particle/z$(z)"
    iteration_dir = "spectral/peaks/particle/z$(z)/n$iteration"
    mkdir(iteration_dir)
    for (m_a,m_a_peaks) in peaks

        N_a,I_a,S_a,r_a = m_a
        multiplet_file = "$(iteration_dir)/$(N_a)_$(I_a)_$(S_a)_$(r_a)"

        sort!( m_a_peaks , by=x->x[1] )
        peaks_matrix = reduce(vcat,transpose.(m_a_peaks))
        writedlm(multiplet_file,peaks_matrix)

    end
end
function compute_spectral_peaks_Tnonzero(
        A::Dict{ IntTripleG , Array{ComplexF64,3} } ,
        oindex2dimensions::Vector{Int64} ,
        irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
        z::Float64 ,
        iteration::Int64 ,
        T::Float64 ,
        iteration_scale::Float64 )

    # multiplets of the A operator
    multiplets_a::Set{IntMultiplet} = Set( (G_a...,r_a) 
                                            for ((_,G_a,_),redmat) in A
                                            for r_a in 1:size(redmat,2) )

    # peak information
    peaks::Dict{ IntMultiplet , Vector{Vector{Float64}} } = Dict( 
        m_a=>Vector{Float64}[] for m_a in multiplets_a
    )

    # partition function
    partition = sum( oindex2dimensions[I]*(S+1)*exp(-e*iteration_scale/T)
                     for ((_,I,S),(E,_)) in irrEU
                     for e in E )

    # avoid extra computation
    scale_over_T = iteration_scale/T

    # iterate over irreducible matrices
    for ((G_u,G_a,G_v),redmat) in A 

        # irrep quantum numbers
        (_,I_u,S_u) = G_u 
        (_,I_a,S_a) = G_a

        # energy vectors
        @views energies_Gu = irrEU[G_u][1]
        @views energies_Gv = irrEU[G_v][1]

        # clebsch-gordan contribution 
        dim_u = oindex2dimensions[I_u]*(S_u+1)
        dim_a = oindex2dimensions[I_a]*(S_a+1)
        cg = dim_u/dim_a

        # iterate over multiplets in irrep combination
        for r_a in axes(redmat,2)

            peaks_m_a = peaks[G_a...,r_a]

            for r_u in axes(redmat,1),
                r_v in axes(redmat,3)

                # coefficient from reduced matrix element
                @inbounds redmatel = abs2(redmat[r_u,r_a,r_v])

                # multiplet energies
                @inbounds e_u = energies_Gu[r_u]
                @inbounds e_v = energies_Gv[r_v]

                # total contribution 
                boltzmann = exp(-e_u*scale_over_T)+exp(-e_v*scale_over_T)
                w = cg*redmatel*boltzmann

                push!( peaks_m_a , [ e_u-e_v , w/partition ] )
            end
        end
    end

    # write to file
    z_dir = "spectral/peaks/particle/z$(z)"
    iteration_dir = "spectral/peaks/particle/z$(z)/n$iteration"
    mkdir(iteration_dir)
    for (m_a,m_a_peaks) in peaks

        N_a,I_a,S_a,r_a = m_a
        multiplet_file = "$(iteration_dir)/$(N_a)_$(I_a)_$(S_a)_$(r_a)"

        sort!( m_a_peaks , by=x->x[1] )
        peaks_matrix = reduce(vcat,transpose.(m_a_peaks))
        writedlm(multiplet_file,peaks_matrix)

    end
end
function compute_spectral_peaks_density_matrix_T0(
        A::Dict{ IntTripleG , Array{ComplexF64,3} } ,
        oindex2dimensions::Vector{Int64} ,
        irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
        z::Float64 ,
        iteration::Int64 )

    # read density matrix
    density_matrix_dir = "spectral/densitymatrices/n$iteration"
    density_matrix = Dict{IntIrrep,Matrix{Float64}}()
    for irrepfile in readdir(density_matrix_dir)
        G = (map(x->parse(Int64,x),split(irrepfile,"_"))...,)
        density_matrix[G] = readdlm("$(density_matrix_dir)/$(irrepfile)")
    end

    # multiplets of the A operator
    multiplets_a::Set{IntMultiplet} = Set( (G_a...,r_a) 
                                            for ((_,G_a,_),redmat) in A
                                            for r_a in 1:size(redmat,2) )

    # peak information
    peaks::Dict{ IntMultiplet , Vector{Vector{Float64}} } = Dict( 
        m_a=>Vector{Float64}[] for m_a in multiplets_a
    )

    # iterate over irreducible matrices
    for ((G_u,G_a,G_v),redmat) in A 

        # irrep quantum numbers
        (_,I_u,S_u) = G_u 
        (_,I_a,S_a) = G_a

        @views energies_Gu = irrEU[G_u][1]
        @views energies_Gv = irrEU[G_v][1]

        u_ground = haskey(density_matrix,G_u)
        v_ground = haskey(density_matrix,G_v)
        (u_ground || v_ground) || continue
        u_ground && (@views dm_u = density_matrix[G_u])
        v_ground && (@views dm_v = density_matrix[G_v])

        # clebsch-gordan contribution 
        dim_u = oindex2dimensions[I_u]*(S_u+1)
        dim_a = oindex2dimensions[I_a]*(S_a+1)
        cg = dim_u/dim_a

        # iterate over multiplets in irrep combination
        for r_a in axes(redmat,2)

            m_a = (G_a...,r_a)

            # particle sum over ground states
            if v_ground
                @inbounds for r_v1 in axes(dm_v,3),
                              r_v2 in axes(dm_v,3)

                    # density matrix element
                    dmel = dm_v[r_v1,r_v2]
                    isapprox(dmel,zero(dmel)) && continue

                    # energies of the kept ("ground") states
                    e_v1 = energies_Gv[r_v1]
                    e_v2 = energies_Gv[r_v2]

                    # iterate over excited (discarded) states
                    for r_u in axes(redmat,1)

                        redmatel = conj(redmat[r_u,r_a,r_v1])*redmat[r_u,r_a,r_v2]
                        isapprox(redmatel,zero(redmatel)) && continue

                        e_u = energies_Gu[r_u]

                        # total weight 
                        w = cg*redmatel*dmel

                        push!( peaks[m_a] , [ e_u-0.5*(e_v1+e_v2) , w ] )

                    end
                end
            end

            # hole sum over ground states
            if u_ground
                @inbounds for r_u1 in axes(dm_u,1),
                              r_u2 in axes(dm_u,1)

                    # density matrix element
                    dmel = dm_u[r_u1,r_u2]
                    isapprox(dmel,zero(dmel)) && continue

                    # energies of the kept ("ground") states
                    e_u1 = energies_Gu[r_u1]
                    e_u2 = energies_Gu[r_u2]

                    # iterate over excited (discarded) states
                    for r_v in axes(redmat,3)

                        redmatel = redmat[r_u1,r_a,r_v]*conj(redmat[r_u2,r_a,r_v])
                        isapprox(redmatel,zero(redmatel)) && continue

                        e_v = energies_Gv[r_v]

                        # total weight 
                        w = cg*redmatel*dmel

                        push!( peaks[m_a] , [ 0.5*(e_u1+e_u2)-e_v , w ] )

                    end
                end
            end

        end
    end

    # write to file
    z_dir = "spectral/peaks/particle/z$(z)"
    iteration_dir = "spectral/peaks/particle/z$(z)/n$iteration"
    mkdir(iteration_dir)
    for (m_a,m_a_peaks) in peaks

        N_a,I_a,S_a,r_a = m_a
        multiplet_file = "$(iteration_dir)/$(N_a)_$(I_a)_$(S_a)_$(r_a)"

        sort!( m_a_peaks , by=x->x[1] )
        peaks_matrix = reduce(vcat,transpose.(m_a_peaks))
        writedlm(multiplet_file,peaks_matrix)

    end
end
function add_correlation_contribution!( 
            correlation_dict::Dict{IntMultiplet,Matrix{Float64}} ,
            A::Dict{ IntTripleG , Array{ComplexF64,3} } ,
            B::Dict{ IntTripleG , Array{ComplexF64,3} } ,
            oindex2dimensions::Vector{Int64} ,
            irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
            iteration::Int64 ,
            broadening_distribution::String ,
            spectral_broadening::Float64 ,
            iteration_scale::Float64 ,
            K_factor::Float64 ;
            correlation_type::String="spectral" ,
            T::Float64=0.0 ,
            limit_shell::Bool=false ,
            extra_iterations::Int64=0 ,
            density_matrix::Dict{IntIrrep,Matrix{Float64}}=Dict{IntIrrep,Matrix{Float64}}() ,
            multiplets_kept::Vector{IntMultiplet}=IntMultiplet[] ,
            multiplets_discarded::Vector{IntMultiplet}=IntMultiplet[] ,
            L::Float64=0.0 ,
            scale::Float64=0.0 ,
            half_weight_idx::Int64=0 ,
            half_weight_energy::Float64=0.0 )

    if (correlation_type=="spectral" && A==B)
        add_spectral_contribution!(correlation_dict,
                                   A,
                                   oindex2dimensions,
                                   irrEU,
                                   iteration,
                                   broadening_distribution,
                                   spectral_broadening,
                                   iteration_scale,
                                   K_factor;
                                   T=T,
                                   limit_shell=limit_shell,
                                   extra_iterations=extra_iterations,
                                   density_matrix=density_matrix,
                                   multiplets_kept=multiplets_kept,
                                   multiplets_discarded=multiplets_discarded,
                                   L=L,
                                   scale=scale,
                                   half_weight_idx=half_weight_idx,
                                   half_weight_energy=half_weight_energy)
    end

end
function add_spectral_contribution!(
            correlation_dict::Dict{IntMultiplet,Matrix{Float64}} ,
            A::Dict{ IntTripleG , Array{ComplexF64,3} } ,
            oindex2dimensions::Vector{Int64} ,
            irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
            iteration::Int64,
            broadening_distribution::String ,
            spectral_broadening::Float64 ,
            iteration_scale::Float64 ,
            K_factor::Float64 ;
            T::Float64=0.0 ,
            limit_shell::Bool=false ,
            extra_iterations::Int64=0 ,
            density_matrix::Dict{IntIrrep,Matrix{Float64}}=Dict{IntIrrep,Matrix{Float64}}() ,
            multiplets_kept::Vector{IntMultiplet}=IntMultiplet[] ,
            multiplets_discarded::Vector{IntMultiplet}=IntMultiplet[] ,
            L::Float64=0.0 ,
            scale::Float64=0.0 ,
            half_weight_idx::Int64=0 ,
            half_weight_energy::Float64=half_weight_energy )

    use_density_matrix = !iszero(length(density_matrix))

    if !use_density_matrix

        if iszero(T)
            add_spectral_contribution_T0!(correlation_dict,
                                          A,
                                          oindex2dimensions,
                                          irrEU,
                                          iteration,
                                          broadening_distribution,
                                          spectral_broadening,
                                          iteration_scale,
                                          K_factor)
        else
            add_spectral_contribution_Tnonzero!(correlation_dict,
                                                A,
                                                oindex2dimensions,
                                                irrEU,
                                                iteration,
                                                broadening_distribution,
                                                spectral_broadening,
                                                iteration_scale,
                                                K_factor,
                                                T,
                                                limit_shell,
                                                extra_iterations)
        end

    elseif use_density_matrix

        if iszero(T)
            add_spectral_contribution_density_matrix_T0!(correlation_dict,
                                                         A,
                                                         oindex2dimensions,
                                                         irrEU,
                                                         iteration,
                                                         broadening_distribution,
                                                         spectral_broadening,
                                                         iteration_scale,
                                                         K_factor,
                                                         density_matrix,
                                                         multiplets_kept,
                                                         multiplets_discarded,
                                                         L,
                                                         scale)
        else
            add_spectral_contribution_density_matrix_Tnonzero!(correlation_dict,
                                                               A,
                                                               oindex2dimensions,
                                                               irrEU,
                                                               iteration,
                                                               broadening_distribution,
                                                               spectral_broadening,
                                                               iteration_scale,
                                                               K_factor,
                                                               density_matrix,
                                                               multiplets_kept,
                                                               multiplets_discarded,
                                                               L,
                                                               scale,
                                                               half_weight_idx,
                                                               half_weight_energy,
                                                               T)

        end

    end
end
function add_spectral_contribution_T0!(
            correlation_dict::Dict{IntMultiplet,Matrix{Float64}} ,
            A::Dict{ IntTripleG , Array{ComplexF64,3} } ,
            oindex2dimensions::Vector{Int64} ,
            irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
            iteration::Int64 ,
            broadening_distribution::String ,
            spectral_broadening::Float64 ,
            iteration_scale::Float64 ,
            K_factor::Float64 )

    # ground multiplets and partition function
    G2R_0::Dict{IntIrrep,Int64} = Dict( 
        G => length(filter(iszero,E))
        for (G,(E,U)) in irrEU 
    )
    partition = sum( R*oindex2dimensions[I]*(S+1) for ((_,I,S),R) in G2R_0 )

    # iterate over irreducible matrices
    for ((G_u,G_a,G_v),redmat) in A 

        # irrep quantum numbers
        (_,I_u,S_u) = G_u 
        (_,I_a,S_a) = G_a

        @views energies_Gu = irrEU[G_u][1]
        @views energies_Gv = irrEU[G_v][1]

        # clebsch-gordan contribution 
        dim_u = oindex2dimensions[I_u]*(S_u+1)
        dim_a = oindex2dimensions[I_a]*(S_a+1)
        cg = dim_u/dim_a

        # iterate over multiplets in irrep combination
        for r_a in axes(redmat,2)

            m_a = (G_a...,r_a)
            correlation_matrix = correlation_dict[m_a]

            for r_u in axes(redmat,1),
                r_v in axes(redmat,3)

                # coefficient from reduced matrix element
                redmatel = abs2(redmat[r_u,r_a,r_v])

                # excitation energy
                e_u = energies_Gu[r_u]
                e_v = energies_Gv[r_v]
                e_diff = e_u-e_v

                # total contribution 
                w = cg*redmatel/partition

                # particle contribution: m_v ground, positive peak
                if iszero(e_v) && e_diff>0
                    #contribution = w*P((K_factor-e_diff)*iteration_scale,
                    #                                           (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*iteration_scale);
                    #                                           distribution=broadening_distribution,
                    #                                           E=e_diff*iteration_scale,
                    #                                           omega=K_factor*iteration_scale)
                    correlation_matrix[end-iteration,2] += w*P((K_factor-e_diff)*iteration_scale,
                                                               (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*iteration_scale);
                                                               distribution=broadening_distribution,
                                                               E=e_diff*iteration_scale,
                                                               omega=K_factor*iteration_scale)
                # hole contribution: m_u ground, negative peak
                elseif iszero(e_u) && e_diff<0
                    correlation_matrix[iteration+1,2] += w*P((-K_factor-e_diff)*iteration_scale,
                                                             (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*iteration_scale);
                                                             distribution=broadening_distribution,
                                                             E=e_diff*iteration_scale,
                                                             omega=K_factor*iteration_scale)
                end
            end
        end
    end

end
function add_spectral_contribution_Tnonzero!(
        correlation_dict::Dict{IntMultiplet,Matrix{Float64}} ,
        A::Dict{ IntTripleG , Array{ComplexF64,3} } ,
        oindex2dimensions::Vector{Int64} ,
        irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
        iteration::Int64 ,
        broadening_distribution::String ,
        spectral_broadening::Float64 ,
        iteration_scale::Float64 ,
        K_factor::Float64 ,
        T::Float64 ,
        limit_shell::Bool ,
        extra_iterations::Int64 )

    # avoid extra computation
    scale_over_T = iteration_scale/T

    # partition function
    partition = sum( oindex2dimensions[I]*(S+1)*exp(-e*scale_over_T)
                     for ((_,I,S),(E,_)) in irrEU
                     for e in E )

    # iterate over irreducible matrices
    for ((G_u,G_a,G_v),redmat) in A 

        # irrep quantum numbers
        (_,I_u,S_u) = G_u 
        (_,I_a,S_a) = G_a

        # energy vectors
        @views energies_Gu = irrEU[G_u][1]
        @views energies_Gv = irrEU[G_v][1]

        # clebsch-gordan contribution 
        dim_u = oindex2dimensions[I_u]*(S_u+1)
        dim_a = oindex2dimensions[I_a]*(S_a+1)
        cg = dim_u/dim_a

        # iterate over multiplets in irrep combination
        for r_a in axes(redmat,2)

            m_a = (G_a...,r_a)
            correlation_matrix = correlation_dict[m_a]

            for r_u in axes(redmat,1),
                r_v in axes(redmat,3)

                # coefficient from reduced matrix element
                @inbounds redmatel = abs2(redmat[r_u,r_a,r_v])

                # multiplet energies
                @inbounds e_u = energies_Gu[r_u]
                @inbounds e_v = energies_Gv[r_v]
                e_diff = e_u-e_v

                # total contribution 
                boltzmann = exp(-e_u*scale_over_T)+exp(-e_v*scale_over_T)
                w = cg*redmatel*boltzmann/partition

                if e_diff>0
                    #correlation_matrix[end-iteration,2] += w*P((K_factor-e_diff)*iteration_scale,
                    #                                           (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*iteration_scale);
                    #                                           distribution=broadening_distribution,
                    #                                           E=e_diff*iteration_scale,
                    #                                           omega=K_factor*iteration_scale)
                    correlation_matrix[end-iteration,2] += w*C( K_factor*iteration_scale ,
                                                               e_diff*iteration_scale ,
                                                               spectral_broadening ,
                                                               broadening_distribution ;
                                                               tq=T)
                elseif e_diff<0
                    #correlation_matrix[iteration+1,2] += w*P((-K_factor-e_diff)*iteration_scale,
                    #                                         (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*iteration_scale);
                    #                                         distribution=broadening_distribution,
                    #                                         E=e_diff*iteration_scale,
                    #                                         omega=K_factor*iteration_scale)
                    correlation_matrix[iteration+1,2] += w*C(   -K_factor*iteration_scale ,
                                                               e_diff*iteration_scale ,
                                                               spectral_broadening ,
                                                               broadening_distribution ;
                                                               tq=T)
                end

                if limit_shell
                    #broad_dist = "lorentzian"
                    #eta_abs = 0.3*iteration_scale
                    #println("-------------")
                    #println()
                    for i in (iteration+1):(iteration+extra_iterations)

                        extra_iteration_energy = -correlation_matrix[i+1,1]
                        relative_scale = extra_iteration_energy/(K_factor*iteration_scale)

                        if e_diff>0
                            correlation_matrix[end-i,2] += w*C(extra_iteration_energy,
                                                               e_diff*iteration_scale,
                                                               spectral_broadening,
                                                               broadening_distribution;
                                                               tq=T)
                            #contribution = C(extra_iteration_energy,
                            #                 e_diff*iteration_scale,
                            #                 eta_abs,
                            #                 broad_dist)
                            #@show e_diff
                            #@show relative_scale
                            #@show contribution
                            #println()
                            #correlation_matrix[end-i,2] += w*P(K_factor*relative_scale-e_diff)*iteration_scale,
                            #                                   iteration_scale;
                            #                                   distribution="lorentzian",
                            #                                   E=e_diff*iteration_scale,
                            #                                   omega=K_factor*relative_scale*iteration_scale)
                        elseif e_diff<0
                            correlation_matrix[i+1,2] += w*C(-extra_iteration_energy,
                                                             e_diff*iteration_scale,
                                                             spectral_broadening,
                                                             broadening_distribution;
                                                             tq=T)
                            #contribution = C(-extra_iteration_energy,
                            #                  e_diff*iteration_scale,
                            #                  eta_abs,
                            #                  broad_dist)
                            #@show e_diff
                            #@show relative_scale
                            #@show contribution
                            #println()
                            #correlation_matrix[i+1,2] += w*P((-K_factor*relative_scale-e_diff)*iteration_scale,
                            #                                 iteration_scale;
                            #                                 distribution="lorentzian",
                            #                                 E=e_diff*iteration_scale,
                            #                                 omega=K_factor*relative_scale*iteration_scale)
                        end
                    end
                end
            end
        end
    end
end
function add_spectral_contribution_density_matrix_T0!(
        correlation_dict::Dict{IntMultiplet,Matrix{Float64}} ,
        A::Dict{ IntTripleG , Array{ComplexF64,3} } ,
        oindex2dimensions::Vector{Int64} ,
        irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
        iteration::Int64 ,
        broadening_distribution::String ,
        spectral_broadening::Float64 ,
        iteration_scale::Float64 ,
        K_factor::Float64 ,
        density_matrix::Dict{IntIrrep,Matrix{Float64}} ,
        multiplets_kept::Vector{IntMultiplet} ,
        multiplets_discarded::Vector{IntMultiplet} ,
        L::Float64 ,
        scale::Float64 )

    total_iterations = size(collect(values(correlation_dict))[1],1)÷2 - 1
    iterscales = map( n->iterscale(scale,L,n) , 0:total_iterations )

    length(multiplets_discarded)==0 && return
    GG_kept      = Set( m[1:3] for m in multiplets_kept )
    GG_discarded = Set( m[1:3] for m in multiplets_discarded )
    G2rr_kept = Dict(
        G=>(1:maximum([ m[4] for m in multiplets_kept if m[1:3]==G ]))
        for G in GG_kept
    )
    G2rr_discarded = Dict(
        G=>(minimum([ m[4] for m in multiplets_discarded if m[1:3]==G ]):maximum([ m[4] for m in multiplets_discarded if m[1:3]==G ]))
        for G in GG_discarded
    )

    if iteration==100
        println( "Iteration $(iteration) with (i+1)=$(iteration+1)" )
        println()
        @show keys(density_matrix)
        #@show collect(values(density_matrix))[1]
        tr_plus = sum(  (N>(iteration+1) ? tr(d) : 0.0) for ((N,_,_),d) in density_matrix )
        tr_minus = sum( (N<(iteration+1) ? tr(d) : 0.0) for ((N,_,_),d) in density_matrix )
        @show tr_plus,tr_minus
        println()
        mults_kept_plus  = length([ m for m in multiplets_kept if m[1]>(iteration+1) ])
        mults_kept_minus = length([ m for m in multiplets_kept if m[1]<(iteration+1) ])
        mults_discarded_plus  = length([ m for m in multiplets_discarded if m[1]>(iteration+1) ])
        mults_discarded_minus = length([ m for m in multiplets_discarded if m[1]<(iteration+1) ])
        @show mults_kept_plus,mults_kept_minus
        @show mults_discarded_plus,mults_discarded_minus
        println()
    end

    particle_sum = 0.0
    hole_sum = 0.0

    e_limit = 100.0

    # iterate over irreducible matrices
    for ((G_u,G_a,G_v),redmat) in A 

        # irrep quantum numbers
        (_,I_u,S_u) = G_u 
        (_,I_a,S_a) = G_a

        @views energies_Gu = irrEU[G_u][1]
        @views energies_Gv = irrEU[G_v][1]

        u_ground = haskey(density_matrix,G_u)
        v_ground = haskey(density_matrix,G_v)
        (u_ground || v_ground) || continue
        u_ground && (@views dm_u = density_matrix[G_u])
        v_ground && (@views dm_v = density_matrix[G_v])

        # clebsch-gordan contribution 
        dim_u = oindex2dimensions[I_u]*(S_u+1)
        dim_a = oindex2dimensions[I_a]*(S_a+1)
        cg = dim_u/dim_a

        if iteration==100
            println( "Redmat contribution:" )
            @show G_u,G_v
            @show cg
            println()
        end

        # iterate over multiplets in irrep combination
        for r_a in axes(redmat,2)

            m_a = (G_a...,r_a)

            correlation_matrix = correlation_dict[m_a]

            # particle sum over ground states
            if (v_ground && haskey(G2rr_discarded,G_u))
                @inbounds for r_v1 in axes(dm_v,2),
                              r_v2 in axes(dm_v,2)

                    # density matrix element
                    dmel = dm_v[r_v2,r_v1]
                    isapprox(dmel,zero(dmel)) && continue

                    if iteration==100
                        @show r_v1,r_v2,dmel
                        println()
                    end

                    # energies of the kept ("ground") states
                    e_v1 = energies_Gv[r_v1]
                    e_v2 = energies_Gv[r_v2]

                    # iterate over excited (discarded) states
                    for r_u in G2rr_discarded[G_u]

                        redmatel = conj(redmat[r_u,r_a,r_v1])*redmat[r_u,r_a,r_v2]
                        iszero(redmatel) && continue

                        e_u = energies_Gu[r_u]
                        e_diff = e_u-0.5*(e_v1+e_v2) 

                        #e_diff>e_limit && continue

                        # total weight 
                        w = real(cg*redmatel*dmel)
                        iszero(w) && continue

                        if iteration==100
                            @show r_u
                            @show e_diff, e_v1, e_v2, 0.5(e_v1+e_v2)
                            @show redmatel
                            @show w
                            println()
                        end

                        @inbounds for (i,i_scale) in enumerate(iterscales)
                            correlation_matrix[(end-(i-1)),2] += w*C( 
                                correlation_matrix[(end-(i-1)),1],
                                e_diff*iteration_scale,
                                (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*correlation_matrix[(end-(i-1)),1]),
                                broadening_distribution 
                            )
                            if i==length(iterscales) && iteration==100
                                particle_sum += w*C( correlation_matrix[(end-(i-1)),1],
                                    e_diff*iteration_scale,
                                    (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*correlation_matrix[(end-(i-1)),1]),
                                    broadening_distribution 
                                )
                            end
                        end
                    end
                end
            end

            # hole sum over ground states if there are discarded v states
            if (u_ground && haskey(G2rr_discarded,G_v))
                @inbounds for r_u1 in axes(dm_u,1),
                              r_u2 in axes(dm_u,1)

                    # density matrix element
                    dmel = dm_u[r_u2,r_u1]
                    isapprox(dmel,zero(dmel)) && continue

                    # energies of the kept ("ground") states
                    e_u1 = energies_Gu[r_u1]
                    e_u2 = energies_Gu[r_u2]

                    if iteration==100
                        @show r_u1,r_u2,dmel
                        println()
                    end

                    # iterate over excited (discarded) states
                    for r_v in G2rr_discarded[G_v]

                        redmatel = redmat[r_u1,r_a,r_v]*conj(redmat[r_u2,r_a,r_v])
                        isapprox(redmatel,zero(redmatel)) && continue

                        e_v = energies_Gv[r_v]
                        e_diff = 0.5*(e_u1+e_u2) - e_v

                        #-e_diff>e_limit && continue

                        # total weight 
                        w = cg*redmatel*dmel
                        iszero(w) && continue

                        if iteration==100
                            @show r_v
                            @show e_diff
                            @show redmatel
                            @show w
                            println()
                        end

                        @inbounds for (i,i_scale) in enumerate(iterscales)
                            correlation_matrix[i,2] += w*C( 
                                correlation_matrix[i,1],
                                e_diff*iteration_scale,
                                (broadening_distribution=="loggaussian" ? spectral_broadening : -spectral_broadening*correlation_matrix[i,1]),
                                broadening_distribution 
                            )
                            if i==length(iterscales) && iteration==100
                                hole_sum += w*C( 
                                    correlation_matrix[i,1],
                                    e_diff*iteration_scale,
                                    (broadening_distribution=="loggaussian" ? spectral_broadening : -spectral_broadening*correlation_matrix[i,1]),
                                    broadening_distribution 
                                )
                            end
                        end

                    end
                end
            end

        end
    end
    @show particle_sum,hole_sum

end
function add_spectral_contribution_density_matrix_Tnonzero!(
        correlation_dict::Dict{IntMultiplet,Matrix{Float64}} ,
        A::Dict{ IntTripleG , Array{ComplexF64,3} } ,
        oindex2dimensions::Vector{Int64} ,
        irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
        iteration::Int64 ,
        broadening_distribution::String ,
        spectral_broadening::Float64 ,
        iteration_scale::Float64 ,
        K_factor::Float64 ,
        density_matrix::Dict{IntIrrep,Matrix{Float64}} ,
        multiplets_kept::Vector{IntMultiplet} ,
        multiplets_discarded::Vector{IntMultiplet} ,
        L::Float64 ,
        scale::Float64 ,
        half_weight_idx::Int64 ,
        half_weight_energy::Float64 ,
        T::Float64 )

    total_iterations = size(collect(values(correlation_dict))[1],1)÷2 - 1
    iterscales = map( n->iterscale(scale,L,n) , 0:total_iterations )
    smearing_parameter = iterscales[half_weight_idx]

    negative_energy_mask = map( x->x<=total_iterations+1 , 1:2*(total_iterations+1) )
    positive_energy_mask = map( x->x>total_iterations+1  , 1:2*(total_iterations+1) )

    high_energy_iterator = 0:half_weight_idx
    low_energy_iterator  = (half_weight_idx+1):total_iterations

    length(multiplets_discarded)==0 && return
    GG_kept      = Set( m[1:3] for m in multiplets_kept )
    GG_discarded = Set( m[1:3] for m in multiplets_discarded )
    G2rr_kept = Dict(
        G=>sort([ m[4] for m in multiplets_kept if m[1:3]==G ])
        for G in GG_kept
    )
    G2rr_discarded = Dict(
        G=>sort([ m[4] for m in multiplets_discarded if m[1:3]==G ])
        for G in GG_discarded
    )
    G2rr_all = mergewith( vcat , G2rr_kept , G2rr_discarded )

    # iterate over irreducible matrices
    for ((G_u,G_a,G_v),redmat) in A 

        # irrep quantum numbers
        (_,I_u,S_u) = G_u 
        (_,I_a,S_a) = G_a

        @views energies_Gu = irrEU[G_u][1]
        @views energies_Gv = irrEU[G_v][1]

        @views dm_u = density_matrix[G_u]
        @views dm_v = density_matrix[G_v]

        # clebsch-gordan contribution 
        dim_u = oindex2dimensions[I_u]*(S_u+1)
        dim_a = oindex2dimensions[I_a]*(S_a+1)
        cg = dim_u/dim_a

        # iterate over multiplets in irrep combination
        for r_a in axes(redmat,2)

            m_a = (G_a...,r_a)

            correlation_matrix = correlation_dict[m_a]

            # particle sum over kept states
            #
            # m_v1 in K => m_v2 in K => m_u in D
            if haskey(G2rr_kept,G_v) && haskey(G2rr_discarded,G_u)

                G2rr_kept_v = G2rr_kept[G_v]
                G2rr_disc_u = G2rr_discarded[G_u]

                @inbounds for r_v1 in G2rr_kept_v,
                              r_v2 in G2rr_kept_v

                    # density matrix element
                    dmel = dm_v[r_v2,r_v1]
                    iszero(dmel) && continue

                    # energies of the kept ("ground") states
                    e_v1 = energies_Gv[r_v1]
                    e_v2 = energies_Gv[r_v2]

                    # iterate over excited (discarded) states
                    for r_u in G2rr_disc_u

                        redmatel = conj(redmat[r_u,r_a,r_v1])*redmat[r_u,r_a,r_v2]
                        iszero(redmatel) && continue

                        e_u = energies_Gu[r_u]
                        e_diff = e_u-0.5*(e_v1+e_v2) 
                        e_diff_scaled = e_diff*iteration_scale

                        #e_diff>e_limit && continue

                        # total weight 
                        w = real(cg*redmatel*dmel)
                        iszero(w) && continue

                        add_excitation_contribution!( 
                            correlation_matrix ,
                            w ,
                            e_diff_scaled ,
                            high_energy_iterator ,
                            low_energy_iterator ,
                            broadening_distribution ,
                            spectral_broadening ,
                            half_weight_energy ;
                            tq=T
                        )
                        #correlation_matrix[positive_energy_mask,2] += map( 
                        #    omega-> omega>=half_weight_energy ? 
                        #            w*C(omega,
                        #                e_diff_scaled,
                        #                (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*e_diff_scaled),
                        #                broadening_distribution) :
                        #            w*C(omega,
                        #                e_diff_scaled,
                        #                (broadening_distribution=="loggaussian" ? spectral_broadening : half_weight_energy),
                        #                broadening_distribution),
                        #    correlation_matrix[positive_energy_mask,1]
                        #)
                        #@inbounds for n in high_energy_iterator
                        #    omega = correlation_matrix[(end-n),1] 
                        #    correlation_matrix[(end-n),2] += w*C( 
                        #        omega,
                        #        e_diff_scaled,
                        #        (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*omega),
                        #        broadening_distribution 
                        #    )
                        #end
                        #@inbounds for n in low_energy_iterator
                        #    omega = correlation_matrix[(end-n),1] 
                        #    correlation_matrix[(end-n),2] += w*C( 
                        #        omega,
                        #        e_diff_scaled,
                        #        half_weight_energy,
                        #        "gaussian" 
                        #    )
                        #end
                        #@inbounds for (i,i_scale) in enumerate(iterscales)
                        #    if n<=half_weight_idx
                        #        correlation_matrix[(end-(i-1)),2] += w*C( 
                        #            correlation_matrix[(end-(i-1)),1],
                        #            e_diff_scaled,
                        #            (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*correlation_matrix[(end-(i-1)),1]),
                        #            broadening_distribution 
                        #        )
                        #    else
                        #        correlation_matrix[(end-(i-1)),2] += w*C( 
                        #            correlation_matrix[(end-(i-1)),1],
                        #            e_diff*iteration_scale,
                        #            smearing_parameter,
                        #            "gaussian" 
                        #        )
                        #    end
                        #end
                    end
                end
            end
            # particle sum over discarded states
            #
            # m_v1 in D => m_v2=m_v1 in D => m_u in (K,D)
            if haskey(G2rr_discarded,G_v)

                G2rr_disc_v = G2rr_discarded[G_v]
                G2rr_all_u = G2rr_all[G_u]

                @inbounds for r_v in G2rr_disc_v

                    # density matrix element
                    dmel = dm_v[r_v,r_v]
                    iszero(dmel) && continue

                    # energies of the "ground" states
                    e_v = energies_Gv[r_v]

                    # iterate over all intermediate states
                    for r_u in G2rr_all_u

                        redmatel = abs2(redmat[r_u,r_a,r_v])
                        iszero(redmatel) && continue

                        e_u = energies_Gu[r_u]
                        e_diff = e_u-e_v
                        e_diff_scaled = e_diff*iteration_scale

                        # total weight 
                        w = real(cg*redmatel*dmel)
                        iszero(w) && continue

                        add_excitation_contribution!( 
                            correlation_matrix ,
                            w ,
                            e_diff_scaled ,
                            high_energy_iterator ,
                            low_energy_iterator ,
                            broadening_distribution ,
                            spectral_broadening ,
                            half_weight_energy ;
                            tq=T
                        )
                        #if sign==1
                        #    @inbounds for n in high_energy_iterator
                        #        omega = correlation_matrix[(end-n),1]
                        #        correlation_matrix[(end-n),2] += w*C( 
                        #            omega,
                        #            e_diff*iteration_scale,
                        #            (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*omega),
                        #            broadening_distribution 
                        #        )
                        #    end
                        #    @inbounds for n in low_energy_iterator
                        #        omega = correlation_matrix[(end-n),1]
                        #        correlation_matrix[(end-n),2] += w*C( 
                        #            omega,
                        #            e_diff*iteration_scale,
                        #            half_weight_energy,
                        #            "gaussian" 
                        #        )
                        #    end
                        #elseif sign==-1
                        #    @inbounds for n in high_energy_iterator
                        #        omega = correlation_matrix[n+1,1]
                        #        correlation_matrix[n+1,2] += w*C( 
                        #            omega,
                        #            e_diff*iteration_scale,
                        #            (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*omega),
                        #            broadening_distribution 
                        #        )
                        #    end
                        #    @inbounds for n in low_energy_iterator
                        #        omega = correlation_matrix[n+1,1]
                        #        correlation_matrix[n+1,2] += w*C( 
                        #            omega,
                        #            e_diff*iteration_scale,
                        #            half_weight_energy,
                        #            "gaussian" 
                        #        )
                        #    end
                        #else
                        #@inbounds for (i,i_scale) in enumerate(iterscales)
                        #    n = i-1
                        #    if n<=half_weight_idx
                        #        correlation_matrix[(end-(i-1)),2] += w*C( 
                        #            sign*correlation_matrix[(end-(i-1)),1],
                        #            e_diff*iteration_scale,
                        #            (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*correlation_matrix[(end-(i-1)),1]),
                        #            broadening_distribution 
                        #        )
                        #    else
                        #        correlation_matrix[(end-(i-1)),2] += w*C( 
                        #            sign*correlation_matrix[(end-(i-1)),1],
                        #            e_diff*iteration_scale,
                        #            smearing_parameter,
                        #            "gaussian" 
                        #        )
                        #    end
                        #end
                    end
                end
            end

            # hole sum over kept states
            if haskey(G2rr_kept,G_u) && haskey(G2rr_discarded,G_v)

                G2rr_kept_u = G2rr_kept[G_u]
                G2rr_disc_v = G2rr_discarded[G_v]

                @inbounds for r_u1 in G2rr_kept_u,
                              r_u2 in G2rr_kept_u

                    # density matrix element
                    dmel = dm_u[r_u2,r_u1]
                    iszero(dmel) && continue

                    # energies of the kept ("ground") states
                    e_u1 = energies_Gu[r_u1]
                    e_u2 = energies_Gu[r_u2]

                    # iterate over excited (discarded) states
                    for r_v in G2rr_disc_v

                        redmatel = redmat[r_u1,r_a,r_v]*conj(redmat[r_u2,r_a,r_v])
                        iszero(redmatel) && continue

                        e_v = energies_Gv[r_v]
                        e_diff = 0.5*(e_u1+e_u2) - e_v
                        e_diff_scaled = e_diff*iteration_scale

                        # total weight 
                        w = cg*redmatel*dmel
                        iszero(w) && continue

                        add_excitation_contribution!( 
                            correlation_matrix ,
                            w ,
                            e_diff_scaled ,
                            high_energy_iterator ,
                            low_energy_iterator ,
                            broadening_distribution ,
                            spectral_broadening ,
                            half_weight_energy ;
                            tq=T
                        )
                        #@inbounds for n in high_energy_iterator
                        #    omega = correlation_matrix[n+1,1]
                        #    correlation_matrix[n+1,2] += w*C(
                        #        omega,
                        #        e_diff_scaled,
                        #        (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*abs(e_diff_scaled)),
                        #        broadening_distribution
                        #    )
                        #end
                        #@inbounds for n in low_energy_iterator
                        #    omega = correlation_matrix[n+1,1]
                        #    correlation_matrix[n+1,2] += w*C(
                        #        omega,
                        #        e_diff_scaled,
                        #        half_weight_energy,
                        #        "gaussian"
                        #    )
                        #end
                        #correlation_matrix[negative_energy_mask,2] += map( 
                        #    omega-> -omega>=half_weight_energy ? 
                        #            w*C(omega,
                        #                e_diff_scaled,
                        #                (broadening_distribution=="loggaussian" ? spectral_broadening : -spectral_broadening*e_diff_scaled),
                        #                broadening_distribution) :
                        #            w*C(omega,
                        #                e_diff_scaled,
                        #                (broadening_distribution=="loggaussian" ? spectral_broadening : half_weight_energy),
                        #                broadening_distribution),
                        #    correlation_matrix[negative_energy_mask,1]
                        #)
                        #@inbounds for (i,i_scale) in enumerate(iterscales)
                        #    n=i-1
                        #    if n<=half_weight_idx
                        #        correlation_matrix[i,2] += w*C( 
                        #            correlation_matrix[i,1],
                        #            e_diff*iteration_scale,
                        #            (broadening_distribution=="loggaussian" ? spectral_broadening : -spectral_broadening*correlation_matrix[i,1]),
                        #            broadening_distribution 
                        #        )
                        #    else
                        #        correlation_matrix[i,2] += w*C( 
                        #            correlation_matrix[i,1],
                        #            e_diff*iteration_scale,
                        #            smearing_parameter,
                        #            "gaussian" 
                        #        )
                        #    end
                        #end

                    end
                end
            end
            # hole sum over discarded states
            if haskey(G2rr_discarded,G_u)

                G2rr_disc_u = G2rr_discarded[G_u]
                G2rr_all_v = G2rr_all[G_v]

                @inbounds for r_u in G2rr_disc_u,
                              r_u in G2rr_disc_u

                    # density matrix element
                    dmel = dm_u[r_u,r_u]
                    iszero(dmel) && continue

                    # energies of the kept ("ground") states
                    e_u = energies_Gu[r_u]

                    # iterate over excited (discarded) states
                    for r_v in G2rr_all_v

                        redmatel = abs2(redmat[r_u,r_a,r_v])
                        iszero(redmatel) && continue

                        e_v = energies_Gv[r_v]
                        e_diff = e_u - e_v
                        e_diff_scaled = e_diff*iteration_scale

                        # total weight 
                        w = cg*redmatel*dmel
                        iszero(w) && continue

                        add_excitation_contribution!( 
                            correlation_matrix ,
                            w ,
                            e_diff_scaled ,
                            high_energy_iterator ,
                            low_energy_iterator ,
                            broadening_distribution ,
                            spectral_broadening ,
                            half_weight_energy ;
                            tq=T
                        )
                        #if sign==-1
                        #    @inbounds for n in high_energy_iterator
                        #        omega = correlation_matrix[n+1,1]
                        #        correlation_matrix[n+1,2] += w*C(
                        #            omega,
                        #            e_diff_scaled,
                        #            (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*abs(omega)),
                        #            broadening_distribution 
                        #        )
                        #    end
                        #    @inbounds for n in low_energy_iterator
                        #        omega = correlation_matrix[n+1,1]
                        #        correlation_matrix[n+1,2] += w*C(
                        #            omega,
                        #            e_diff_scaled,
                        #            half_weight_energy,
                        #            "gaussian" 
                        #        )
                        #    end
                        #elseif sign==1
                        #    @inbounds for n in high_energy_iterator
                        #        omega = correlation_matrix[end-n,1]
                        #        correlation_matrix[end-n,2] += w*C(
                        #            omega,
                        #            e_diff_scaled,
                        #            (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*abs(omega)),
                        #            broadening_distribution 
                        #        )
                        #    end
                        #    @inbounds for n in low_energy_iterator
                        #        omega = correlation_matrix[end-n,1]
                        #        correlation_matrix[end-n,2] += w*C(
                        #            omega,
                        #            e_diff_scaled,
                        #            half_weight_energy,
                        #            "gaussian" 
                        #        )
                        #    end
                        #end
                        #@inbounds for (i,i_scale) in enumerate(iterscales)
                        #    n=i-1
                        #    if n<=half_weight_idx
                        #        correlation_matrix[i,2] += w*C( 
                        #            sign*correlation_matrix[i,1],
                        #            e_diff*iteration_scale,
                        #            (broadening_distribution=="loggaussian" ? spectral_broadening : -spectral_broadening*correlation_matrix[i,1]),
                        #            broadening_distribution 
                        #        )
                        #    else
                        #        correlation_matrix[i,2] += w*C( 
                        #            sign*correlation_matrix[i,1],
                        #            e_diff*iteration_scale,
                        #            smearing_parameter,
                        #            "gaussian" 
                        #        )
                        #    end
                        #end
                    end
                end
            end

        end
    end

end
function add_excitation_contribution!( correlation_matrix::Matrix{Float64} ,
                                       w ,
                                       e_diff_scaled::Float64 ,
                                       high_energy_iterator ,
                                       low_energy_iterator ,
                                       broadening_distribution::String ,
                                       spectral_broadening::Float64 ,
                                       half_weight_energy::Float64 ;
                                       tq::Float64=0.0 )

    open( "spectral/peaks.dat" , "a") do f
        write( f , "$(e_diff_scaled) $w\n" )
    end

    high_energy_broadening = broadening_distribution in ("loggaussian","interpolated") ? 
                             spectral_broadening : 
                             spectral_broadening*abs(e_diff_scaled)
    low_energy_distribution = broadening_distribution=="interpolated" ?
                              "interpolated" :
                              "gaussian"
    low_energy_broadening = broadening_distribution=="interpolated" ?
                            spectral_broadening :
                            half_weight_energy

    if e_diff_scaled>0
        @inbounds for n in high_energy_iterator
            omega = correlation_matrix[end-n,1]
            correlation_matrix[end-n,2] += w*C(
                omega,
                e_diff_scaled,
                high_energy_broadening,
                broadening_distribution;
                tq=tq
            )
        end
        @inbounds for n in low_energy_iterator
            omega = correlation_matrix[end-n,1]
            correlation_matrix[end-n,2] += w*C(
                omega,
                e_diff_scaled,
                low_energy_broadening,
                low_energy_distribution;
                tq=tq
            )
        end
    else
        @inbounds for n in high_energy_iterator
            omega = correlation_matrix[n+1,1]
            correlation_matrix[n+1,2] += w*C(
                omega,
                e_diff_scaled,
                high_energy_broadening,
                broadening_distribution;
                tq=tq
            )
        end
        @inbounds for n in low_energy_iterator
            omega = correlation_matrix[n+1,1]
            correlation_matrix[n+1,2] += w*C(
                omega,
                e_diff_scaled,
                low_energy_broadening,
                low_energy_distribution;
                tq=tq
            )
        end
    end
end
function save_correlation_spectral_decomposition(
            spectral_functions::Dict{String,Dict{IntMultiplet,Matrix{Float64}}} ,
            label::String ,
            z::Float64 ;
            spectraldir::String="spectral"
    )

    spectral_function = spectral_functions["spectral"]

    # iterate over excitation multiplets
    for (orbital_idx,(m_a,spectral_a)) in enumerate(spectral_function)

        N = size(spectral_a,1)÷2
        iterations = N - 1

        # even and odd spectral functions
        start_even = iseven(iterations) ? iterations+3 : iterations+2
        start_odd  = iseven(iterations) ? iterations+2 : iterations+3
        spectral_function_even = vcat(
            spectral_a[1:2:iterations+1,:],
            spectral_a[start_even:2:end,:]
        )
        spectral_function_odd = vcat(
            spectral_a[2:2:iterations+1,:],
            spectral_a[start_odd:2:end,:]
        )

        # interpolation
        omegas_evenodd = vcat(
            spectral_a[2:iterations,1] ,
            spectral_a[iterations+3:end-1,1]
        )
        interpolator_even = linear_interpolation(spectral_function_even[:,1],spectral_function_even[:,2])
        interpolator_odd  = linear_interpolation(spectral_function_odd[:,1] ,spectral_function_odd[:,2])
        spectral_function_evenodd = zeros(Float64,length(omegas_evenodd),2)
        spectral_function_evenodd[:,1] .= omegas_evenodd
        spectral_function_evenodd[:,2] .= map( x->0.5*(interpolator_even(x)+interpolator_odd(x)) , omegas_evenodd )

        # write to file
        orbital_header = "# Excitation multiplet: $m_a , index=$orbital_idx\n"
        write_spectral_function( 
            spectral_filename(label,z=z,tail="_even",orbital=orbital_idx,spectraldir=spectraldir) ,
            spectral_function_even ,
            orbitalresolved_header=orbital_header
        )
        write_spectral_function(
            spectral_filename(label,z=z,tail="_odd",orbital=orbital_idx,spectraldir=spectraldir) ,
            spectral_function_odd ,
            orbitalresolved_header=orbital_header
        )
        write_spectral_function(
            spectral_filename(label,z=z,orbital=orbital_idx,spectraldir=spectraldir) ,
            spectral_function_evenodd ,
            orbitalresolved_header=orbital_header
        )

    end

end

iterscale( scale_0 , L , n ) = scale_0 * L^(-(n-2)/2.0)

function update_operator( redmat_iaj::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} },
                          multiplets_a_block::Vector{NTuple{4,Int64}}, 
                          Ksum_o_array::Array{ComplexF64,6} ,
                          Ksum_s_array::Array{ComplexF64,6} ,
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
                iszero(K) && continue

                # total prefactor
                sign_times_K = sign * K

                # < G_i || f^\dagger_{G_a} || G_j >
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

function compute_correlation_spectral_decomposition(
            L::Float64 ,
            iterations::Int64 ,
            scale::Float64 ,
            spectral_broadening::Float64 ,
            label::String ,
            z::Float64 ;
            broadening_distribution::String="loggaussian",
            K_factor::Float64=2.0 ,
            width_two::Bool=false ,
            correlation_type::String="spectral" ,
            triple_excitation::Bool=false ,
            T::Float64=0.0 ,
            betabar::Float64=1.0 )
    if iszero(T)
        compute_correlation_spectral_decomposition_T0(
            L,
            iterations,
            scale,
            spectral_broadening,
            label,
            z;
            broadening_distribution=broadening_distribution,
            K_factor=K_factor,
            width_two=width_two,
            correlation_type=correlation_type,
            triple_excitation=triple_excitation)
    elseif !iszero(T)
        compute_correlation_spectral_decomposition_Tnonzero(
            L,
            iterations,
            scale,
            spectral_broadening,
            label,
            z,
            T,
            betabar;
            broadening_distribution=broadening_distribution,
            K_factor=K_factor,
            width_two=width_two,
            correlation_type=correlation_type,
            triple_excitation=triple_excitation)
    end
end
function compute_correlation_spectral_decomposition_T0( 
            L::Float64 ,
            iterations::Int64 ,
            scale::Float64 ,
            spectral_broadening::Float64 ,
            label::String ,
            z::Float64 ;
            broadening_distribution::String="loggaussian",
            K_factor::Float64=2.0 ,
            width_two::Bool=false ,
            correlation_type::String="spectral" ,
            triple_excitation::Bool=false )

    # even and odd iteration ranges
    even_iterator = 2:2:iterations
    odd_iterator  = 1:2:iterations
    # the chosen energies are 
    #
    #   e = K_factor * omega_N
    #
    omegas_odd = sort([ 
        ( K_factor * sign * scale * L^(-(N-2)/2.0) )
        for sign in [-1.0,1.0]
        for N in odd_iterator
    ])
    omegas_even = sort([ 
        ( K_factor * sign * scale * L^(-(N-2)/2.0) )
        for sign in [-1.0,1.0]
        for N in even_iterator
    ])
    scalings = [ iterscale(scale,L,n) for n in 1:iterations ]

    # directory where the data is stored and names of multiplets
    correlation_type=="spectral" && (data_dir = "spectral/peaks/particle/z$z")
    multiplets_names = readdir("$(data_dir)/n0")
    parsemultiplet(multiplet_name) = (map(x->parse(Int64,x),split(multiplet_name,"_"))...,)
    excitation_multiplets = Set(map(parsemultiplet,multiplets_names))

    # spectral function
    spectral_odd  = [0.0 for _ in omegas_odd]
    spectral_even = [0.0 for _ in omegas_even]
    spectral_odd  = Dict{IntMultiplet,Vector{Float64}}( 
        multiplet=>copy(spectral_odd) for multiplet in excitation_multiplets
    )
    spectral_even  = Dict{IntMultiplet,Vector{Float64}}( 
        multiplet=>copy(spectral_even) for multiplet in excitation_multiplets
    )

    iteration_energy_rescaled = K_factor
    for n in 1:iterations

        iteration_scale = scalings[n]

        iteration_data_dir = "$(data_dir)/n$n"
        for multiplet_name in readdir(iteration_data_dir)

            excitation_multiplet = parsemultiplet(multiplet_name)
            multiplet_file = "$(iteration_data_dir)/$multiplet_name"
            multiplet_data = readdlm(multiplet_file)

            for i in axes(multiplet_data,1)

                peak_energy_rescaled = multiplet_data[i,1]
                peak_weight = multiplet_data[i,2]

                contribution_positive = 0.0
                contribution_negative = 0.0
                if broadening_distribution=="gaussian"
                    if peak_energy_rescaled>0
                        contribution_positive = peak_weight*P(iteration_energy_rescaled-peak_energy_rescaled,spectral_broadening)/iteration_scale
                    elseif peak_energy_rescaled<0
                        contribution_negative = peak_weight*P(-iteration_energy_rescaled-peak_energy_rescaled,spectral_broadening)/iteration_scale
                    end
                elseif broadening_distribution=="loggaussian" && !iszero(peak_weight)
                    if peak_energy_rescaled>0
                        contribution_positive = peak_weight*P(iteration_energy_rescaled-peak_energy_rescaled,
                                                              spectral_broadening;
                                                              distribution="loggaussian",
                                                              E=peak_energy_rescaled,
                                                              omega=iteration_energy_rescaled)/iteration_scale
                    elseif peak_energy_rescaled<0
                        contribution_negative = peak_weight*P(-iteration_energy_rescaled-peak_energy_rescaled,
                                                              spectral_broadening;
                                                              distribution="loggaussian",
                                                              E=peak_energy_rescaled,
                                                              omega=-iteration_energy_rescaled)/iteration_scale
                    end
                end
                if isodd(n)
                    j = 1 + Int64((n-1)/2)
                    spectral_odd[excitation_multiplet][end-(j-1)] += contribution_positive
                    spectral_odd[excitation_multiplet][j] += contribution_negative
                else
                    j = Int64(n/2)
                    spectral_even[excitation_multiplet][end-(j-1)] += contribution_positive
                    spectral_even[excitation_multiplet][j] += contribution_negative
                end


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
    # compute average 
    spectral_evenodd = Dict(
        excitation_multiplet => 0.5*( spectral_even_interpolated[excitation_multiplet] + 
                                      spectral_odd_interpolated[excitation_multiplet]   )
        for excitation_multiplet in excitation_multiplets
    )

    # compute spline-interpolated smooth spectral function
    omegas_evenodd_spline,spectral_evenodd_spline = spline_interpolate_spectral_function( 
        omegas_evenodd ,
        spectral_evenodd ,
        orbitalresolved=true
    )

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
function compute_correlation_spectral_decomposition_Tnonzero( 
            L::Float64 ,
            iterations::Int64 ,
            scale::Float64 ,
            spectral_broadening::Float64 ,
            label::String ,
            z::Float64 ,
            T::Float64 ,
            betabar::Float64 ;
            broadening_distribution::String="loggaussian",
            K_factor::Float64=2.0 ,
            width_two::Bool=false ,
            correlation_type::String="spectral" ,
            triple_excitation::Bool=false )

    # energy scales
    scalings = [ iterscale(scale,L,n) for n in 1:iterations ]

    # find limit shell
    _,n_limit = findmin(abs.(scalings.-betabar*T))
    println()
    println( "Last shell used for spectral function calculations: $n_limit" )

    # even and odd iteration ranges
    even_iterator = 2:2:iterations
    odd_iterator  = 1:2:iterations
    # the chosen energies are 
    #
    #   e = K_factor * omega_N
    #
    omegas_odd = sort([ 
        ( K_factor * sign * scale * L^(-(N-2)/2.0) )
        for sign in [-1.0,1.0]
        for N in odd_iterator
    ])
    omegas_even = sort([ 
        ( K_factor * sign * scale * L^(-(N-2)/2.0) )
        for sign in [-1.0,1.0]
        for N in even_iterator
    ])


    # directory where the data is stored and names of multiplets
    correlation_type=="spectral" && (data_dir = "spectral/peaks/particle/z$z")
    multiplets_names = readdir("$(data_dir)/n0")
    parsemultiplet(multiplet_name) = (map(x->parse(Int64,x),split(multiplet_name,"_"))...,)
    excitation_multiplets = Set(map(parsemultiplet,multiplets_names))

    # spectral function
    spectral_odd  = [0.0 for _ in eachindex(omegas_odd)]
    spectral_even = [0.0 for _ in eachindex(omegas_even)]
    spectral_odd  = Dict{IntMultiplet,Vector{Float64}}( 
        multiplet=>copy(spectral_odd) for multiplet in excitation_multiplets
    )
    spectral_even  = Dict{IntMultiplet,Vector{Float64}}( 
        multiplet=>copy(spectral_even) for multiplet in excitation_multiplets
    )

    iteration_energy_rescaled = K_factor
    # normal calculation 
    for n in 1:(n_limit-1)

        iteration_scale = scalings[n]

        iteration_data_dir = "$(data_dir)/n$n"
        for multiplet_name in readdir(iteration_data_dir)

            excitation_multiplet = parsemultiplet(multiplet_name)
            multiplet_file = "$(iteration_data_dir)/$multiplet_name"
            multiplet_data = readdlm(multiplet_file)

            for i in axes(multiplet_data,1)

                peak_energy_rescaled = multiplet_data[i,1]
                peak_weight = multiplet_data[i,2]

                contribution_positive = 0.0
                contribution_negative = 0.0
                if broadening_distribution=="gaussian"
                    if peak_energy_rescaled>0
                        contribution_positive = peak_weight*P(iteration_energy_rescaled-peak_energy_rescaled,spectral_broadening)/iteration_scale
                    elseif peak_energy_rescaled<0
                        contribution_negative = peak_weight*P(-iteration_energy_rescaled-peak_energy_rescaled,spectral_broadening)/iteration_scale
                    end
                elseif broadening_distribution=="loggaussian" && !iszero(peak_weight)
                    if peak_energy_rescaled>0
                        contribution_positive = peak_weight*P(iteration_energy_rescaled-peak_energy_rescaled,
                                                              spectral_broadening;
                                                              distribution="loggaussian",
                                                              E=peak_energy_rescaled,
                                                              omega=iteration_energy_rescaled)/iteration_scale
                    elseif peak_energy_rescaled<0
                        contribution_negative = peak_weight*P(-iteration_energy_rescaled-peak_energy_rescaled,
                                                              spectral_broadening;
                                                              distribution="loggaussian",
                                                              E=peak_energy_rescaled,
                                                              omega=-iteration_energy_rescaled)/iteration_scale
                    end
                end
                if isodd(n)
                    j = 1 + Int64((n-1)/2)
                    spectral_odd[excitation_multiplet][end-(j-1)] += contribution_positive
                    spectral_odd[excitation_multiplet][j] += contribution_negative
                else
                    j = Int64(n/2)
                    spectral_even[excitation_multiplet][end-(j-1)] += contribution_positive
                    spectral_even[excitation_multiplet][j] += contribution_negative
                end


            end
        end
    end
    # limit n calculation
    iteration_data_dir = "$(data_dir)/n$n_limit"
    limit_iteration_scale = scalings[n_limit]
    for multiplet_name in readdir(iteration_data_dir)

        excitation_multiplet = parsemultiplet(multiplet_name)
        multiplet_file = "$(iteration_data_dir)/$multiplet_name"
        multiplet_data = readdlm(multiplet_file)

        for n in n_limit:iterations

            iteration_scale = scalings[n]
            iteration_energy = iteration_scale*K_factor
            iteration_energy_rescaled = iteration_energy/limit_iteration_scale

            # cut excitations below limit
            if (iteration_energy-minimum(abs.(multiplet_data[:,1].*limit_iteration_scale)))<spectral_broadening*iteration_scale
                @show minimum(abs.(multiplet_data[:,1]))
                println( "Lowest n and shell energy achieved with it: n=$n , e=$(iteration_energy_rescaled)" )
                println()
                if iseven(n)
                    j_even = Int64((n-2)/2)
                    j_odd = 1 + Int64((n-2)/2)
                else
                    j_even = Int64((n-1)/2)
                    j_odd = 1 + Int64((n-3)/2)
                end
                mask_even = [ (j<=j_even || j>=(length(omegas_even)-(j_even-1))) for j in eachindex(omegas_even) ]
                mask_odd  = [ (j<=j_odd  || j>=(length(omegas_odd) -(j_odd -1))) for j in eachindex(omegas_odd)  ]
                spectral_even[excitation_multiplet] = spectral_even[excitation_multiplet][mask_even]
                spectral_odd[excitation_multiplet]  = spectral_odd[ excitation_multiplet][mask_odd]
                omegas_even = omegas_even[mask_even]
                omegas_odd  = omegas_odd[mask_odd]
                break
            end

            for i in axes(multiplet_data,1)

                peak_energy_rescaled = multiplet_data[i,1]
                peak_weight = multiplet_data[i,2]

                contribution_positive = 0.0
                contribution_negative = 0.0
                if broadening_distribution=="gaussian"
                    if peak_energy_rescaled>0
                        contribution_positive = peak_weight*P(iteration_energy-peak_energy_rescaled*limit_iteration_scale,
                                                              spectral_broadening*iteration_scale)
                    elseif peak_energy_rescaled<0
                        contribution_negative = peak_weight*P(-iteration_energy-peak_energy_rescaled*limit_iteration_scale,
                                                              spectral_broadening*iteration_scale)
                    end
                elseif broadening_distribution=="loggaussian" && !iszero(peak_weight)
                    if peak_energy_rescaled>0
                        contribution_positive = peak_weight*P(iteration_energy-peak_energy_rescaled*limit_iteration_scale,
                                                              spectral_broadening;
                                                              distribution="loggaussian",
                                                              E=peak_energy_rescaled*limit_iteration_scale,
                                                              omega=iteration_energy)
                    elseif peak_energy_rescaled<0
                        contribution_negative = peak_weight*P(-iteration_energy-peak_energy_rescaled*limit_iteration_scale,
                                                              spectral_broadening;
                                                              distribution="loggaussian",
                                                              E=peak_energy_rescaled*limit_iteration_scale,
                                                              omega=-iteration_energy)
                    end
                end
                if isodd(n)
                    j = 1 + Int64((n-1)/2)
                    spectral_odd[excitation_multiplet][end-(j-1)] += contribution_positive
                    spectral_odd[excitation_multiplet][j] += contribution_negative
                else
                    j = Int64(n/2)
                    spectral_even[excitation_multiplet][end-(j-1)] += contribution_positive
                    spectral_even[excitation_multiplet][j] += contribution_negative
                end


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
    # compute average 
    spectral_evenodd = Dict(
        excitation_multiplet => 0.5*( spectral_even_interpolated[excitation_multiplet] + 
                                      spectral_odd_interpolated[excitation_multiplet]   )
        for excitation_multiplet in excitation_multiplets
    )

    # compute spline-interpolated smooth spectral function
    omegas_evenodd_spline,spectral_evenodd_spline = spline_interpolate_spectral_function( 
        omegas_evenodd ,
        spectral_evenodd ,
        orbitalresolved=true
   )

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

# ================= #
# NON-SIMPLE GROUPS #
# ================= #
function add_correlation_contribution_nonsimple!( 
            correlation_dict::Dict{IntMultiplet,Matrix{Float64}} ,
            A::Dict{ IntTripleG , Array{ComplexF64,4} } ,
            B::Dict{ IntTripleG , Array{ComplexF64,4} } ,
            oindex2dimensions::Vector{Int64} ,
            irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
            iteration::Int64 ,
            broadening_distribution::String ,
            spectral_broadening::Float64 ,
            iteration_scale::Float64 ,
            K_factor::Float64 ;
            correlation_type::String="spectral" ,
            T::Float64=0.0 ,
            limit_shell::Bool=false ,
            extra_iterations::Int64=0 ,
            density_matrix::Dict{IntIrrep,Matrix{Float64}}=Dict{IntIrrep,Matrix{Float64}}() ,
            multiplets_kept::Vector{IntMultiplet}=IntMultiplet[] ,
            multiplets_discarded::Vector{IntMultiplet}=IntMultiplet[] ,
            L::Float64=0.0 ,
            scale::Float64=0.0 ,
            half_weight_idx::Int64=0 ,
            half_weight_energy::Float64=0.0 
    )

    if (correlation_type=="spectral" && A==B)
        add_spectral_contribution_nonsimple!(
            correlation_dict,
            A,
            oindex2dimensions,
            irrEU,
            iteration,
            broadening_distribution,
            spectral_broadening,
            iteration_scale,
            K_factor;
            T=T,
            limit_shell=limit_shell,
            extra_iterations=extra_iterations,
            density_matrix=density_matrix,
            multiplets_kept=multiplets_kept,
            multiplets_discarded=multiplets_discarded,
            L=L,
            scale=scale,
            half_weight_idx=half_weight_idx,
            half_weight_energy=half_weight_energy
        )
    end

end
function add_spectral_contribution_nonsimple!(
            correlation_dict::Dict{IntMultiplet,Matrix{Float64}} ,
            A::Dict{ IntTripleG , Array{ComplexF64,4} } ,
            oindex2dimensions::Vector{Int64} ,
            irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
            iteration::Int64,
            broadening_distribution::String ,
            spectral_broadening::Float64 ,
            iteration_scale::Float64 ,
            K_factor::Float64 ;
            T::Float64=0.0 ,
            limit_shell::Bool=false ,
            extra_iterations::Int64=0 ,
            density_matrix::Dict{IntIrrep,Matrix{Float64}}=Dict{IntIrrep,Matrix{Float64}}() ,
            multiplets_kept::Vector{IntMultiplet}=IntMultiplet[] ,
            multiplets_discarded::Vector{IntMultiplet}=IntMultiplet[] ,
            L::Float64=0.0 ,
            scale::Float64=0.0 ,
            half_weight_idx::Int64=0 ,
            half_weight_energy::Float64=half_weight_energy )

    use_density_matrix = !iszero(length(density_matrix))

    if !use_density_matrix

        if iszero(T)
            add_spectral_contribution_T0_nonsimple!(
                correlation_dict,
                A,
                oindex2dimensions,
                irrEU,
                iteration,
                broadening_distribution,
                spectral_broadening,
                iteration_scale,
                K_factor
            )
        else
            add_spectral_contribution_Tnonzero!(
                correlation_dict,
                A,
                oindex2dimensions,
                irrEU,
                iteration,
                broadening_distribution,
                spectral_broadening,
                iteration_scale,
                K_factor,
                T,
                limit_shell,
                extra_iterations
            )
        end

    elseif use_density_matrix

        if iszero(T)
            add_spectral_contribution_density_matrix_T0!(correlation_dict,
                                                         A,
                                                         oindex2dimensions,
                                                         irrEU,
                                                         iteration,
                                                         broadening_distribution,
                                                         spectral_broadening,
                                                         iteration_scale,
                                                         K_factor,
                                                         density_matrix,
                                                         multiplets_kept,
                                                         multiplets_discarded,
                                                         L,
                                                         scale)
        else
            add_spectral_contribution_density_matrix_Tnonzero!(correlation_dict,
                                                               A,
                                                               oindex2dimensions,
                                                               irrEU,
                                                               iteration,
                                                               broadening_distribution,
                                                               spectral_broadening,
                                                               iteration_scale,
                                                               K_factor,
                                                               density_matrix,
                                                               multiplets_kept,
                                                               multiplets_discarded,
                                                               L,
                                                               scale,
                                                               half_weight_idx,
                                                               half_weight_energy,
                                                               T)

        end

    end
end
function add_spectral_contribution_T0_nonsimple!(
            correlation_dict::Dict{IntMultiplet,Matrix{Float64}} ,
            A::Dict{ IntTripleG , Array{ComplexF64,4} } ,
            oindex2dimensions::Vector{Int64} ,
            irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
            iteration::Int64 ,
            broadening_distribution::String ,
            spectral_broadening::Float64 ,
            iteration_scale::Float64 ,
            K_factor::Float64 
    )

    # ground multiplets and partition function
    G2R_0::Dict{IntIrrep,Int64} = Dict(
        G => length(filter(iszero,E))
        for (G,(E,U)) in irrEU
    )
    partition = sum( R*oindex2dimensions[I]*(S+1) for ((_,I,S),R) in G2R_0 )

    # iterate over irreducible matrices
    for ((G_u,G_a,G_v),redmat) in A

        # irrep quantum numbers
        (_,I_u,S_u) = G_u
        (_,I_a,S_a) = G_a

        energies_Gu = irrEU[G_u][1]
        energies_Gv = irrEU[G_v][1]

        # clebsch-gordan contribution 
        dim_u = oindex2dimensions[I_u]*(S_u+1)
        dim_a = oindex2dimensions[I_a]*(S_a+1)
        cg = dim_u/dim_a

        # iterate over creation operator multiplets
        for r_a in axes(redmat,3)

            m_a = (G_a...,r_a)
            correlation_matrix = correlation_dict[m_a]

            # iterate over eigenstate multiplets
            for r_u in axes(redmat,2),
                r_v in axes(redmat,4)

                # excitation energy
                e_u = energies_Gu[r_u]
                e_v = energies_Gv[r_v]
                e_diff = e_u-e_v

                # iterate over cg multiplicity index
                for α in axes(redmat,1)

                    # coefficient from reduced matrix element
                    redmatel = abs2(redmat[α,r_u,r_a,r_v])

                    # total contribution 
                    w = cg*redmatel/partition

                    # particle contribution: m_v ground, positive peak
                    if iszero(e_v) && e_diff>0
                        #contribution = w*P((K_factor-e_diff)*iteration_scale,
                        #                                           (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*iteration_scale);
                        #                                           distribution=broadening_distribution,
                        #                                           E=e_diff*iteration_scale,
                        #                                           omega=K_factor*iteration_scale)
                        correlation_matrix[end-iteration,2] += w*P((K_factor-e_diff)*iteration_scale,
                                                                   (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*iteration_scale);
                                                                   distribution=broadening_distribution,
                                                                   E=e_diff*iteration_scale,
                                                                   omega=K_factor*iteration_scale)
                    # hole contribution: m_u ground, negative peak
                    elseif iszero(e_u) && e_diff<0
                        correlation_matrix[iteration+1,2] += w*P((-K_factor-e_diff)*iteration_scale,
                                                                 (broadening_distribution=="loggaussian" ? spectral_broadening : spectral_broadening*iteration_scale);
                                                                 distribution=broadening_distribution,
                                                                 E=e_diff*iteration_scale,
                                                                 omega=K_factor*iteration_scale)
                    end
                end
            end
        end
    end

end

# update excitation operator
function update_operator_nonsimple(
        redmat_iaj::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,4} } ,
        multiplets_a_block::Vector{NTuple{4,Int64}} ,
        fsum::FSum ,
        combinations_Gu_muiualpha::Dict{NTuple{3,Int64}, Vector{Tuple{IntMultiplet,IntMultiplet,IntMultiplet,Int64}}} ,
        irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
        cg_o_comb2A::Dict{ NTuple{3,Int64} , Int64 }
    )::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,4} }

    # G => multiplicity of G
    G2R_uv::Dict{IntIrrep,Int64} = Dict( G_u=>size(U_u,1) for (G_u,(_,U_u)) in irrEU )
    # a
    G2R_a::Dict{NTuple{3,Int64},Int64} = Dict( 
        G=>R for (G,R) in get_irreps( Set(multiplets_a_block) ; multiplicity=true ) 
    )

    # irrep -> decomposition -> multiplet
    Gu2GmuGi2rmuiualpha::Dict{IntIrrep,Dict{NTuple{2,IntIrrep},Vector{NTuple{4,Int64}}}} = Dict(
        G_u => Dict(
        (G_mu,G_i) => [(m_mu[4],m_i[4],m_u[4],α_u) for (m_mu,m_i,m_u,α_u) in multiplet_combinations_Gu if (m_i[1:3]==G_i && m_mu[1:3]==G_mu)]
            for (G_mu,G_i) in Set( (m_mu[1:3],m_i[1:3]) for (m_mu,m_i,_,_) in multiplet_combinations_Gu )
        )
        for (G_u,multiplet_combinations_Gu) in combinations_Gu_muiualpha
    )
    # Gu2GmuGi2rmuiualpha = (
    #     G_u::IntIrrep => (
    #     (G_mu,G_i)::NTuple{2,IntIrrep} => [(m_mu[4],m_i[4],m_u[4],α_u)::NTuple{4,Int64} for (m_mu,m_i,m_u,α_u) in multiplet_combinations_Gu if (m_i[1:3]==G_i && m_mu[1:3]==G_mu)]
    #             for (G_mu,G_i)::NTuple{2,IntIrrep} in Set( (m_mu[1:3],m_i[1:3]) for (m_mu,m_i,_,_) in multiplet_combinations_Gu )::Set{NTuple{2,IntIrrep}}
    #     )
    #     for (G_u,multiplet_combinations_Gu) in combinations_Gu_muiualpha
    # )

    # full-sized matrices for avoiding allocations
    R_uv_max = maximum(values(G2R_uv))
    R_a_max = maximum(values(G2R_a))
    A_max = maximum(values(cg_o_comb2A))
    tmp_full = zeros(ComplexF64,R_uv_max,R_uv_max)
    uav_matrix_full = zeros(ComplexF64,A_max,R_uv_max,R_a_max,R_uv_max)

    # initialize result dictionary
    redmat_uav::Dict{ IntTripleG , Array{ComplexF64,4} } = Dict()

    # G_u, G_v iteration
    for (G_u,GmuGi2rmuiualpha) in Gu2GmuGi2rmuiualpha,
        (G_v,GnuGj2rnujvalpha) in Gu2GmuGi2rmuiualpha

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

            # maximum outer multiplicity β=1,…,B_uav
            B_uav = get( cg_o_comb2A , (I_u,I_a,I_v) , zero(Int64) )
            iszero(B_uav) && continue

            # < Γ_u || f†_{Γ_a} || Γ_v >_β
            # @views uav_matrix = uav_matrix_full[1:B_uav,1:R_u,1:R_a,1:R_v]
            @views uav_matrix = uav_matrix_full[1:B_uav,1:R_u,1:R_a,1:R_v]
            uav_matrix .= zero(ComplexF64)

            # G_i,G_mu,G_j,G_nu iteration
            for ((G_mu,G_i),rmuiualphas) in GmuGi2rmuiualpha,
                ((G_nu,G_j),rnujvalphas) in GnuGj2rnujvalpha

                # irrep quantum numbers 
                N_i,I_i,S_i = G_i
                N_j,I_j,S_j = G_j
                N_munu,I_munu,S_munu = G_mu

                # early discard
                G_mu==G_nu     || continue
                N_i==(N_j+N_a) || continue
                G_munu = G_mu

                # clebsch-gordan sum array
                f_spin_and_sign,f_orbital_array = fsum[(G_u,G_v,G_a,G_i,G_j,G_munu)]
                (iszero(f_spin_and_sign) && iszero(length(f_orbital_array))) && continue

                # ⟨ Γ_i || f†_{Γ_a} || Γ_j ⟩_α
                iaj_matrix::Array{ComplexF64,4} = get( redmat_iaj , (G_i,G_a,G_j) , zeros(ComplexF64,0,0,0,0) )
                iszero(length(iaj_matrix)) && continue

                # maximum outer multiplicity α=1,…,A_iaj
                A_iaj = cg_o_comb2A[I_i,I_a,I_j]

                # u,i,mu,α_u ; v,j,nu,α_v  multiplet iteration
                for (r_mu,r_i,r_u,α_u) in rmuiualphas,
                    (r_nu,r_j,r_v,α_v) in rnujvalphas

                    r_mu==r_nu || continue

                    # outer multiplicites of ⟨i||a||j⟩^[n-1]_α, ⟨u||a||v⟩^[n]_β
                    for α in 1:A_iaj,
                        β in 1:B_uav

                        f = f_spin_and_sign * f_orbital_array[α_u,α_v,α,β]

                        @inbounds for r_a in 1:R_a
                            uav_matrix[β,r_u,r_a,r_v] += f * iaj_matrix[α,r_i,r_a,r_j]
                        end
                    end

                end # u,i,mu , v,j,nu , a multiplet iteration

            end # G_i,G_mu,G_j,G_nu iteration

            # final discard
            is_matrix_zero(uav_matrix) && continue

            # transform matrix
            @inbounds @views for r_a in 1:G2R_a[G_a],
                                 β   in 1:B_uav

                mul!( tmp , uav_matrix[β,:,r_a,:] , U_v )
                mul!( uav_matrix[β,:,r_a,:] , U_u' , tmp )

            end

            # redmat_uav
            # @inbounds push!( redmat_uav , (G_u,G_a,G_v) => uav_matrix[:,:,:,:] )
            res::Array{ComplexF64,4} = copy(uav_matrix)
            @inbounds push!( redmat_uav , (G_u,G_a,G_v) => res )
            #@inbounds redmat_uav[(G_u,G_a,G_v)] = uav_matrix_full[1:R_u,1:R_a,1:R_v]
            #@inbounds redmat_uav[(G_u,G_a,G_v)] = uav_matrix[:,:,:]

        end # G_a iteration
    end # G_u, G_v iteration

    return redmat_uav
end
