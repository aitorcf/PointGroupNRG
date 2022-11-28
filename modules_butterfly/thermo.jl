# ##########################################################
# This module contains functions that compute thermodynamic
# averages and other related quantities with the information 
# of a given NRG step. 
# ##########################################################

include( "discretization.jl" )

function temperature( N::Int64 , L::Float64 , betabar::Float64 ; z::Float64=0.0 , discretization="standard" ) 
    if discretization=="standard"
        return 0.5 * (1+L^(-1)) * L^(z-(N-2)/2) / (betabar)#*(L^(-z)) # compensate for the L^(z) on betabar 
    elseif discretization=="co2005" 
        eco0_z = compute_epsilon_co2005(8,z,L)[9]*L^4
        return eco0_z*L^(N-2)/betabar 
    end
end

function partition( 
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            betabar::Float64 ,
            oirreps2dimensions::Dict )

    #   Z = Tr exp( - betabar * H_N )
    #     = Sum_G( D * ( Sum_i exp(-betabar*E_i) ))

    part = 0
    for (G,(E,U)) in irrEU
        # orbital and spin
        (I,S) = G[2:3] 

        # degeneracy/dimensionality
        Do = oirreps2dimensions[I]
        Ds = 2*S+1
        D = Do*Ds
        
        # weighted sum over G 
        part += D * sum( exp( - betabar * e ) for e in E )
    end
    return part
end

# int method
function partition( 
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            betabar::Float64 ,
            oindex2dimensions::Vector{Int64} )
    #   Z = Tr exp( - betabar * H_N )
    #     = Sum_G( D * ( Sum_i exp(-betabar*E_i) ))

    part = 0
    for (G,(E,U)) in irrEU
        # orbital and spin
        (I,S) = G[2:3] 

        # degeneracy/dimensionality
        Do = oindex2dimensions[I]::Int64
        Ds = S+1
        D = Do*Ds
        
        # weighted sum over G 
        part += D * sum( exp( - betabar * e ) for e in E )
    end
    return part::Float64
end

function magsusc( irrEU , betabar , oirreps2dimensions::Dict ; verbose=false )

    verbose && println( "MAGNETIC SUSCEPTIBILITY CALCULATION" )

    part = partition( irrEU , betabar , oirreps2dimensions )

    mag = 0
    for (G,(E,U)) in irrEU
        # orbital and spin
        (I,S) = G[2:3] 

        # degeneracy/dimensionality
        Do = oirreps2dimensions[I]
        Ds = 2*S+1
        D = Do*Ds
        
        # weighted sum over G 
        # S^2 = S(S+1)
        # Sz^2 = S^2/3
        contrib = D * (S*(S+1)/3.0) * sum( exp( - betabar * e ) for e in E ) 
        mag += contrib

        if verbose 
            println( "G = $G, D = $D, E = $E" )
            println( "contribution = $(contrib/part)" )
            println()
        end
    end
    return mag/part 
end

# int method
function magsusc( 
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            betabar::Float64 ,
            oindex2dimensions::Vector{Int64} ;
            verbose=false )

    verbose && println( "MAGNETIC SUSCEPTIBILITY CALCULATION" )

    part::Float64 = partition( irrEU , betabar , oindex2dimensions )

    mag::Float64 = 0
    for (G,(E,U)) in irrEU
        # orbital and spin
        (I,S) = G[2:3] 

        # degeneracy/dimensionality
        Do = oindex2dimensions[I]::Int64
        Ds = S+1
        D = Do*Ds
        
        # weighted sum over G 
        # S^2 = S(S+1)
        # Sz^2 = S^2/3
        contrib = D * (S/2.0*(S/2.0+1)/3.0) * sum( exp( - betabar * e ) for e in E ) 
        mag += contrib

        if verbose 
            println( "G = $G, D = $D, E = $E" )
            println( "contribution = $(contrib/part)" )
            println()
        end
    end
    return mag/part 
end

# string method
function number( irrEU , betabar , oirreps2dimensions::Dict ; verbose=false )

    part = partition( irrEU , betabar , oirreps2dimensions )

    N = 0
    for (G,(E,U)) in irrEU
        # orbital and spin
        (N,I,S) = G

        # degeneracy/dimensionality
        Do = oirreps2dimensions[I]
        Ds = 2*S+1
        D = Do*Ds
        
        # weighted sum over G 
        contrib = D * N * sum( exp( - betabar * e ) for e in E ) 
        N += contrib
    end
    return N/part 
end

# int method 
function number( 
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            betabar::Float64 ,
            oindex2dimensions::Vector{Int64} ;
            verbose=false )

    part::Float64 = partition( irrEU , betabar , oindex2dimensions )

    N::Float64 = 0.0
    for (G,(E,U)) in irrEU
        # orbital and spin
        (N,I,S) = G

        # degeneracy/dimensionality
        Do = oindex2dimensions[I]
        Ds = S+1
        D = Do*Ds
        
        # weighted sum over G 
        contrib = D * N * sum( exp( - betabar * e ) for e in E ) 
        N += contrib
    end
    return N/part 
end

function energy( irrEU , betabar , oirreps2dimensions::Dict ; verbose=false )

    part = partition( irrEU , betabar , oirreps2dimensions )

    energy = 0
    for (G,(E,U)) in irrEU
        # orbital and spin
        (N,I,S) = G

        # degeneracy/dimensionality
        Do = oirreps2dimensions[I]
        Ds = 2*S+1
        D = Do*Ds
        
        # weighted sum over G 
        contrib = D * sum( (e*exp( - betabar * e )) for e in E ) 
        energy += contrib
    end
    return energy/part 
end

# int method
function energy( 
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            betabar::Float64 ,
            oindex2dimensions::Vector{Int64} ;
            verbose=false )

    part::Float64 = partition( irrEU , betabar , oindex2dimensions )

    energy::Float64 = 0
    for (G,(E,U)) in irrEU
        # orbital and spin
        (N,I,S) = G

        # degeneracy/dimensionality
        Do = oindex2dimensions[I]
        Ds = S+1
        D = Do*Ds
        
        # weighted sum over G 
        contrib = D * sum( (e*exp( - betabar * e )) for e in E ) 
        energy += contrib
    end
    return energy/part 
end

# int method 
function entropy( 
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            betabar::Float64 ,
            oindex2dimensions::Vector{Int64} ;
            verbose=false )

    part::Float64 = partition( irrEU , betabar , oindex2dimensions ) 
    return log(part)

end

# ----------------------- #
# IMPURITY THERMODYNAMICS #
# ----------------------- #
function imp_qnums( 
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            oindex2dimensions::Vector{Int64} ,
            combinations_uprima::Dict{NTuple{3,Int64}, Vector{NTuple{3,NTuple{4,Int64}}}},
            ss_ip ,
            nn_ip )
    # compute averages of impurity quantum numbers
    # in order to keep track of the atomic low-energy
    # states

    # average impurity spin^2 and N per multiplet
    ss_u = Dict()
    nn_u = Dict()

    # iterate over irrep blocks
    for (G,(E,U)) in irrEU 

        # iterate over multiplets in irrep block
        # (store the block-shell combination)
        for (m_u,m_mu,m_i) in combinations_uprima[G]

            # qnum setup
            ss_u[m_u] = 0.0
            nn_u[m_u] = 0.0

            # outer multiplicity
            r_u = m_u[4]

            ## iterate over primed multiplets
            for c_up in combinations_uprima[G]

                # U(Γ) term
                m_up = c_up[1]
                r_up = m_up[4]
                uterm = abs2(U[r_up,r_u])

                # primed term 
                m_ip = c_up[3]

                # quantum number averages
                ss_u[m_u] += uterm*ss_ip[m_ip]
                nn_u[m_u] += uterm*nn_ip[m_ip]

            end
        end
    end

    return (ss_u,nn_u)
end

function impspin( 
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            betabar::Float64 ,
            oindex2dimensions::Vector{Int64} ,
            ss_u ; 
            verbose=false )

    s_imp = 0.0

    part::Float64 = partition( irrEU , betabar , oindex2dimensions )

    energy = 0
    for (G,(E,U)) in irrEU

        # number of multiplets in Γ
        R_G = length(E)

        # orbital and spin
        (N,I,S) = G

        # degeneracy/dimensionality
        Do = oindex2dimensions[I]
        Ds = S+1
        D = Do*Ds
        
        # weighted sum over Γ
        for r_g in 1:R_G 
            m_u = (G...,r_g)
            s_imp += D*ss_u[m_u]*exp(-betabar*E[r_g])
        end
        
    end
    return s_imp/part
end

function impnum( 
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            betabar::Float64 ,
            oindex2dimensions::Vector{Int64} ,
            nn_u ; 
            verbose=false )

    n_imp = 0.0

    part::Float64 = partition( irrEU , betabar , oindex2dimensions )

    energy = 0
    for (G,(E,U)) in irrEU

        # number of multiplets in Γ
        R_G = length(E)

        # orbital and spin
        (N,I,S) = G

        # degeneracy/dimensionality
        Do = oindex2dimensions[I]
        Ds = S+1
        D = Do*Ds
        
        # weighted sum over Γ
        for r_g in 1:R_G 
            m_u = (G...,r_g)
            n_imp += D * nn_u[m_u] * exp(-betabar*E[r_g]) 
        end
        
    end
    return n_imp/part
end

# ------------------- #
# IMPURITY MULTIPLETS #
# ------------------- #
function imp_mults( 
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            oindex2dimensions::Vector{Int64} ,
            combinations_uprima::Dict{NTuple{3,Int64}, Vector{NTuple{3,NTuple{4,Int64}}}},
            mm_ip::Dict{NTuple{4,Int64},Vector{Float64}} )

    # average impurity spin^2 and N per multiplet
    M = length(collect(values(mm_ip))[1])
    mm_u::Dict{NTuple{4,Int64},Vector{Float64}} = Dict()

    # iterate over irrep blocks
    for (G,(E,U)) in irrEU 

        # iterate over multiplets in irrep block
        # (store the block-shell combination)
        for (m_u,m_mu,m_i) in combinations_uprima[G]

            # qnum setup
            mm_u[m_u] = Vector{Float64}([0.0 for i in 1:M])

            # outer multiplicity
            r_u = m_u[4]

            ## iterate over primed multiplets
            for c_up in combinations_uprima[G]

                # U(Γ) term
                m_up = c_up[1]
                r_up = m_up[4]
                uterm = abs2(U[r_up,r_u])

                # primed term 
                m_ip = c_up[3]

                # quantum number averages
                mm_u[m_u] .+= uterm*mm_ip[m_ip]

            end
        end
    end

    return mm_u
end
function mult_thermo( 
            irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
            betabar::Float64 ,
            oindex2dimensions::Vector{Int64} ,
            mm_u::Dict{NTuple{4,Int64},Vector{Float64}} ; 
            verbose=false )

    M::Int64 = length(collect(values(mm_u))[1])
    mult_imp::Vector{Float64} = [0.0 for i in 1:M]

    part::Float64 = partition( irrEU , betabar , oindex2dimensions )

    for (G,(E,U)) in irrEU

        # number of multiplets in Γ
        R_G = length(E)

        # orbital and spin
        (N,I,S) = G

        # degeneracy/dimensionality
        Do = oindex2dimensions[I]
        Ds = S+1
        D = Do*Ds
        
        # weighted sum over Γ
        for r_g in 1:R_G 
            m_u = (G...,r_g)
            mult_imp .+= D*mm_u[m_u]*exp(-betabar*E[r_g])
        end
        
    end
    return mult_imp./part
end
