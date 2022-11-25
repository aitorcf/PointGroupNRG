# ===============================================
#
# This module contains functions that can be used 
# in order to compute the hopping coefficients of 
# an Anderson Hamiltonian. 
#
# For the moment, IT ASSUMES THAT THE DISPERSION
# IS LINEAR, of the form ϵ_k = k.
#
# The calculation of ϵ and ξ coefficients is done 
# using the numeric formula until N=7, and then 
# the asymptotic formula for N≥8.
#
# ===============================================


# computes the values of z given N_z
get_Z( Nz::Int64 ) = collect( Float64(i)/Nz for i=0:(Nz-1) )

# inter-shell hopping as in KWW (1980)
function xi( n::Int64 , L::Float64 ; z::Float64=0.0 ) 
    nn = n - 2 
    return n>1 ? (1-L^(-nn-1)) * (1-L^(-2*nn-1))^(-0.5) * (1-L^(-2*nn-3))^(-0.5)::Float64 : Float64(0.0)
end

# computes ξ_n coefficients, where 
#       ξ_n = ϵ_n / ebar_n,
# ebar_n being the values obtained using the
# asymptotic formula.
function compute_xi_vector( N::Int64 , 
                            z::Float64 , 
                            L::Float64 ; 
                            discretization::String="standard" ,
                            verbose::Bool=false ) 

    zz = BigFloat(z)
    LL = BigFloat(L)
    N7 = minimum((7,N))

    ϵ = compute_epsilon_vector(N7,
                               zz,
                               LL;
                               verbose=verbose,
                               discretization=discretization)

    if discretization=="standard"

        return vcat( [ϵ[n]/ebar(z,L,n-1) for n=1:N7] , 
                     [1 for _=(N7+1):N] )

    elseif discretization=="co2005" 

        eco0 = compute_epsilon_co2005(8,z,L)[9]*L^4
        ebars = [eco0*L^(-n/2) for n=0:N]

        return vcat( [ϵ[i]/ebars[i] for i=1:N7] , 
                     [1 for _=(N7+1):N] )

    end
end

function compute_epsilon_co2005( 
            N::Int64 , 
            z::R ,
            L::R ; 
            verbose=false ) where {R<:Real}

    # ϵ_i = ϵ_{i+1} for there not to be a 0 index
    epsilon::Vector{BigFloat} = zeros( BigFloat , N+1 )

    # convert to bigfloat
    z,L = BigFloat.((z,L))

    prefactor::BigFloat = 1.0
    e::BigFloat = 0.0

    N7 = minimum((N,7))
    
    # numeric part 
    for M=0:N

        if verbose
            println( "Computing ϵ_$M" )
            println( "==============" )
        end

        f = F(M,z,L ; discretization="co2005" )
        h = H11(M+1,epsilon[1:M]) 
        e = prefactor*sqrt( f - h )

        if verbose 
            @show M
            @show f 
            @show h
            @show prefactor
            @show e
            println()
        end

        prefactor /= e
        epsilon[M+1] = e

    end
    return epsilon

end

# computes the ϵ_n coefficients
function compute_epsilon_vector( 
                N::Int64 , 
                z::BigFloat , 
                L::BigFloat ; 
                discretization="standard" ,
                verbose=true ) 

    # ϵ_i -> ϵ_{i+1} for there not to be a 0 index
    epsilon::Vector{BigFloat} = zeros( BigFloat , N+1 )

    prefactor::BigFloat = 1.0
    e::BigFloat = 0.0

    N7 = minimum((N,7))
    
    if discretization=="standard"

        # numeric part 
        for M=0:N7
            if verbose
                println( "Computing ϵ_$M" )
                println( "==============" )
            end
            f = F(M,z,L ; discretization=discretization )
            h = H11(M+1,epsilon[1:M]) 
            e = prefactor*sqrt( f - h )
            if verbose 
                @show M
                @show f 
                @show h
                @show prefactor
                @show e
                println()
            end
            prefactor /= e
            epsilon[M+1] = e
        end

        # asymptotic part
        epsilon[N7+1:N+1] = collect( ebar(z,L,n-1) for n in (N7+1):(N+1) ) 

    elseif discretization=="co2005"
        # numeric part 
        for M=0:N
            if verbose
                println( "Computing ϵ_$M" )
                println( "==============" )
            end
            f = F(M,z,L ; discretization=discretization )
            h = H11(M+1,epsilon[1:M]) 
            e = prefactor*sqrt( f - h )
            if verbose 
                @show M
                @show f 
                @show h
                @show prefactor
                @show e
                println()
            end
            prefactor /= e
            epsilon[M+1] = e
        end

        # asymptotic part
        #epsilon[N7+1:N+1] = collect( ebar(z,L,n-1) for n in (N7+1):(N+1) ) 
    end

    return epsilon

end

# m-th k-point in a z-displaced discretization
function k_m( z::BigFloat , L::BigFloat , m::Int64 )
    return m==0 ? BigFloat(1.0) : BigFloat(L^(-m+z))
end

# distance between adjacent k-points 
# k_m and k_{m+1}
function d( z::BigFloat , L::BigFloat , m::Int64 ) 
    return BigFloat( k_m(z,L,m) - k_m(z,L,m+1) )
end

# F_M(z,Λ) for ϵ_k = k
function F( N::Int64 , z::BigFloat , L::BigFloat ; discretization="standard" )
    N2::Int64 = 2*N+2 
    N3::Int64 = 2*N+3
    f0::BigFloat = 0.0
    fr::BigFloat = 0.0
    # original discretization, the one appearing in yoshida, whitaker, oliveira
    if discretization=="standard"
        Lm1::BigFloat = L^(-1) 
        Lm1z::BigFloat = L^(-1+z)
        f0 = (1-Lm1z) * ((1+Lm1z)/2.0)^N2
        fr = (1-Lm1) * ((1+Lm1)/2.0)^N2 * Lm1z^N3 / (1-L^(-N3))
        return f0 + fr
    # alternative discretization value from campo & oliveira 2005
    elseif discretization=="co2005" 
        f0 = (1-L^(z-1))^N3 / ((1-z)*log(L))^N2 
        fr = (1-L^(-1))^N3 * L^(N3*(z-1)) / (log(L))^N2 / (1-L^(-N3)) 
        return f0 + fr
    end
end

# E_M(z,Λ): energy average in the interval
function E( z::BigFloat , L::BigFloat , m::Int64 ; discretization="standard" ) 
    if discretization=="standard"
        k_h::BigFloat = k_m(z,L,m+1)
        k_a::BigFloat = k_m(z,L,m)
        return BigFloat( (k_a^2-k_h^2) / 2.0 / d(z,L,m) )
    elseif discretization=="co2005"
        return E_CO2005(z,L,m)
    end
end

# E_M(z,Λ) as defined in the alternative discretization 
# proposed by Campo & Oliveira (2005)
function E_CO2005( z::BigFloat , L::BigFloat , m::Int64 )
    if m==0
        return BigFloat( (1-L^(z-1))/((1-z)*log(L)) )
    elseif m>0 
        return BigFloat( (1-L^(-1))*L^(z-m)/log(L) )
    end
end

# [(H_N)^(2N)]_11
function H11( M1::Int64 , epsilon::Vector{BigFloat} )

    M1==1 && return 0.0

    H::Matrix{BigFloat} = zeros(BigFloat,M1,M1)
    # assign terms (first and last aside)
    H[1,2] = epsilon[1]
    H[M1,M1-1] = epsilon[M1-1]
    @inbounds for i=2:(M1-1)
        H[i,i+1] = epsilon[i]
        H[i,i-1] = epsilon[i-1]
    end
    # take the power (2N+2)=2M1
    Hpow::Matrix{BigFloat} = H^(2*M1)
    return Hpow[1,1]
end

function ebar( z::R , L::R , n::Int64 ) where {R<:Real}
    return convert( R , 0.5*(1+L^(-1))*L^(z-n/2.0) )
end

function compute_ebar0_z( z::R , L::R ; discretization="standard" ) where {R<:Real} 
    if discretization=="standard" 
        return ebar(z,L,0)
    elseif discretization=="co2005" 
        return compute_epsilon_co2005(8,z,L)[9]*(L^4)
    end
end
