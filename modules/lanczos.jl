#!/usr/bin/env julia 

using LinearAlgebra 

function lanczos( H::Matrix{BigFloat} , v_1::Vector{BigFloat} ) 

    M = size(H)[1]

    U::Matrix{BigFloat} = similar(H)
    U[:,1] .= v_1

    vv::Vector{Vector{BigFloat}} = [v_1]

    # set first seed vector
    v_n = v_1
    # iterate
    for i=1:(M-1)

        # produce non-normalized vector 
        v = H*v_n - sum( (vv[j]'*H*v_n)*vv[j] for j=1:i )

        # normalize vector 
        v ./= sqrt( v'*v )

        # update seed vector
        v_n = v

        # introduce new vector in list
        push!( vv , v )

        # introduce new vector in U
        U[:,i+1] .= v

    end

    return U' * H * U

end

function integrate( f::Function , 
                    interval::Tuple{BigFloat,BigFloat} ;
                    Npoints::Int64=Int64(1e4) )::BigFloat


    d::BigFloat = interval[2]-interval[1]
    dx::BigFloat = d/Npoints

    # xpoints 
    xx = BigFloat[ interval[1]+n*dx for n in 0:(Npoints-1) ]
     
    integral::BigFloat = 0.0 
    average::BigFloat = 0.0
    @inbounds for i in 1:(length(xx)-1)

        average = (f(xx[i])+f(xx[i+1]))/2.0

        integral += average*dx
        
    end

    return integral

end

function representative_energy( interval ; scheme="alternative" , eta::Function=x->1.0 ) 

    if scheme=="alternative" 

        return integrate( x->eta(x)   , interval ) /
               integrate( x->eta(x)/x , interval )

    elseif scheme=="standard" 

        return integrate( x->x   , interval ) /
               integrate( x->1.0 , interval )

    end

end

function get_f1_E( eta::Function , 
                   L::BigFloat , 
                   z::BigFloat , 
                   N::Int64 ;
                   scheme="alternative" )

    # diagonal conduction matrix elements
    E  = BigFloat[ 0 for n in 1:(2*N) ]

    # coefficients of f_0 expressed in the a_n basis
    f1 = BigFloat[ 0 for n in 1:(2*N) ]

    # etabar: average of eta
    etabar::BigFloat = integrate( x->eta(x) , 
                        (BigFloat(-1.0),BigFloat(1.0)) ; 
                        Npoints=Int64(1e6) )

    # type-define variables 
    idx::Int64 = 0
    x_pn::BigFloat = 0.0
    x_pn1::BigFloat = 0.0
    interval::NTuple{2,BigFloat} = (0.0,0.0)

    for n=1:N

        # negative part
        idx = n
        x_pn  = n==1 ? -1.0 : -L^(-n+1+z) 
        x_pn1 = -L^(-n+z)
        interval = ( x_pn , x_pn1 )
        E[idx] = representative_energy( interval ; scheme=scheme , eta=eta )
        f1[idx] = sqrt( integrate(eta,interval) / etabar )
        
        # positive part
        idx = 2*N - (n-1)
        x_pn  = n==1 ? 1.0 : L^(-n+1+z) 
        x_pn1 = L^(-n+z)
        interval = ( x_pn1 , x_pn )
        E[idx] = representative_energy( interval ; scheme=scheme , eta=eta )
        f1[idx] = sqrt( integrate(eta,interval) / etabar )

    end

    return (f1,E)

end

function get_hoppings( N::Int64 , 
                       L::Float64 , 
                       z::Float64 , 
                       eta::Function ;
                       scheme="alternative" )

    # convert Float64 to BigFloat
    LL = BigFloat(L) 
    zz = BigFloat(z)

    # obtain innermost shell and diagonal basis elements
    (f1,E) = get_f1_E( eta , LL , zz ,  N ; scheme=scheme )

    # construct diagonal hamiltonian 
    H::Matrix{BigFloat} = diagm( E )

    # construct off-diagonal hamiltonian
    Hp = lanczos( H , f1 )

    # off-diagonal entries ϵ
    hops = [Hp[i,i+1] for i in 1:N]

    ## diagonal entries 
    #diags = [Hp[i,i] for i in 1:N]
    #@show diags

    # asymptotic off-diagonals ̄ϵ
    hopscaling = [abs(hops[i]*sqrt(L)-hops[i+1]) for i in 7:12]
    C = findfirst( x->x==minimum(hopscaling) , hopscaling )+6
    C = 7
    @show C
    hop0 = hops[C]*LL^((C-1)/2)
    hops_asymptotic = [hop0*LL^(-(n-1)/2) for n in 1:N] 

    # fractions ξ
    xis = hops ./ hops_asymptotic
    xis[C:end] .= 1.0

    return (Float64.(hops),Float64.(hops_asymptotic),Float64.(xis))

end

test = false
if test
    # -----------------------
    # TESTING 

    N = 40
    L = 10.0
    eta = x->1.0
    scheme = "standard" 

    z = 0.0
    println( "======" )
    println( "z = $z" )
    println( "======" )
    h,ha,xi = get_hoppings(N,L,z,eta;scheme=scheme) 
    @show h 
    println()
    @show ha 
    println()
    @show xi
    println()
    analytic = [(n>1 ? (1-L^(-(n-2)-1)) * (1-L^(-2*(n-2)-1))^(-0.5) * (1-L^(-2*(n-2)-3))^(-0.5) : 0.0) for n in 1:N]
    @show analytic
    println()


    z = 0.5
    println( "======" )
    println( "z = $z" )
    println( "======" )
    h2,ha2,xi2 = get_hoppings(N,L,z,eta;scheme=scheme ) 
    @show h2 
    println()
    @show ha2 
    println()
    @show xi2
    println()
    @show analytic

    println( "==========" )
    println( "COMPARISON" )
    println( "==========" )
    @show ha[1]*L^z 
    @show ha2[1]
    analytic = [(n>1 ? (1-L^(-(n-2)-1)) * (1-L^(-2*(n-2)-1))^(-0.5) * (1-L^(-2*(n-2)-3))^(-0.5) : 0.0) for n in 1:N]



    # -----------------------

end
