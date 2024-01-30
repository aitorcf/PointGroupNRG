#!/usr/bin/env julia 


# IMPLEMENTED DISCRETIZATION METHODS
#
# Discretization scheme
# - yoshida1990
# - campo2005
#
# Tridiagonalization 
# - yoshida1990
# - chen1995

using LinearAlgebra 
using Printf

const discretization_default = "campo2005"
const tridiagonalization_default = "chen1995"
const enforce_particle_hole_symmetry_default = true
const enforce_asymptotic_behavior_default = true

const discretizations = Set([
    "yoshida1990" , 
    "campo2005" 
])
const tridiagonalizations = Set([
    "matrix" ,
    "chen1995" ,
    "gonzalez-buxton1998" # brute-force chen1995
])
const integrate_methods = Set([
    "left" ,
    "middle"
])

# ========= #
# Interface #
# ========= #
function discretize_bands( dos_channels::Dict{String,Vector{Function}} ,
                           L::Float64 ,
                           z::Float64 ,
                           number_of_shells::Int64 ;
                           discretization::String=discretization_default ,
                           tridiagonalization::String=tridiagonalization_default ,
                           enforce_particle_hole_symmetry::Bool=enforce_particle_hole_symmetry_default ,
                           enforce_asymptotic_behavior::Bool=enforce_particle_hole_symmetry_default )

    @assert (discretization in discretizations) "Discretization not implemented."
    @assert (tridiagonalization in tridiagonalizations) "Tridiagonalization not implemented."

    # obtain diagonal and codiagonal elements of the tridiagonal hamiltonians
    println( "Computing discretized channels..." )
    channels_tridiagonal::Dict{String,Vector{Tuple{Vector{Float64},Vector{Float64}}}} =
        Dict( 
            orbital_irrep => [
                discretize_band( dos,
                                 L,
                                 z,
                                 number_of_shells;
                                 discretization=discretization,
                                 tridiagonalization=tridiagonalization ,
                                 enforce_particle_hole_symmetry=enforce_particle_hole_symmetry ,
                                 enforce_asymptotic_behavior=enforce_particle_hole_symmetry_default )
                for dos in dos_functions
            ] for (orbital_irrep,dos_functions) in dos_channels
        )

    return channels_tridiagonal

end

# ============= #
# Main function #
# ============= #

# Hybridization as Γ(ϵ)
function discretize_band( 
            dos::Function ,
            L::Float64 ,
            z::Float64 ,
            number_of_shells::Int64 ;
            discretization::String=discretization_default ,
            tridiagonalization::String=tridiagonalization_default ,
            enforce_particle_hole_symmetry::Bool=false ,
            enforce_asymptotic_behavior::Bool=enforce_particle_hole_symmetry_default ,
            verbose::Bool=false ,
            check::Bool=false )

    if verbose
        println( "BAND DISCRETIZATION" )
        println( "L = $L" )
        println( "z = $z" )
        println( "particle-hole symmetry enforced? $enforce_particle_hole_symmetry" )
    end

    number_of_shells_needed = number_of_shells
    number_of_shells_to_compute = number_of_shells<100 ? 100 : number_of_shells

    # maximum discretization interval to be considered
    maximum_j = iseven(number_of_shells_to_compute) ? number_of_shells_to_compute÷2 : number_of_shells_to_compute÷2+1

    # convert all parameters to BigFloat
    ddooss(x) = BigFloat(dos(x))
    LL = BigFloat(L)
    zz = BigFloat(z)

    # hamiltonian in diagonal basis according to discretization
    Hcon_diagonal = construct_nontransformed_hamiltonian(
                        ddooss ,
                        LL ,
                        maximum_j ,
                        zz ,
                        discretization )

    # innermost shell vector
    f0 = innermost_shell_vector( ddooss ,
                                 LL ,
                                 maximum_j ,
                                 zz )

    # lanczos tridiagonalization according to chosen method
    is_dos_symmetric = all( isapprox(dos(x),dos(-x)) for x in 0.0:0.01:1.0 )
    if (enforce_particle_hole_symmetry && !is_dos_symmetric && verbose) 
        println( "DOS not symmetric: cannot enforce p-h symmetry." )
    end
    # enforce_ph: if DOS is ph-symmetric, diagonals are identically 0
    enforce_ph = enforce_particle_hole_symmetry && is_dos_symmetric
    diagonals,codiagonals = lanczos( Hcon_diagonal , f0 , LL ; method=tridiagonalization , enforce_particle_hole_symmetry=enforce_ph )

    # check
    check && check_tridiagonalization( diagonals , 
                                       codiagonals , 
                                       Hcon_diagonal ; 
                                       enforce_asymptotic_behavior=enforce_asymptotic_behavior ,
                                       L=LL )

    # rescaled parameters
    #diagonals_rescaled   = Float64.(diagonals.*collect( sqrt(L)^(m-2) for m in eachindex(diagonals) ))
    #codiagonals_rescaled = Float64.(codiagonals.*collect( sqrt(L)^(m-1) for m in eachindex(codiagonals) ))
    diagonals_rescaled   = Float64.(diagonals.*collect( sqrt(L)^(m-2) for m in eachindex(diagonals) ))[1:number_of_shells_needed]
    codiagonals_rescaled = Float64.(codiagonals.*collect( sqrt(L)^(m-1) for m in eachindex(codiagonals) ))[1:number_of_shells_needed-1]

    # enforce asymptotic behavior
    if enforce_asymptotic_behavior
        asymptotize_codiagonals!(codiagonals_rescaled)
    end

    if verbose
        println( "DIAGONALS" )
        @printf "lanczos | rescaled\n"
        for j in 1:number_of_shells
            @printf "%7.3e | %8.3e \n" diagonals[j] diagonals_rescaled[j]
        end
        println()
        println( "CODIAGONALS" )
        @printf "lanczos |  rescaled\n"
        for j in 1:(number_of_shells-1)
            @printf "%7.3e | %8.3e \n"  codiagonals[j] codiagonals_rescaled[j]
        end
        println()
    end

    return ( diagonals_rescaled , codiagonals_rescaled )
end

# ========= #
# Utilities #
# ========= #
function integrate( f::Function , 
                    interval::NTuple{2,R} ;
                    number_of_subintervals::Int64=Int64(1e2) ,
                    method::String="middle" )::R where {R<:Real}

    @assert (method in integrate_methods) "Integration method not implemented."

    # finite value step
    step::BigFloat = (interval[2]-interval[1])/number_of_subintervals

    # interval starting point
    start::BigFloat = interval[1]

    # integrate using chosen method
    if method=="left"
        return integrate_left(f,number_of_subintervals,start,step)
    elseif method=="middle"
        return integrate_middle(f,number_of_subintervals,start,step)
    else 
        error( "Integration method not implemented." )
    end
end
function integrate_left( f::Function ,
                         number_of_subintervals::Int64 ,
                         start::BigFloat ,
                         step::BigFloat )

    integral::BigFloat = 0.0 

    x::BigFloat = start

    for _ in 1:number_of_subintervals
        integral += f(x)*step
        x += step
    end

    return integral
end
function integrate_middle( f::Function ,
                           number_of_subintervals::Int64 ,
                           start::BigFloat ,
                           step::BigFloat )

    integral::BigFloat = 0.0 

    x::BigFloat = start+step/2.0

    for _ in 1:number_of_subintervals
        integral += f(x)*step
        x += step
    end

    return integral
end

function jiter( maximum_j::Int64 )
    return Iterators.filter( !iszero , -maximum_j:maximum_j )
end

# ============= #
# Hybridization #
# ============= #
#
# We DEFINE the hybridization to be
#
#   Γ(ϵ) = π ρ(ϵ) V² ≡ 
#
# with 
#
#   ∫ρ(ϵ)dϵ = 1.
#
# This, in turn, defines ρ(ϵ) as
# an effective DOS for the system 
# and 
#   t≡V 
# as the hopping parameter
# that enters the Hamiltonian.


# The hybridization function Γ(ϵ) is decomposed
# into its constituents, t and ρ(ϵ), by the 
# following relations:
#   V = √(∫Γ(ϵ)dϵ/π),
#   ρ(ϵ) = (Γ(ϵ)/π) / V²
function decompose_hybridization( G::Function )::Tuple{Float64,Function}

    # hopping parameter t=V
    integration_interval = (-1.0,1.0)
    V = sqrt(integrate(G,integration_interval))/pi

    # effective DOS ρ(ϵ)
    effective_dos(x) = G(x) / ( V^2 * pi )

    return ( V , effective_dos )
end


# ============== #
# Discretization #
# ============== #
function discretization_energy( L::BigFloat , j::Int64 , z::BigFloat )::BigFloat

    absolute_energy::BigFloat = abs(j)==1 ? 1.0 : L^(1-abs(j)+z)

    return j<0 ? -absolute_energy : absolute_energy
end

# energy interval 
function energy_interval( L::BigFloat , j::Int64 , z::BigFloat )::Tuple{BigFloat,BigFloat}

    # positive interval
    j>0 && return ( discretization_energy(L,j+1,z) , discretization_energy(L,j,z) )

    # negative interval
    return ( discretization_energy(L,j,z) , discretization_energy(L,j-1,z) )

end

# representative energy 
function representative_energy( dos::Function , 
                                L::BigFloat ,
                                j::Int64 , 
                                z::BigFloat ,
                                discretization::String )::BigFloat

    interval = energy_interval(L,j,z)

    if discretization=="yoshida1990"

        numerator   = integrate( x->dos(x)*x , interval )
        denominator = integrate( dos , interval )
        if iszero(numerator) # use flat DOS
            numerator   = integrate( x->1.0*x , interval )
            denominator = integrate( x->1.0 , interval )
        end

        return numerator/denominator 

    elseif discretization=="campo2005"

        numerator = integrate( dos , interval )
        denominator = integrate( x->dos(x)/x , interval ) 
        if iszero(numerator) # use flat DOS
            numerator = integrate( x->1.0 , interval )
            denominator = integrate( x->1.0/x , interval ) 
        end

        return numerator/denominator
    end
end

# ============================================== #
# Preparation of Hamiltonian and innermost shell #
# ============================================== #
function construct_nontransformed_hamiltonian(
            dos::Function ,
            L::BigFloat ,
            maximum_j::Int64 ,
            z::BigFloat ,
            discretization::String )::Matrix{BigFloat}

    representative_energies = collect( representative_energy(dos,L,j,z,discretization) for j in jiter(maximum_j) )

    return diagm(representative_energies)

end

function innermost_shell_vector( dos::Function ,
                                 L::BigFloat ,
                                 maximum_j::Int64 ,
                                 z::BigFloat )::Vector{BigFloat}

    f0::Vector{BigFloat} = [
        sqrt(integrate( dos , energy_interval(L,j,z) )) 
        for j in jiter(maximum_j)
    ]

    @views LinearAlgebra.normalize!(f0)

    return f0
end

# ================== #
# Tridiagonalization #
# ================== #
function lanczos( H::Matrix{BigFloat} , 
                  v_1::Vector{BigFloat} ,
                  L::BigFloat ; 
                  method="chen1995" ,
                  enforce_particle_hole_symmetry::Bool=false )::Tuple{Vector{BigFloat},Vector{BigFloat}}

    if method=="matrix"

        return lanczos_matrix(H,v_1)

    elseif method=="chen1995"

        return lanczos_chen1995(H,v_1,L,enforce_particle_hole_symmetry)

    elseif method=="gonzalez-buxton1998"

        return lanczos_gonzalezbuxton1998(H,v_1,L;enforce_particle_hole_symmetry=enforce_particle_hole_symmetry)

    end
end

# Yoshida, Whitaker, Oliveira (1990)
function lanczos_matrix( H::Matrix{BigFloat} , v_1::Vector{BigFloat} ) 

    M = size(H,1)

    # unitary diagonalization matrix
    U::Matrix{BigFloat} = similar(H)
    U[1,:] .= v_1

    # iterate
    for n=1:(M-1)

        # produce non-normalized vector 
        U[n+1,:] .= H*U[n,:] - sum( dot(U[m,:],H,U[n,:])*U[m,:] for m in 1:n )

        # normalize vector 
        @views LinearAlgebra.normalize!( U[n+1,:] )

    end

    # obtain tridiagonal hamiltonian
    H_tridiag = similar(H)
    tmp = similar(H)
    mul!( tmp , H , transpose(U) )
    mul!( H_tridiag , U , tmp )

    # extract diagonals and codiagonals
    diagonals,codiagonals = obtain_diagonals_codiagonals( H_tridiag )

    return diagonals,codiagonals

end

# Chen, Jayaprakash (1995)
function lanczos_chen1995( H::Matrix{BigFloat} , 
                           v_1::Vector{BigFloat} ,
                           L::BigFloat ,
                           enforce_particle_hole_symmetry::Bool ;
                           verbose::Bool=false )::Tuple{Vector{BigFloat},Vector{BigFloat}}

    # maximum shell
    N = Int64(size(H,1))
    representative_energies = collect( H[i,i] for i in axes(H,1) )

    U::Matrix{BigFloat} = similar(H)
    U[1,:] .= v_1

    # diagonal and codiagonal entries
    diagonals = zeros(BigFloat,N)
    codiagonals = zeros(BigFloat,N-1)
    # intitialization
    diagonals[1] = enforce_particle_hole_symmetry ? 0.0 : sum( representative_energies[:].*(U[1,:].^2) )
    codiagonals[1] = sqrt(dot( (representative_energies.-diagonals[1]).^2 , U[1,:].^2 ))

    # step 1
    U[2,:] = (representative_energies[:].-diagonals[1]).*U[1,:]/codiagonals[1]
    @views LinearAlgebra.normalize!(U[2,:])
    codiagonals[1] = sum( U[1,:].*U[2,:].*representative_energies[:] )
    diagonals[2] = enforce_particle_hole_symmetry ? 0.0 : sum( representative_energies[:].*(U[2,:].^2) )
    codiagonals[2] = sqrt( sum( (representative_energies[:].^2).*(U[2,:].^2) ) - codiagonals[1]^2 - diagonals[2]^2 )

    # iterate
    for n=2:(N-1)

        # safety thresholds in j
        upper_threshold_j = iseven(n) ? n÷2 : (n÷2+1)
        lower_threshold_j = -upper_threshold_j
        # thresholds as indices
        lower_threshold_idx = Int64(lower_threshold_j+1+N÷2)
        upper_threshold_idx = Int64(upper_threshold_j+N÷2)
        # safety mask
        mask_safe = map( x->(x<=lower_threshold_idx || x>=upper_threshold_idx) , 1:N )
        mask_unsafe = @. !mask_safe

        # produce non-normalized vector f_{n+1} in the safe region
        U[n+1,mask_safe] .= (( representative_energies[mask_safe] .- diagonals[n] ).*U[n,mask_safe] - codiagonals[n-1].*U[n-1,mask_safe] )./codiagonals[n]

        # for the unsafe region, solve orthogonality problem
        unsafe_j_size = sum(mask_unsafe)
        @views A = U[1:unsafe_j_size,mask_unsafe]
        b = - reduce( + , U[1:unsafe_j_size,mask_safe]*U[n+1,mask_safe] , dims=2 )
        if !isapprox(det(A),0.0) 
            U[n+1,mask_unsafe] .= A\b
        else
            println("WARNING: Singular matrix in calculation of U from orthogonality. Applying brute force instead.")
            U[n+1,mask_unsafe] .= (( representative_energies[mask_unsafe] .- diagonals[n] ).*U[n,mask_unsafe] - codiagonals[n-1].*U[n-1,mask_unsafe] )./codiagonals[n]
        end

        # check norm
        if iszero(norm(U[n+1,:]))
            println( "Vector norm = 0. Band decouples from impurity at iteration $(n+1)." )
            codiagonals[n,:] .= 0.0
            diagonals[n+1,:] .= 0.0
            return diagonals, codiagonals
        else
            # normalize vector 
            @views LinearAlgebra.normalize!(U[n+1,:])
        end

        # correct previous codiagonals
        codiagonals[n] = sum( U[n,:].*U[n+1,:].*representative_energies[:] )
        if iszero(codiagonals[n])
            codiagonals[n:end] = 0
            return diagonals,codiagonals
        end
        
        # new diagonal and codiagonals
        diagonals[n+1] = enforce_particle_hole_symmetry ? 0.0 : sum( representative_energies[:].*(U[n+1,:].^2) )
        if n!==(N-1)
            codiagonal_square = dot(representative_energies.^2,U[n+1,:].^2) - codiagonals[n]^2 - diagonals[n+1]^2
            if (!(codiagonal_square>=0) || isnan(codiagonal_square))
                println( "Numerical error in discretization at step $(n): square root of $codiagonal_square" )
                if enforce_particle_hole_symmetry
                    print( "Looking for asymptotic behavior instead... " )
                    for nn in 3:n 
                        if (isapprox(codiagonals[nn],codiagonals[nn-2]/L;rtol=1e-2) && isapprox(codiagonals[nn-2],codiagonals[nn-4]/L;rtol=1e-2))
                            println( "found! Computing codiagonals from asymptotic behavior" )

                            for m in eachindex(codiagonals[(n+1):end])
                                if iseven(m) 
                                    codiagonals[n+m] = codiagonals[n]/(sqrt(L)^(m))
                                else
                                    codiagonals[n+m] = codiagonals[n-1]/sqrt(L)^(m+1)
                                end
                            end

                            return diagonals,codiagonals
                        end
                    end
                    println( "not found. Quitting" )
                end
                error( "Discretization crashed at step $n" )
            end
            codiagonals[n+1] = sqrt(codiagonal_square)
        end

    end

    if verbose 
        println( "DIAGONALS" )
        for i in eachindex(diagonals)
            @printf "%5i %.5e %.5e" i diagonals[i] diagonals[i]*3.0^(i-1)
            println()
        end
        println( "CODIAGONALS" )
        for i in eachindex(codiagonals)
            @printf "%5i %.5e %.5e" i codiagonals[i] codiagonals[i]*sqrt(3.0)^(i-1)
            println()
        end
    end

    return diagonals,codiagonals

end

function lanczos_gonzalezbuxton1998( H::Matrix{BigFloat} , 
                                     v_1::Vector{BigFloat} ,
                                     L::BigFloat ;
                                     verbose::Bool=false ,
                                     enforce_particle_hole_symmetry::Bool=false )::Tuple{Vector{BigFloat},Vector{BigFloat}}

    # maximum shell
    N = Int64(size(H,1))
    representative_energies = collect( H[i,i] for i in axes(H,1) )

    U::Matrix{BigFloat} = similar(H)
    U[1,:] .= v_1

    # diagonal and codiagonal entries
    diagonals = zeros(BigFloat,N)
    codiagonals = zeros(BigFloat,N-1)
    # intitialization
    diagonals[1] = enforce_particle_hole_symmetry ? 0.0 : sum( representative_energies[:].*(U[1,:].^2) )
    codiagonals[1] = sqrt(dot( (representative_energies.-diagonals[1]).^2 , U[1,:].^2 ))

    # step 1
    U[2,:] = (representative_energies[:].-diagonals[1]).*U[1,:]/codiagonals[1]
    @views LinearAlgebra.normalize!(U[2,:])
    codiagonals[1] = sum( U[1,:].*U[2,:].*representative_energies[:] )
    diagonals[2] = enforce_particle_hole_symmetry ? 0.0 : sum( representative_energies[:].*(U[2,:].^2) )
    codiagonals[2] = sqrt( sum( (representative_energies[:].^2).*(U[2,:].^2) ) - codiagonals[1]^2 - diagonals[2]^2 )

    # iterate
    for n=2:(N-1)

        # produce non-normalized vector f_{n+1} in the safe region
        U[n+1,:] .= (( representative_energies[:] .- diagonals[n] ).*U[n,:] - codiagonals[n-1].*U[n-1,:] )./codiagonals[n]

        # check norm
        if iszero(norm(U[n+1,:]))
            println( "Vector norm = 0. Band decouples from impurity at iteration $(n+1)." )
            codiagonals[n,:] .= 0.0
            diagonals[n+1,:] .= 0.0
            return diagonals, codiagonals
        else
            # normalize vector 
            @views LinearAlgebra.normalize!(U[n+1,:])
        end

        # correct previous codiagonals
        codiagonals[n] = sum( U[n,:].*U[n+1,:].*representative_energies[:] )
        if iszero(codiagonals[n])
            codiagonals[n:end] = 0
            return diagonals,codiagonals
        end
        
        # new diagonal and codiagonals
        diagonals[n+1] = enforce_particle_hole_symmetry ? 0.0 : sum( representative_energies[:].*(U[n+1,:].^2) )
        if n!==(N-1)
            codiagonal_square = dot(representative_energies.^2,U[n+1,:].^2) - codiagonals[n]^2 - diagonals[n+1]^2
            if (!(codiagonal_square>=0) || isnan(codiagonal_square))
                println( "Numerical error in discretization at step $(n): square root of $codiagonal_square" )
                if enforce_particle_hole_symmetry
                    print( "Looking for asymptotic behavior instead... " )
                    for nn in 3:n 
                        if (isapprox(codiagonals[nn],codiagonals[nn-2]/L;rtol=1e-2) && isapprox(codiagonals[nn-2],codiagonals[nn-4]/L;rtol=1e-2))
                            println( "found! Computing codiagonals from asymptotic behavior" )

                            for m in eachindex(codiagonals[(n+1):end])
                                if iseven(m) 
                                    codiagonals[n+m] = codiagonals[n]/(sqrt(L)^(m))
                                else
                                    codiagonals[n+m] = codiagonals[n-1]/sqrt(L)^(m+1)
                                end
                            end

                            return diagonals,codiagonals
                        end
                    end
                    println( "not found. Quitting" )
                end
                error( "Discretization crashed at step $n" )
            end
            codiagonals[n+1] = sqrt(codiagonal_square)
        end

    end

    if verbose 
        println( "DIAGONALS" )
        for i in eachindex(diagonals)
            @printf "%5i %.5e %.5e" i diagonals[i] diagonals[i]*3.0^(i-1)
            println()
        end
        println( "CODIAGONALS" )
        for i in eachindex(codiagonals)
            @printf "%5i %.5e %.5e" i codiagonals[i] codiagonals[i]*sqrt(3.0)^(i-1)
            println()
        end
    end

    return diagonals,codiagonals
end


function obtain_diagonals_codiagonals( H::Matrix )::Tuple{Vector{BigFloat},Vector{BigFloat}}

    N = size(H,1)

    diagonals = collect( H[n,n] for n in axes(H,1) )
    codiagonals = collect( H[n,n+1] for n in 1:(N-1) )

    return (diagonals,codiagonals)

end

# ======================= #
# Asymptotic codiagonals  #
# ======================= #
function asymptotize_codiagonals!( codiagonals::Vector{<:Real} ;
                                   rtol::Float64=1e-4 )

    # if insufficient codiagonals, do nothing
    length(codiagonals)<10 && return

    for n in 10:length(codiagonals)

        # double check standard asymptotic behavior
        if (isapprox(codiagonals[n],codiagonals[n-1];rtol=rtol) && isapprox(codiagonals[n-1],codiagonals[n-2];rtol=rtol))
            codiagonals[n:end] .= codiagonals[n]

        # double check asymptotic behavior with even-odd alternation
        elseif (isapprox(codiagonals[n],codiagonals[n-2]) && isapprox(codiagonals[n-2],codiagonals[n-4]) &&
                isapprox(codiagonals[n-1],codiagonals[n-3]) && isapprox(codiagonals[n-3],codiagonals[n-5]))

            codiagonals[n:2:end] .= codiagonals[n]
            codiagonals[(n-1):2:end] .= codiagonals[n-1]

        end
    end
end


# ======================== #
# Check by rediagonalizing #
# ======================== #
function check_tridiagonalization( diagonals::Vector{BigFloat} ,
                                   codiagonals::Vector{BigFloat} ,
                                   H_diagonal::Matrix{BigFloat} ;
                                   enforce_asymptotic_behavior::Bool=enforce_asymptotic_behavior_default ,
                                   L::BigFloat=1.0 ) #L is needed if enforce_particle_hole_symmetry

    codiags = copy(codiagonals)
    if enforce_asymptotic_behavior
        codiags_rescaled = [c*L^((n-1)/2.0) for (n,c) in enumerate(codiagonals)]
        asymptotize_codiagonals!(codiags_rescaled)
        codiags = [cr*L^(-(n-1)/2.0) for (n,cr) in enumerate(codiags_rescaled)]
    end

    H_tridiagonal = similar(H_diagonal)
    H_tridiagonal .= 0.0
    for n in axes(H_tridiagonal,1)
        H_tridiagonal[n,n] = diagonals[n]
        if n!==size(H_tridiagonal,1) 
            H_tridiagonal[n,n+1] = codiags[n]
            H_tridiagonal[n+1,n] = codiags[n]
        end
    end
    F = eigen(Float64.(H_tridiagonal))
    E_rediagonal = sort(F.values)
    E_representative = sort([H_diagonal[n,n] for n in axes(H_diagonal,1)])
    @printf "representative energies | tridiagonal-rediagonal | relative difference\n"
    for (e_representative,e_rediagonal) in zip(E_representative,E_rediagonal)
        @printf "%-23.4e | %-22.4e | %-.4e\n" e_representative e_rediagonal abs(e_representative-e_rediagonal)/e_representative
    end
    println()

end
