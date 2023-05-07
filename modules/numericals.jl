include( "symbols.jl" )

using Combinatorics
using LinearAlgebra

import Base.:* 
import Base.:+
import Base.:-
import Base.:^
import Base.:(==)
import Base.string
import Base.println 
import Base: adjoint


# ##### #
# BASIS #
# ..... #

abstract type Basis end

# CANONICAL BASIS: Basis on which states and operators 
#   are going to be constructed for the first time.
struct CanonicalBasis{SFS<:SymbolFockState} <: Basis
    states::Vector{SFS}
end
function sortfunction( state::SFS ) where {SFS<:SymbolFockState}

    hilbert = state.hilbert 

    len = sum( state.occ )
    orbitals = []
    for idx in findall(state.occ) 
        push!( orbitals , hilbert.states[idx][3] )
    end
    orbitals = (orbitals...,)

    spindict = Dict( "u"=>1 , "d"=>2 )
    spins = []
    for idx in findall(state.occ) 
        push!( spins , spindict[hilbert.states[idx][4]] )
    end
    spins = (spins...,)

    return (len,orbitals,spins)
end

function CanonicalBasis( h::HS ) where {HS<:HilbertSpace}
    states =  map( x->SymbolFockState(x,h) , combinations(h.states) )
    sort!( states , by=x->sortfunction(x) )
    prepend!( states , [SymbolFockState((),h)] )
    return CanonicalBasis( states )
end
function CanonicalBasis( h::HS , n::Int64 ) where {HS<:HilbertSpace}
    # Canonical basis with N=n
    if n==0
        states = [SymbolFockState( () , h )]
    else
        states =  map( x -> SymbolFockState(x,h) , combinations(h.states,n) )
    end
    return CanonicalBasis( states )
end

function *( b1::CanonicalBasis{SFS} , b2::CanonicalBasis{SFS} )::CanonicalBasis{SFS} where {SFS<:SymbolFockState}

    hilbert = b1.states[1].hilbert * b2.states[1].hilbert
    states::Vector{SFS} = []

    occlength1::Int64 = length(b1.states[1].occ)
    occlength2::Int64 = length(b2.states[1].occ)
    occ::Vector{Bool} = [false for _ in 1:(occlength1+occlength2)]

    for scf1::SFS in b1.states, 
        scf2::SFS in b2.states

        occ[1:occlength1]                           .= scf1.occ
        occ[(occlength1+1):(occlength1+occlength2)] .= scf2.occ
        #occ::Vector{Bool} = vcat( scf1.occ::Vector{Bool} , scf2.occ::Vector{Bool} )
        scf::SFS = SymbolFockState( copy(occ) , hilbert )
        push!( states , scf )

    end
    return CanonicalBasis( states )
end

function +( b1::CanonicalBasis{SFS} , b2::CanonicalBasis{SFS} )::CanonicalBasis{SFS} where {SFS<:SymbolFockState}
    return CanonicalBasis([ b1.states... , b2.states... ])
end

@inbounds function (==)( b1::CanonicalBasis{SFS} , b2::CanonicalBasis{SFS} )::Bool where {SFS<:SymbolFockState}
    return all( b1.states.==b2.states )
end

function println( b::CanonicalBasis{SFS} ) where {SFS<:SymbolFockState}
    for (i,s) in enumerate(b.states)
        println( "$i: " , string(s) )
    end 
end



# #####################################
# STATE 
#
# - Second quantization object.
# - Represented by a vector in a basis. 
# .....................................

struct State{B<:Basis} 
    vector::Vector{ComplexF64}
    basis::B
end
function State( s::SFS , cb::CanonicalBasis{SFS} )::State{CanonicalBasis{SFS}} where {SFS<:SymbolFockState}
    i = findfirst( x -> x==s , cb.states )
    return State( Vector{ComplexF64}([ j==i ? ComplexF64(1.0) : zero(ComplexF64) for j=1:length(cb.states) ]) , cb )
end
function State( st::T , 
                hilbert::HilbertSpace{T} , 
                basis::CanonicalBasis{SymbolFockState{N,HilbertSpace{T}}} )::State{CanonicalBasis{SymbolFockState{N,HilbertSpace{T}}}} where {T<:Tuple,N<:Number}
    return State( SymbolFockState(st,hilbert) , basis )
end
function State( b::B )::State{B} where {B<:Basis} 
    return State( zeros(ComplexF64,length(b.states)) , b )
end


function +( s1::S , s2::S )::S where {S<:State} 
    return State( s1.vector .+ s2.vector , s1.basis )
end
function -( s1::S , s2::S )::S where {S<:State} 
    return s1 + (-1)*s2
end
function *( a::N , s::S )::S where {N<:Number,S<:State} 
    return State( a.*s.vector , s.basis )
end
function *( s::S , a::N )::S where {N<:Number,S<:State} 
    return a*s
end
@inbounds function *( bra::S , ket::S ) where {S<:State} 
    return bra.vector' * ket.vector
end
# tensor product
function x( s1::State{CB} , s2::State{CB} )::State{CB} where {CB<:CanonicalBasis}
    v1 = @view s1.vector[:] 
    v2 = @view s2.vector[:] 
    v = Vector{ComplexF64}(undef,length(v1)*length(v2))
    i::Int64=0
    for e1::ComplexF64 in v1, 
        e2::ComplexF64 in v2 

        i+=1
        v[i] = e1*e2
        #append!( v , e1*e2 )

    end
    return State( vec(v) , s1.basis*s2.basis )
end

function extend_basis( s::State , 
                       newbasis::CanonicalBasis )
    firstindex = findfirst(map(x->x==s.basis.states[1] , newbasis.states ))
    lastindex  = findfirst(map(x->x==s.basis.states[end] , newbasis.states ))
    leading  = firstindex==1 ? [] : [0 for _=1:(firstindex-1)]
    trailing = lastindex==length(newbasis.states) ? [] : 
                [0 for _=(lastindex+1):length(newbasis.states)] 
    v = [ leading... , s.vector... , trailing... ]
    return State( v , newbasis )
end

vacuum( b::CanonicalBasis ) = State( b.states[1] , b )

function normalize!( s::State ) 
    norm = sqrt(sum(s.vector.^2))
    s.vector ./= norm
end


function string( state::State )
    sum([ abs(x) for x in state.vector])<1e-6 && return "0"
    basis = state.basis
    strings = []
    for (i,x) in enumerate(state.vector)
        abs(x)<1e-2 && continue
        push!( strings , "($x) $(string(basis.states[i]))" )
    end
    return join( strings , " + " )
end
function println( state::State )
    println( string(state) )
end

# * GENERIC BASIS: Basis in which states and operators 
# * may be represented after a basis transformation.
struct GenericBasis{S<:State} <: Basis 
    states::Array{S}
end
@inbounds function GenericBasis( cb::CanonicalBasis{SFS} , u::Matrix{ComplexF64} ) where {SFS<:SymbolFockState}
    states::Vector{State} = [State(cb) for _=1:length(cb.states)]
    u = u'
    for i=1:size(u,1),
        j=1:size(u,2)

        states[i] += u[i,j]*State(cb.states[j],cb)

    end
    return GenericBasis( states )
end
function GenericBasis( cb::CanonicalBasis{SFS} ) where {SFS<:SymbolFockState} 
    return GenericBasis( cb , Matrix{ComplexF64}(LinearAlgebra.I(length(cb.states))) )
end
function println( b::GenericBasis{S} ) where {S<:State}
    for (i::Int64,s::S) in enumerate(b.states)
        print( "$i: " )
        println( s )
    end 
end


# ######## #
# OPERATOR #
# ........ #

mutable struct Operator{B<:Basis}
    matrix::Matrix{ComplexF64}
    basis::B
end
@inbounds function Operator( o::SO , basis::CB ; herm=false )::Operator{CB} where {SO<:SymbolOperator,CB<:CanonicalBasis}
    mat = Matrix{ComplexF64}(undef,length(basis.states),length(basis.states))
    M = length(basis.states)
    if herm 
        for i=1:length(basis.states), j=i:length(basis.states)
            mat[i,j] = sandwich( basis.states[i] , o , basis.states[j] )
            i!==j && (mat[j,i] = conj( mat[i,j] ))
        end
    else
        @inbounds for j=1:M, 
                      i=1:M
            mat[i,j] = sandwich(basis.states[i],o,basis.states[j])
        end
    end
    return Operator( mat , basis )
end
function Operator( a::N , basis::CB )::Operator{CB} where {N<:Number,CB<:CanonicalBasis} 
    return a * Operator( I , basis )
end
function Operator( o::Operator{CB} , basis::GB , herm=false )::Operator{GB} where {CB<:CanonicalBasis,GB<:GenericBasis}
    mat = Matrix{ComplexF64}(length(basis.states),length(basis.states))
    if herm 
        @inbounds for i=1:length(basis.states),
                      j=i:length(basis.states)

            mat[i,j] = basis.states[i] * o * basis.states[j]
            i!==j && (mat[j,i] = conj( mat[i,j] ))

        end
    else
        @inbounds for i=1:length(basis.states),
                      j=1:length(basis.states)

            mat[i,j] = basis.states[i] * o * basis.states[j]

        end
    end
    return Operator( mat , basis )
end

@inbounds function +( o1::O , o2::O )::O where {O<:Operator} 
    return Operator( o1.matrix + o2.matrix , o1.basis )
end
@inbounds function -( o1::O , o2::O )::O where {O<:Operator} 
    return o1 + (-1)*o2
end
@inbounds function *( o1::O , o2::O )::O where {O<:Operator} 
    return Operator( o1.matrix*o2.matrix , o1.basis )
end
@inbounds function *( o::Operator{B} , s::State{B} )::State{B} where {B<:Basis} 
    return State( o.matrix*s.vector , o.basis )
end
@inbounds function *( s::State{B} , o::Operator{B} )::State{B} where {B<:Basis} 
    return State( o.matrix'*s.vector , o.basis )
end
@inbounds function *( a::N , o::O )::O where {N<:Number,O<:Operator} 
    return Operator( a*o.matrix , o.basis )
end
function *( o::O , a::N )::O where {N<:Number,O<:Operator} 
    return a*o
end
@inbounds function +( o::O , a::N )::O where {N<:Number,O<:Operator} 
    return o + Operator( a , o.basis )
end
function +( a::N , o::O )::O where {N<:Number,O<:Operator} 
    return o + a
end
function -( o::O , a::N )::O where {N<:Number,O<:Operator} 
    return o + Operator( -a , o.basis )
end
@inbounds function -( a::N , o::O )::O where {N<:Number,O<:Operator} 
    return Operator( a , o.basis ) - o
end
function ^( o::O , p::Int64 )::O where {O<:Operator} 
    return p==2 ? o*o : o*(o^(p-1))
end
function adjoint( o::O )::O where {O<:Operator} 
    return Operator( Matrix{ComplexF64}(o.matrix') , o.basis )
end
function length( o::O )::Int64 where {O<:Operator} 
    return length(o.matrix)
end

using LinearAlgebra: eigen
function diagonalize( o::O )::Tuple{Vector{Float64},Matrix{ComplexF64}} where {O<:Operator}
    o::O = 0.5*( o + o' )
    F = eigen( o.matrix )
    (e, u) = ( real(F.values) , F.vectors )
    return ( e , u )
end

function hermitize!( o::O ) where {O<:Operator}
    o.matrix = 0.5*( o.matrix + o.matrix' )
end


function string( o::O ) where {O<:Operator} 
    return string( o.matrix )
end
function println( o::O ) where {O<:Operator}
    println( string(o) )
end

function block( o::O , gb::GB ) where {O<:Operator,GB<:GenericBasis}
    n = length(gb.states)
    mat = Matrix{ComplexF64}( undef , n , n )
    for i=1:n, 
        j=1:n 

        mat[i,j] = gb.states[i] * o * gb.states[j]

    end
    return Operator( mat , gb )
end

function print_diagonals( o::O ) where {O<:Operator}
    for (i,sfs) in enumerate(o.basis.states)
        println( sfs )
        println( o.matrix[i,i] )
        println()
    end
end
