include( "symbols.jl" )
using Combinatorics
using LinearAlgebra
using SharedArrays
using Distributed

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
struct CanonicalBasis <: Basis
    states::Array{SymbolFockState,1}
end
function CanonicalBasis( h::HilbertSpace )
    states =  map( x -> SymbolFockState(x,h) , combinations(h.states) )
    prepend!( states , [SymbolFockState((),h)] )
    return CanonicalBasis( states )
end
function CanonicalBasis( h::HilbertSpace , n::Integer )
    # Canonical basis with N=n
    if n==0
        states = [SymbolFockState( () , h )]
    else
        states =  map( x -> SymbolFockState(x,h) , combinations(h.states,n) )
    end
    return CanonicalBasis( states )
end

function *( b1::CanonicalBasis , b2::CanonicalBasis )
    hilbert = b1.states[1].hilbert * b2.states[1].hilbert
    states = []
    for scf1 in b1.states, scf2 in b2.states
        occ = [ scf1.occ... , scf2.occ... ]
        scf = SymbolFockState( occ , hilbert )
        push!( states , scf )
    end
    return CanonicalBasis( states )
end

+( b1::CanonicalBasis , b2::CanonicalBasis ) = 
    CanonicalBasis([ b1.states... , b2.states... ])

(==)( b1::CanonicalBasis , b2::CanonicalBasis ) = 
    all( b1.states.==b2.states )

function println( b::CanonicalBasis ) 
    for (i,s) in enumerate(b.states)
        println( "$i: " , string(s) )
    end 
end

# GENERIC BASIS: Basis in which states and operators 
#   may be represented after a basis transformation.
struct GenericBasis <: Basis 
    states::Array
end
@inbounds function GenericBasis( cb::CanonicalBasis , u::Matrix )
    states = [State(cb) for i=1:length(cb.states)]
    u = u'
    for i=1:size(u,1)
        for j=1:size(u,2)
            states[i] += u[i,j]*State(cb.states[j],cb)
        end
    end
    return GenericBasis( states )
end
GenericBasis( cb::CanonicalBasis ) = GenericBasis( cb , 
                                     Matrix(LinearAlgebra.I(length(cb.states))) )

function println( b::GenericBasis ) 
    for (i,s) in enumerate(b.states)
        print( "$i: " )
        println( s )
    end 
end


# #####################################
# STATE 
#
# - Second quantization object.
# - Represented by a vector in a basis. 
# .....................................

struct State 
    vector::Array{ComplexF64}
    basis::CanonicalBasis
end
function State( s::SymbolFockState , cb::CanonicalBasis )
    i = findfirst( x -> x==s , cb.states )
    return State( Array{ComplexF64}([ j==i ? 1 : 0 for j=1:length(cb.states) ]) , cb )
end
function State( st::Tuple , hilbert::HilbertSpace , basis::CanonicalBasis )
    return State( SymbolFockState(st,hilbert) , basis )
end
State( b::Basis ) = State( zeros(length(b.states)) , b )


+( s1::State , s2::State ) = State( s1.vector + s2.vector , s1.basis )
-( s1::State , s2::State ) = s1 + (-1)*s2
*( a::Number , s::State ) = State( a*s.vector , s.basis )
*( s::State , a::Number ) = a*s
*( bra::State , ket::State ) = bra.vector' * ket.vector
function x( s1::State , s2::State )
    v1 = s1.vector 
    v2 = s2.vector 
    v  = []
    for e1 in v1, e2 in v2 
        append!( v , e1*e2 )
    end
    return State( vec(v) , s1.basis*s2.basis )
end

function extend_basis( s::State , newbasis::CanonicalBasis )
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


# ######## #
# OPERATOR #
# ........ #

mutable struct Operator 
    matrix::Array{ComplexF64}
    basis::Basis
end
@inbounds function Operator( o::SymbolOperator , basis::CanonicalBasis ; herm=false )
    mat = Matrix{ComplexF64}(undef,length(basis.states),length(basis.states))
    if herm 
        for i=1:length(basis.states), j=i:length(basis.states)
            mat[i,j] = sandwich( basis.states[i] , o , basis.states[j] )
            i!==j && (mat[j,i] = conj( mat[i,j] ))
        end
    else
        mat = [ sandwich(basis.states[i],o,basis.states[j])
                for i=1:length(basis.states), j=1:length(basis.states) ]
    end
    return Operator( mat , basis )
end
Operator( a::Number , basis::CanonicalBasis ) = a * Operator( I , basis )
@inbounds function Operator( o::Operator , basis::GenericBasis , herm=false )
    mat = Matrix{ComplexF64}(length(basis.states),length(basis.states))
    if herm 
        for i=1:length(basis.states)
            for j=i:length(basis.states)
                mat[i,j] = basis.states[i] * o * basis.states[j]
                i!==j && (mat[j,i] = conj( mat[i,j] ))
            end
        end
    else
        mat = [ ComplexF64(basis.states[i]*o*basis.states[j]) 
                for i=1:length(basis.states), j=1:length(basis.states)]
    end
    return Operator( mat , basis )
end

+( o1::Operator , o2::Operator ) = Operator( o1.matrix + o2.matrix , o1.basis )
-( o1::Operator , o2::Operator ) = o1 + (-1)*o2
*( o1::Operator , o2::Operator ) = Operator( o1.matrix*o2.matrix , o1.basis )
*( o::Operator , s::State ) = State( o.matrix*s.vector , o.basis )
*( s::State , o::Operator ) = State( o.matrix'*s.vector , o.basis )
*( a::Number , o::Operator ) = Operator( a*o.matrix , o.basis )
*( o::Operator , a::Number ) = a*o
+( o::Operator , a::Number ) = o + Operator( a , o.basis )
+( a::Number , o::Operator ) = o + a
-( o::Operator , a::Number ) = o + Operator( -a , o.basis )
-( a::Number , o::Operator ) = Operator( a , o.basis ) - o
^( o::Operator , p::Integer ) = p==2 ? o*o : o*(o^(p-1))
adjoint( o::Operator ) = Operator( o.matrix' , o.basis )
length( o::Operator ) = length(o.matrix)

using LinearAlgebra: eigen
function diagonalize( o::Operator )
    o = 0.5*( o + o' )
    F = eigen( o.matrix )
    (e, u) = ( real(F.values) , F.vectors )
    return ( e , u )
end

function hermitize!( o::Operator )
    o.matrix = 0.5*( o.matrix + o.matrix' )
end


string( o::Operator ) = string( o.matrix )
function println( o::Operator )
    println( string(o) )
end

function block( o::Operator , gb::GenericBasis )
    n = length(gb.states)
    mat = Array{ComplexF64,2}( undef , n , n )
    for i=1:n, j=1:n 
        mat[i,j] = gb.states[i] * o * gb.states[j]
    end
    return Operator( mat , gb )
end

function print_diagonals( o::Operator )
    for (i,sfs) in enumerate(o.basis.states)
        println( sfs )
        println( o.matrix[i,i] )
        println()
    end
end
    


