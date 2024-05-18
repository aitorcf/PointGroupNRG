import Base.:(==)
import Base.:*
import Base.string
import Base.println
import Base.findfirst
import Base.length

# -------------------------------------------------------- 
# This file contains structures and functions used 
# in computations involving symbolic states and operators,
# i.e. those that are not yet numericaly represented by 
# vectors and matrices, respectively.
# --------------------------------------------------------


# ##########################################
# HILBERT SPACE
#
# - Provides (i) the one-particle states and
#            (ii) their ordering
# ...........................................

struct HilbertSpace{T<:Tuple}
    states::Vector{T}
end
HilbertSpace( ss... ) = HilbertSpace([ ss... ])
HilbertSpace() = HilbertSpace(Tuple{}[])


*( h1::HS , h2::HS ) where {HS<:HilbertSpace} = HilbertSpace([ h1.states... , h2.states... ])
length( h::HS ) where {HS<:HilbertSpace} = length( h.states )


function println( hilbert::HS ) where {HS<:HilbertSpace}
    for s in hilbert.states 
        println(s)
    end 
end



# #####################################
# SYMBOL FOCK STATE 
#
# - A single canonical basis state.
# - Operators symbolically act on this.
# .....................................

struct SymbolFockState{N<:Number,HS<:HilbertSpace}
    # occ = ( gamma , I , mu , sigma )
    occ::Vector{N}
    hilbert::HS
end
## general constructor (overriden by specific dispatch methods below)
#function SymbolFockState( state::T , h::HilbertSpace{T} ) where {T<:Tuple}
#    # zero state
#    if state==(0,)
#        boolv = Array{Bool}([])
#    # empty state 
#    elseif state==()
#        boolv = zeros( Bool , length(h.states) )
#    # singly occupied state
#    else
#       boolv = map( x -> x==state , h.states )
#    end
#    return SymbolFockState( boolv , h )
#end
# constructor for singly occupied state
function SymbolFockState( state::T , h::HilbertSpace{T} ) where {T<:Tuple}
    boolv = map( x -> x==state , h.states )
    return SymbolFockState( boolv , h )
end
# constructor for zero state
function SymbolFockState( state::Tuple{Int64} , h::HilbertSpace{T} ) where {T<:Tuple}
    boolv = Bool[]
    return SymbolFockState( boolv , h )
end
# constructor for empty state
function SymbolFockState( state::Tuple{} , h::HilbertSpace{T} ) where {T<:Tuple}
    boolv = zeros( Bool , length(h.states) )
    return SymbolFockState( boolv , h )
end
# constructor for a list of occupied states 
function SymbolFockState( states::Vector{T} , h::HilbertSpace{T} ) where T<:Tuple
    ss = Set(states)
    return SymbolFockState( map( x -> (x in states) , h.states ) , h )
end

@inbounds function (==)( s1::SFS , s2::SFS ) where {SFS<:SymbolFockState}
    return s1.occ==s2.occ
end

function string( sfs::SFS ) where {SFS<:SymbolFockState}
    # zero
    length(sfs.occ)==0 && return "0"
    # empty
    sum( sfs.occ )==0 && return "|)"
    s = []
    for (i,stuple) in enumerate(sfs.hilbert.states)
        sfs.occ[i]==0 && continue
        push!( s , "$(stuple[1]) $(stuple[2]) $(stuple[3]) $(stuple[4])" )
    end 
    return "| " * join( s , " , " ) * " )"
end
function println( sfs::SymbolFockState )
    println( string( sfs ) )
end



# ###########################################
# SYMBOL OPERATORS
#
# - Single creation/annihilation operator.
# - They act symbolically on symbolic states.
# ...........................................

abstract type SymbolOperator end


# CREATION OPERATOR #

struct SymbolCreationOperator{T<:Tuple} <: SymbolOperator 
    state::T
end

function *( o::SymbolCreationOperator{T}, s::SymbolFockState{N,HilbertSpace{T}} ) where {N<:Number,T<:Tuple}

    # if state=0, return 0
    length(s.occ)==0 && return (0,SymbolFockState((0,),s.hilbert))

    # index of state
    idx = findfirst( x->x==o.state , s.hilbert.states ) 

    # if state occupied, return 0
    s.occ[idx]==1 && return ( 0, SymbolFockState((0,),s.hilbert) ) 

    # sign due to permutations
    p = sum( s.occ[1:(idx-1)] )
    sign = (-1)^p

    # occupation to add 
    #newpart = map( x -> x==idx , 1:length(s.occ) )
    newocc = copy( s.occ ) 
    newocc[idx] = 1

    #return ( sign , SymbolFockState(s.occ+newpart,s.hilbert) )
    return ( sign , SymbolFockState(newocc,s.hilbert) )
end


# ANNIHILATION OPERATOR #

struct SymbolAnnihilationOperator{T<:Tuple} <: SymbolOperator 
    state::T
end

function *( o::SymbolAnnihilationOperator{T} , s::SymbolFockState{N,HilbertSpace{T}} ) where {N<:Number,T<:Tuple}

    # if state=0, return 0
    length(s.occ)==0 && return (0,SymbolFockState((0,),s.hilbert))

    # index of state 
    idx = findfirst( map(x->x==o.state , s.hilbert.states ) )

    # if state empty, return 0
    s.occ[idx]==0 && return (0,SymbolFockState((0,),s.hilbert)) 

    # sign due to permutations 
    p = sum( s.occ[1:(idx-1)] )
    sign = (-1)^p

    # occupation to subtract
    #newpart = map( x -> x==idx , 1:length(s.occ) )
    newocc = copy( s.occ ) 
    newocc[idx] = 0

    return ( sign , SymbolFockState(s.occ-newpart,s.hilbert) )
end


# IDENTITY OPERATOR #

struct SymbolIdentityOperator <: SymbolOperator end
function *( o::SymbolIdentityOperator , s::SFS )::Tuple{ComplexF64,SymbolFockState} where {SFS<:SymbolFockState}
    return (1,s)
end
#I = SymbolIdentityOperator()

# printing
string( sco::SCO ) where {SCO<:SymbolCreationOperator} = 
        "[c| $(sco.state[1]) $(sco.state[2])"* 
        " $(sco.state[3]) $(sco.state[4]) ]"
function println( sco::SCO ) where {SCO<:SymbolCreationOperator}
    println( string(sco) )
end
string( sco::SAO ) where {SAO<:SymbolAnnihilationOperator} = 
        "[a| $(sco.state[1]) $(sco.state[2])"* 
        " $(sco.state[3]) $(sco.state[4]) ]"
function println( sco::SymbolAnnihilationOperator )
    println( string(sco) )
end
string( sco::SymbolIdentityOperator ) = "I"
function println( sco::SIO ) where {SIO<:SymbolIdentityOperator}
    println( string(sco) )
end



# ############# #
# MISCELLANEOUS #
# ............. #

@inbounds function braket( bra::SFS , ket::SFS )::Int64 where {SFS<:SymbolFockState} 
    return bra==ket ? 1 : 0
end
@inbounds function sandwich( bra::SFS , op::SO , ket::SFS )::ComplexF64 where {SFS<:SymbolFockState,SO<:SymbolOperator}
    opket::Tuple{ComplexF64,SFS} = op*ket
    braket_nocoeff::Int64 = braket( bra, opket[2] )
    return braket_nocoeff==zero(braket_nocoeff) ? zero(ComplexF64) : opket[1]
end
