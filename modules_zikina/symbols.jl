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

struct HilbertSpace 
    states::Array{Tuple}
end
HilbertSpace( ss... ) = HilbertSpace([ ss... ])


*( h1::HilbertSpace , h2::HilbertSpace ) = 
    HilbertSpace( h1.states... , h2.states... )
length( h::HilbertSpace ) = length( h.states )


function println( hilbert::HilbertSpace )
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

struct SymbolFockState
    # occ = ( gamma , I , mu , sigma )
    occ::Vector{Bool}
    hilbert::HilbertSpace
end
# constructor for singly occupied, empty or zero state
function SymbolFockState( state::Tuple , h::HilbertSpace )
    # zero state
    if state==(0,)
        boolv = Array{Bool}([])
    # empty state 
    elseif state==()
        boolv = zeros( Bool , length(h.states) )
    # singly occupied state
    else
        boolv = map( x -> x==state , h.states )
    end
    return SymbolFockState( boolv , h )
end
# constructor for a list of occupied states 
function SymbolFockState( states::Vector{T} , h::HilbertSpace ) where T<:Tuple
    boolv = map( x -> (x in states) , h.states )
    return SymbolFockState( boolv , h )
end

function (==)( s1::SymbolFockState , s2::SymbolFockState )
    return s1.occ==s2.occ
end

function string( sfs::SymbolFockState )
    # zero
    sfs.occ==[] && return "0"
    # empty
    all( sfs.occ.==0 ) && return "|)"
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

struct SymbolCreationOperator <: SymbolOperator 
    state::Tuple
end

function *( o::SymbolCreationOperator , s::SymbolFockState )
    # if state=0, return 0
    s.occ==[] && return (0,SymbolFockState((0,),s.hilbert))
    # index of state
    idx = findfirst( map(x->x==o.state , s.hilbert.states ) )
    # if state occupied, return 0
    s.occ[idx]==1 && return (0,SymbolFockState((0,),s.hilbert)) 
    # sign due to permutations
    p = sum( s.occ[1:(idx-1)] )
    sign = (-1)^p
    # occupation to add 
    newpart = map( x -> x==idx , 1:length(s.occ) )
    return ( sign , SymbolFockState(s.occ+newpart,s.hilbert) )
end


# ANNIHILATION OPERATOR #

struct SymbolAnnihilationOperator <: SymbolOperator 
    state::Tuple
end


function *( o::SymbolAnnihilationOperator , s::SymbolFockState )
    # if state=0, return 0
    s.occ==[] && return (0,SymbolFockState((0,),s.hilbert))
    # index of state 
    idx = findfirst( map(x->x==o.state , s.hilbert.states ) )
    # if state empty, return 0
    s.occ[idx]==0 && return (0,SymbolFockState((0,),s.hilbert)) 
    # sign due to permutations 
    p = sum( s.occ[1:(idx-1)] )
    sign = (-1)^p
    # occupation to subtract
    newpart = map( x -> x==idx , 1:length(s.occ) )
    return ( sign , SymbolFockState(s.occ-newpart,s.hilbert) )
end


# IDENTITY OPERATOR #

struct SymbolIdentityOperator <: SymbolOperator end
function *( o::SymbolIdentityOperator , s::SymbolFockState )
    return (1,s)
end
I = SymbolIdentityOperator()

# printing
string( sco::SymbolCreationOperator ) = "[c| $(sco.state[1]) $(sco.state[2])"* 
                                           " $(sco.state[3]) $(sco.state[4]) ]"
function println( sco::SymbolCreationOperator )
    println( string(sco) )
end
string( sco::SymbolAnnihilationOperator ) = "[a| $(sco.state[1]) $(sco.state[2])"* 
                                               " $(sco.state[3]) $(sco.state[4]) ]"
function println( sco::SymbolAnnihilationOperator )
    println( string(sco) )
end
string( sco::SymbolIdentityOperator ) = "I"
function println( sco::SymbolIdentityOperator )
    println( string(sco) )
end



# ############# #
# MISCELLANEOUS #
# ............. #

braket( bra::SymbolFockState , ket::SymbolFockState ) = bra==ket ? 1 : 0
@inbounds function sandwich( bra::SymbolFockState , op::SymbolOperator , ket::SymbolFockState )
    opket = op*ket
    return braket( bra, opket[2] )==0 ? 0 : opket[1]
end
