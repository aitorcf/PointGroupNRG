#%%
using Combinatorics 
using LinearAlgebra
import Base.show
import Base.*
import Base.+
import Base.-
import Base.(==)
import Base.getindex
import Base.size
import Base.display
import Base.length
import Base.sort
import Base.iszero
import LinearAlgebra.:dot

# ##############################################
# BASIS ELEMENT 
#
# - This structure represents a canonical 
#   basis element, which is a simple ket 
#   of the form
#       | 1 2 1 )   [orbital]
#       | ↑ ↑ ↓ )   [spin]
#
# - It only contains a vector of the appropriate
#   type, but it can be string-represented
#   in a visual manner as a ket.
#
# - Structure:
#       BasisElement 
#
#       SimpleBasisElement   <:BasisElement
#       CompoundBasisElement <: BasisElement
#
#       OrbitalBasisElement <: SimpleBasisElement
#       SpinBasisElement    <: SimpleBasisElement
#
# ...............................................

# *************
# BASIS ELEMENT
# .............
abstract type BasisElement end

# CONSTRUCTORS
BasisElement( occ::Vector{Int64} ) = OrbitalBasisElement(occ)
BasisElement( occ::Vector{Char}  ) = SpinBasisElement(occ)
BasisElement( occ::Vector{Nothing} ; target="spin") = target=="spin" ? SpinZeroBasisElement(occ) : OrbitalZeroBasisElement(occ)

# OPERATIONS
# equality
function (==)( be1::T , be2::T ) where {T<:BasisElement}
    return be1.occupations==be2.occupations
end
# dot (scalar) product or braket
dot( be1::T , be2::T ) where {T<:BasisElement} = be1==be2 ? 1 : 0

# PRINTING
function show( io::IO , be::T ) where {T<:BasisElement}
    print( io , string(be) )
end

# ********************
# SIMPLE BASIS ELEMENT
# ....................
abstract type SimpleBasisElement <: BasisElement end

# OPERATIONS
# order and sorting
# |↑ ↓ ) -> sorted; |1 2) -> sorted;
# |↓ ↑ ) -> not sorted,|2 1) -> not sorted.
is_ordered( be::T ) where {T<:SimpleBasisElement} = be==sort(be)
sort( be::T ) where {T<:SimpleBasisElement} = BasisElement( sort(be.occupations) )
# length (particle number)
length( be::T ) where {T<:SimpleBasisElement} = length(be.occupations)
# getindex
getindex( be::T , i::Int64 ) where {T<:SimpleBasisElement} = be.occupations[i]

# DISPLAY
function string( be::T ) where {T<:SimpleBasisElement}
    return "| " * reduce(*,map(x->"$x ",be.occupations))*")"
end

# *********************
# ORBITAL BASIS ELEMENT
# .....................
struct OrbitalBasisElement <: SimpleBasisElement
    occupations::Vector{Int64} 
end

# **************************
# ORBITAL-NULL BASIS ELEMENT
# ..........................
struct OrbitalNullBasisElement <: SimpleBasisElement
    occupations::Vector{Nothing} 
end

function show( io::IO , be::OrbitalNullBasisElement )
    print(io, "| "*"- "^length(be.occupations)*")")
end

# ******************
# SPIN BASIS ELEMENT
# ..................
struct SpinBasisElement <: SimpleBasisElement
    occupations::Vector{Char}
end

# OPERATIONS
# get z-component of spin
function spinz( sbe::SpinBasisElement ) 
    return sum(map(x->x=='↑' ? 0.5 : -0.5,[sbe[i] for i in 1:length(sbe)]))
end

# ********************
# SPIN-0 BASIS ELEMENT
# ....................
struct SpinZeroBasisElement <: SimpleBasisElement
    occupations::Vector{Nothing}
end

sort( be::SpinZeroBasisElement ) = be
function string( be::T ) where {T<:SpinZeroBasisElement}
    return "| " * reduce(*,map(x->"- ",be.occupations))*")"
end

# ******************
# SPIN BASIS ELEMENT
# ..................
struct JBasisElement <: SimpleBasisElement
    occupations::Vector{Float64}
end

# OPERATIONS
# get z-component of total angular momentum
function jz( sbe::JBasisElement ) 
    return sum(sbe.occupations)
end

# *******************************
# ABSTRACT COMPOUND BASIS ELEMENT
# ...............................
abstract type CompoundBasisElementGeneral <: BasisElement end

# order and sorting (see above)
sort( be::CBEG ) where {CBEG<:CompoundBasisElementGeneral} = CompoundBasisElement( sort(be.orbital) , sort(be.spin) )
isorted( be::CBEG ) where {CBEG<:CompoundBasisElementGeneral} = be.orbital==sort(be.orbital) || be.spin==sort(be.spin)

# **********************
# COMPOUND BASIS ELEMENT
# ......................
struct CompoundBasisElement <: CompoundBasisElementGeneral 
    orbital::OrbitalBasisElement
    spin::SpinBasisElement
end

# OPERATIONS
# equality
(==)( be1::CompoundBasisElement , be2::CompoundBasisElement ) =
    ( be1.orbital==be2.orbital && be1.spin==be2.spin )
(==)( be_c::CompoundBasisElement , be_os::Tuple{OrbitalBasisElement,SpinBasisElement} ) =
    ( be_c.orbital==be_os[1] && be_c.spin==be_os[2] )

# DISPLAY
function string( be::CompoundBasisElement )
    opart = "| " * reduce(*,map(x->"$x ",be.orbital.occupations))*")" 
    spart = "| " * reduce(*,map(x->"$x ",be.spin.occupations))*")" 
    return opart * spart
end

# ***********************************
# DOUBLE-GROUP COMPOUND BASIS ELEMENT
# ...................................
struct DoubleGroupCompoundBasisElement <: CompoundBasisElementGeneral 
    orbital::OrbitalBasisElement
    spin::SpinZeroBasisElement
end

# OPERATIONS
# equality
(==)( be1::DoubleGroupCompoundBasisElement , be2::DoubleGroupCompoundBasisElement ) =
    ( be1.orbital==be2.orbital )
(==)( be_c::DoubleGroupCompoundBasisElement , be_os::Tuple{OrbitalBasisElement,SpinZeroBasisElement} ) =
    ( be_c.orbital==be_os[1] )

# DISPLAY
function string( be::DoubleGroupCompoundBasisElement )
    opart = "| " * reduce(*,map(x->"$x ",be.orbital.occupations))*")" 
    spart = "| " * reduce(*,map(x->"- ",be.spin.occupations))*")" 
    return opart * spart
end

# ************************
# J COMPOUND BASIS ELEMENT
# ........................
struct JCompoundBasisElement <: CompoundBasisElementGeneral 
    orbital::OrbitalNullBasisElement
    spin::JBasisElement
end

# OPERATIONS
# equality
(==)( be1::JCompoundBasisElement , be2::JCompoundBasisElement ) =
    ( be1.spin==be2.spin )
(==)( be_c::JCompoundBasisElement , be_os::Tuple{OrbitalNullBasisElement,JBasisElement} ) =
    ( be_c.spin==be_os[2] )

# DISPLAY
function string( be::JCompoundBasisElement )
    opart = "| " * reduce(*,map(x->"- ",be.orbital.occupations))*")" 
    spart = "| " * reduce(*,map(x->"$x ",be.spin.occupations))*")" 
    return opart * spart
end

# #########################################################
# CANONICAL BASIS 
# 
# - It is the basis where states will be represented
#   as vectors.
#
# - It is represented as a vector of BasisElement 
#   components.
#
# - Structure:
#
#       CanonicalBasis 
#
#       SimpleCanonicalBasis   [abstract] <: CanonicalBasis
#       CompoundCanonicalBasis [concrete] <: CanonicalBasis
#
#       OrbitalCanonicalBasis  <: SimpleCanonicalBasis
#
# #########################################################

# ***************
# CANONICAL BASIS 
# ...............
abstract type CanonicalBasis end

# OPERATIONS
# length (number of basis elements)
length( cb::T ) where {T<:CanonicalBasis} = length( cb.states )
# getindex (array interface)
function getindex( cb::T , i::Int64 ) where {T<:CanonicalBasis} 
    return cb.states[i]
end

# **********************
# SIMPLE CANONICAL BASIS
# ......................
abstract type SimpleCanonicalBasis <: CanonicalBasis end

# DISPLAY
function Base.show( io::IO, 
                    ocb::T ) where {T<:SimpleCanonicalBasis}
    for (i,s) in enumerate(ocb.states)
        println( io , "$i: $s" )
    end
end
function Base.show( io::IO, 
                    ::MIME"text/plain", 
                    ocb::T ) where {T<:SimpleCanonicalBasis}
    show( ocb )
end


# ***********************
# ORBITAL CANONICAL BASIS
# .......................
struct OrbitalCanonicalBasis <: SimpleCanonicalBasis
    states::Vector{OrbitalBasisElement}
end

# CONSTRUCTOR
function OrbitalCanonicalBasis( n::Int64 , m::Int64 )
    combs = map( x->[i for i in x] , Base.Iterators.product([1:m for _ in 1:n]...))
    states::Vector{OrbitalBasisElement} = 
        [OrbitalBasisElement(comb) 
         for comb in sort(reshape(combs,(m^n)))]
    return OrbitalCanonicalBasis( states )
end


# ********************
# SPIN CANONICAL BASIS
# ....................
#%%
struct SpinCanonicalBasis <: SimpleCanonicalBasis
    states::Vector{SpinBasisElement}
end

# CONSTRUCTOR
function SpinCanonicalBasis( n::Int64 )
    m = 2
    combs = map( x->[i for i in x] , 
                 Base.Iterators.product( [('↑','↓') for _ in 1:n]... ) )
    states::Vector{SpinBasisElement} = 
        [SpinBasisElement(comb) 
         for comb in sort(reshape(combs,(m^n)))]
    return SpinCanonicalBasis( states )
end

# *****************
# J CANONICAL BASIS
# .................
#%%
struct JCanonicalBasis <: SimpleCanonicalBasis
    states::Vector{JBasisElement}
end

# CONSTRUCTOR
function JCanonicalBasis( n::Int64 , J::Float64 )
    combs = map( x->[i for i in x] , Base.Iterators.product([-J:J for _ in 1:n]...))
    states::Vector{JBasisElement} = [
        JBasisElement(comb) 
        for comb in sort(reshape(combs,(Int64(2J+1)^n)))
    ]
    return JCanonicalBasis( states )
end

# **********************
# SPIN-0 CANONICAL BASIS
# ......................
#%%
struct SpinZeroCanonicalBasis <: SimpleCanonicalBasis
    states::Vector{SpinZeroBasisElement}
end

# CONSTRUCTOR
function SpinZeroCanonicalBasis( n::Int64 )
    states::Vector{SpinZeroBasisElement} = [SpinZeroBasisElement([nothing for _ in 1:n])]
    return SpinZeroCanonicalBasis( states )
end

# ****************************
# ORBITAL-NULL CANONICAL BASIS
# ............................
#%%
struct OrbitalNullCanonicalBasis <: SimpleCanonicalBasis
    states::Vector{OrbitalNullBasisElement}
end

# CONSTRUCTOR
function OrbitalNullCanonicalBasis( n::Int64 )
    states::Vector{OrbitalNullBasisElement} = [OrbitalNullBasisElement([nothing for _ in 1:n])]
    return OrbitalNullCanonicalBasis( states )
end

# ********************************
# COMPOUND CANONICAL BASIS GENERAL
# ................................
abstract type CompoundCanonicalBasisGeneral <: CanonicalBasis end

# ************************
# COMPOUND CANONICAL BASIS
# ........................
struct CompoundCanonicalBasis <: CompoundCanonicalBasisGeneral
    states::Vector{CompoundBasisElement}
end

# CONSTRUCTOR
function CompoundCanonicalBasis( n , m ) 
    obasis = OrbitalCanonicalBasis( n , m ) 
    sbasis = SpinCanonicalBasis( n )
    cstates = [ CompoundBasisElement(os,ss) for os in obasis.states 
                                            for ss in sbasis.states ]
    return CompoundCanonicalBasis( cstates )
end

# DISPLAY
function Base.show( io::IO, 
                    ccb::CompoundCanonicalBasis )
    for (i,s) in enumerate(ccb.states)
        println( io , "$i: $(s.orbital)$(s.spin)" )
    end
end
function Base.show( io::IO, 
                    ::MIME"text/plain", 
                    ccb::CompoundCanonicalBasis )
    show( ccb )
end

# *************************************
# DOUBLE-GROUP COMPOUND CANONICAL BASIS
# .....................................
struct DoubleGroupCompoundCanonicalBasis <: CompoundCanonicalBasisGeneral
    states::Vector{DoubleGroupCompoundBasisElement}
end

# CONSTRUCTOR
function DoubleGroupCompoundCanonicalBasis( n , m ) 
    obasis = OrbitalCanonicalBasis( n , m ) 
    sbasis = SpinZeroCanonicalBasis( n )
    cstates = [ DoubleGroupCompoundBasisElement(os,ss) for os in obasis.states 
                                            for ss in sbasis.states ]
    return DoubleGroupCompoundCanonicalBasis( cstates )
end

# DISPLAY
function Base.show( io::IO, 
                    ccb::DoubleGroupCompoundCanonicalBasis )
    for (i,s) in enumerate(ccb.states)
        println( io , "$i: $(s.orbital)$(s.spin)" )
    end
end
function Base.show( io::IO, 
                    ::MIME"text/plain", 
                    ccb::DoubleGroupCompoundCanonicalBasis )
    show( ccb )
end

# **************************
# J COMPOUND CANONICAL BASIS
# ..........................
struct JCompoundCanonicalBasis <: CompoundCanonicalBasisGeneral
    states::Vector{JCompoundBasisElement}
end

# CONSTRUCTOR
function JCompoundCanonicalBasis( n::Int64 , J::Float64 ) 
    obasis = OrbitalNullCanonicalBasis(n)
    sbasis = JCanonicalBasis(n,J)
    cstates = [ 
        JCompoundBasisElement(os,ss) 
        for os in obasis.states 
        for ss in sbasis.states 
    ]
    return JCompoundCanonicalBasis( cstates )
end

# DISPLAY
function Base.show( io::IO, 
                    ccb::JCompoundCanonicalBasis )
    for (i,s) in enumerate(ccb.states)
        println( io , "$i: $(s.orbital)$(s.spin)" )
    end
end
function Base.show( io::IO, 
                    ::MIME"text/plain", 
                    ccb::JCompoundCanonicalBasis )
    show( ccb )
end

# ######################################
# STATE
#
# - It is a Fock state represented in a 
#   CanonicalBasis as a complex vector.
#
# ######################################
mutable struct State
    vec::Vector{ComplexF64}
    basis::CanonicalBasis
end

# CONSTRUCTOR
function State( be::BE , b::CB ) where {BE<:BasisElement,CB<:CanonicalBasis}
    for (i,bs) in enumerate(b.states)
        if bs==be
            v = [ j==i ? ComplexF64(1.0) : ComplexF64(0.0) for j=1:length(b.states) ]
            return State( v , b )
        end
    end
end

# OPERATIONS
# equality
function (==)( s::State , zer::T ) where {T<:Number}
    if zer==zero(zer)
        return all(s.vec.==zero(ComplexF64))
    end
end
iszero(s::State) = iszero(s.vec)
# getindex (array interface)
getindex( s::State , i::Int64 ) = s.vec[i]
# cleaning: keep only vector components with norm above tolerance (1e-6)
function clean_state!( s::State )
    s.vec = ComplexF64[ isapprox(LinearAlgebra.norm(x),0.0,atol=1e-6) ? 0.0im : x for x in s.vec ] 
end
function clean_symstates!( symstates::Dict{ T , State } ) where {T<:Tuple}
    for k in keys(symstates) 
        clean_state!( symstates[k] ) 
    end
end
#function clean_symstates!( symstates::Dict{ Tuple{Tuple,Float64,Int64,Float64,Int64}, State } )
#    for k in keys(symstates) 
#        clean_state!( symstates[k] ) 
#    end
#end
#function clean_symstates!( symstates::Dict{ Tuple{String, Float64, Int64, Float64, Int64}, State} )
#    for k in keys(symstates) 
#        clean_state!( symstates[k] ) 
#    end
#end

# tensor multiplication of an orbital state and a spin state
function tensormult( s_o::State ,
                     s_s::State ,
                     basis_c::CCBG ) where {CCBG<:CompoundCanonicalBasisGeneral}

    s_c = State( zeros(ComplexF64,length(basis_c)) , basis_c )

    basis_o = s_o.basis
    basis_s = s_s.basis

    for (n_o,e_o) in enumerate(s_o.vec),
        (n_s,e_s) in enumerate(s_s.vec)

        reftuple = (basis_o[n_o],basis_s[n_s])

        @inbounds for i in 1:length(basis_c) 
            if basis_c[i]==reftuple
                s_c.vec[i] += e_o*e_s
                break
            end
        end
    end

    return s_c
end
#function tensormult( s_o::State ,
#                     s_s::State ,
#                     basis_c::DoubleGroupCompoundCanonicalBasis )
#
#    s_c = State( zeros(ComplexF64,length(basis_c)) , basis_c )
#
#    basis_o::OrbitalCanonicalBasis = s_o.basis
#    basis_s::SpinZeroCanonicalBasis = s_s.basis
#
#    for (n_o,e_o) in enumerate(s_o.vec),
#        (n_s,e_s) in enumerate(s_s.vec)
#
#        reftuple = (basis_o[n_o],basis_s[n_s])
#
#        @inbounds for i in 1:length(basis_c) 
#            if basis_c[i]==reftuple
#                s_c.vec[i] += e_o*e_s
#                break
#            end
#        end
#    end
#    return s_c
#end
# norm and normalization
@inline norm( s::State ) = sqrt(sum(abs.(s.vec).^2))
@inline norm( v::Vector{T} ) where {T<:Number} = sqrt(sum(abs.(v).^2))
@inline function normalize!( s::State ) 
    s.vec ./= norm(s)
end
@inline function normalize!( v::Vector{T} ) where {T<:Number} 
    v ./= norm(v) 
end

# ALGEBRA
# coefficient * state 
*( c::ComplexF64 , s::State ) = State( c.*s.vec , s.basis )
function mul!( tmp::State ,
               c::ComplexF64 ,
               s::State ) 

    tmp.vec .*= c.*s.vec

end
*( s::State , c::ComplexF64 ) = c*s
function mul!( tmp::State ,
               s::State ,
               c::ComplexF64 )

    tmp.vec .*= c.*s.vec

end

# state +/- state
+( s1::State , s2::State ) = State( s1.vec.+s2.vec , s1.basis )
-( s1::State , s2::State ) = State( s1.vec.-s2.vec , s1.basis )

# DISPLAY
function Base.show( io::IO, 
                    s::State )
    printels = []
    for (i,bs) in enumerate(s.basis.states)
        coef = s.vec[i]
        coef==0 && continue
            el = "($coef)$bs"
        push!( printels , el )
    end
    if printels==[]
        print( io , 0 )
    else
        toprint = reduce( (x,y)->x*" + "*y  , printels ) 
        print( io , toprint )
    end
end
function Base.show( io::IO, 
                    ::MIME"text/plain", 
                    s::State )
    show(s)
end


# #######################################
# PARTITIONS AND TABLEAUX
#
# - Structures and functions related with 
#   the manipulation and generation of 
#   partitions and the tableaux related 
#   with them
#
# #######################################
 
# **********
# PARTITIONS
# ..........
struct Partition 
    vec::Vector{Int64}
end
Partition( p::Tuple ) = Partition(collect(p))

# OPERATIONS
# getindex (array interface)
function getindex( partition::Partition , n::Int64 )
    return partition.vec[n]
end
# length (number of elements to permute)
length( p::Partition ) = length( p.vec )
# equality
function (==)( p1::Partition , p2::Partition )
    length(p1)!==length(p2) && return false 
    for i in 1:length(p1) 
        p1[i]!==p2[i] && return false
    end
    return true 
end
# obtain all possible partitions for n elements
function partitions( n::Int64 ; max_i::Int64=-1 )
    outer = false
    if max_i==-1
        max_i=n 
        outer = true
    end
    if n==0
        return [[]]
    end
    pp = []
    for i=1:max_i
        for p in partitions( (n-i) ; max_i=minimum((n-i,i)) )
            push!( pp , [i,p...] )
        end
    end
    if outer
        return [Partition(p) for p in pp]
    else
        return pp
    end
end
# all possible bipartitions for n elements
function bipartitions( n )
    return [ p for p in partitions(n) if length(p.vec)<=2 ]
end
# obtain the spin corresponding to a bipartition
function bipart2spin( bipart ) 
    spin::Float64 = 0.0
    if length(bipart.vec)==2
        spin = (bipart.vec[1]-bipart.vec[2])/2.0  
    else
        spin = bipart.vec[1]/2.0
    end
    return spin 
end
# get length of a given partition column
function columnlength( part::Partition , col::Int64 ) 
    return sum(map( row->row>=col ? 1 : 0 , part.vec ))
end
# get length of a given partition row
function complementary( part::Partition ) 
    vec = [columnlength(part,i) for i in 1:part[1]]
    return Partition(vec)
end

# DISPLAY
function show( io::IO , partition::Partition )
    numberstring = reduce(*,map(x->"$x,",partition.vec))
    numberstring = numberstring[1:(end-1)]
    print( io , "(" * numberstring * ")" )
end
show( io::IO, ::MIME"text/plain", p::Partition ) = show(p)
# pretty printing
function pp( partition::Partition )
    println( "·" * "---·"^partition.vec[1] )
    println( "|" * "   |"^partition.vec[1] )
    println( "·" * "---·"^partition.vec[1] )
    for i in 2:length(partition)
        println( "|" * "   |"^partition.vec[i] )
        println( "·" * "---·"^partition.vec[i] )
    end
end


# *******
# TABLEAU
# .......
#%% 
struct Tableau
    partition::Partition
    numbers::Vector{Int64}
end

# OPERATIONS
# getindex (array-like interface)
function getindex( tableau::Tableau , i::Int64 , j::Int64 )
    part::Vector{Int64} = tableau.partition.vec
    i==1 && return tableau.numbers[j]
    idx::Int64 = sum( part[ii] for ii=1:(i-1) ) + j
    return tableau.numbers[idx]
end
function getindex( tableau::Tableau , i::Int64 ) 
    return tableau.numbers[i]
end
function getindex( t::Tableau , i::Int64 , J::UnitRange{Int64} )
    return [ t[i,j] for j in J ]::Vector{Int64}
end
function getindex( t::Tableau , I::UnitRange{Int64} , j::Int64 )
    return [ t[i,j] for i in I ]
end
# size 
function size( tableau::Tableau , dim::Int64 )
    pvec = tableau.partition.vec
    dim==1 && return length(pvec)
    dim==2 && return maximum(pvec)
end
# length (number of rows) 
length( t::Tableau ) = length( t.partition )

# PARTITION=>TABLEAUX GENERATORS
# partition => tableaux
function get_tableaux( partition::Partition ) 
    maxnum::Int64 = sum( partition.vec ) 
    return [Tableau(partition,perm) 
            for perm in permutations(1:maxnum)]
end
# partition => standard tableaux
function standard_tableaux( partition::Partition )
    tableaux::Vector{Tableau} = []
    for tableau in get_tableaux(partition)
        istandard = true
        for row=1:size(tableau,1)
            for col=1:partition[row]
                if row!==1 
                    if tableau[row,col]<tableau[row-1,col]
                        istandard=false 
                        break
                    end
                end
                if col!==1 
                    if tableau[row,col]<tableau[row,col-1]
                        istandard=false 
                        break
                    end
                end
            end
        end
        istandard || continue
        push!( tableaux , deepcopy(tableau) )
    end
    return tableaux
end
# partition => partitions2tableaux 
function get_partitions2tableaux( partitions::Vector{Partition} )
    return Dict( p=>standard_tableaux(p) for p in partitions )
end

# DISPLAY
function show( io::IO , tableau::Tableau )
    accum = 1
    for part in tableau.partition.vec
        tabrow = tableau.numbers[accum:(accum+part-1)]
        for i=1:part
            print( io , "$(tabrow[i]) " )
        end
        print( io , "\n" )
        accum += part
    end
end
show( io::IO, ::MIME"text/plain", t::Tableau ) = show(t)
# pretty printing
function pp( tableau::Tableau )
    n = tableau.numbers[1]
    println( "·" * "---·"^tableau.partition.vec[1] )
    accum = 1
    for part in tableau.partition.vec 
        tabrow = tableau.numbers[accum:(accum+part-1)]
        for i=1:part 
            n = tabrow[i]
            print( "|" * " $n " )
        end
        println( "|" )
        println( "·" * "---·" ^part )
        accum += part
    end
end


# ######################################################
# GENERAL OPERATOR INTERFACE 
#
# - Operator structure and basis multiplication
#   operation for operators acting on states.
#
# - This basic interface is inherited by all
#   operator-like structs (permutation, symmetrizer...)
#
# ######################################################
abstract type Operator end

function *( o::T , s::State ) where {T<:Operator} 
    return State( o.mat*s.vec , s.basis )
end
function mul!( tmp::State , 
               o::T , 
               s::State ) where {T<:Operator} 
    LinearAlgebra.mul!( tmp.vec , o.mat , s.vec )
end


# #########################################################
# SYMMETRIC GROUP AND PERMUTATIONS 
#
# - Permutations and their action on state-like objects
#   and also on tableaux:
#       Permutation 
#
#       AbstractPermutation <: Permutation
#       MatrixPermutation   <: Permutation
#
#       SymmetricGroup
#
# - The Symmetric Group for obtaining all possible relevant 
#   permutations.
#
# #########################################################

# ***********
# PERMUTATION
# ...........
abstract type Permutation end

# ********************
# ABSTRACT PERMUTATION
# ....................
# defined abstractly, not on a basis.
# this means it can act on any object,
# regardless of whether it is defined 
# on a basis or not. it works e.g. for
# basiselement or tableau objects.
struct AbstractPermutation 
    vec::Vector{Int64}
end

# ALGEBRA
# simple basis element permutation 
function *( p::AbstractPermutation , bs::T ) where {T<:SimpleBasisElement} 
    return BasisElement( p*bs.occupations )
end
function mul!( tmp::T , p::AbstractPermutation , bs::T ) where {T<:SimpleBasisElement}
    mul!( tmp.occupations , p , bs.occupations ) 
end
# compound basis element permutation
function *( p::AbstractPermutation , bs::CBSG ) where {CBSG<:CompoundCanonicalBasisGeneral}
    return CompoundBasisElement( p*bs.orbital , p*bs.spin )
end
function mul!( tmp::CBEG ,
               p::AbstractPermutation ,
               bs::CBEG ) where {CBEG<:CompoundBasisElementGeneral}
    mul!( tmp.orbital , p , bs.orbital )
    mul!( tmp.spin    , p , bs.spin    )
end
# state permutation 
function *( p::AbstractPermutation , s::State ) 
    return State( MatrixPermutation(p,s.basis)*s.vec , s.basis )
end
function mul!( tmp::State , 
               p::AbstractPermutation ,
               s::State ) 
    mp = MatrixPermutation(p,s.basis).mat
    tmp.vec .= zero(tmp.vec[1])
    @inbounds for j=1:length(vec), i=1:length(vec)
        tmp.vec[i] .+= mp[i,j]*s.vec[j]
    end
end

# OPERATIONS
# getindex (array-like interface)
getindex( p::AbstractPermutation , n::Int64 ) = p.vec[n]
# length (number of elements to permute)
length( p::AbstractPermutation ) = length( p.vec )
#parity 
function parity( v::Vector{Int64} )
    par = 1 
    temp = copy( v )
    for i in 1:length(v)
        if temp[i]!==i 
            idx = findfirst(x->x==i,temp)
            temp[i],temp[idx] = i,temp[i]
            par *= -1
        end
    end
    return par
end
parity( p::AbstractPermutation ) = parity( p.vec )
# vector permutation (permutes by index)
function *( p::AbstractPermutation , v::Vector{T} ) where {T} 
    rv::Vector{T} = Vector{T}(undef,length(v))
    @inbounds for i=1:length(v) 
        rv[i] = v[p[i]]::T
    end
    return rv
end
function mul!( tmp::Vector{T} , 
               p::AbstractPermutation , 
               v::Vector{T} ) where T 
    @inbounds for i=1:length(v) 
        tmp[i] = v[p[i]]::T
    end
end
# tableau permutation (permutes by value)
function *( p::AbstractPermutation , t::Tableau )
    rv::Vector{Int64} = Vector{Int64}(undef,length(p.vec))
    @inbounds for i=1:length(p)
        rv[i] = p[t[i]]::Int64
    end
    return Tableau( t.partition , rv )
end
function mul!( tmp::Tableau ,
               p::AbstractPermutation , 
               t::Tableau ) 
    @inbounds for i=1:length(p) 
        tmp.numbers[i] = p[t[i]]::Int64 
    end
end

# DISPLAY
show( io::IO , p::AbstractPermutation ) = print( io , p.vec )
show( io::IO, ::MIME"text/plain", p::AbstractPermutation ) = show( p ) 

# ******************
# MATRIX PERMUTATION
# ..................
# defined concretely for a given basis.
# it acts on states defined on that basis.
struct MatrixPermutation <: Operator
    mat::Matrix{Int64}
end

# CONSTRUCTOR
function get_bs_idx( bs::BS , basis::CB ) where {BS<:BasisElement,CB<:CanonicalBasis}
    @inbounds for i=1:length(basis) 
        bs==basis[i] && return i
    end
    return nothing
end
function MatrixPermutation( ap::AbstractPermutation , cb::T ) where {T<:CanonicalBasis}
    N = length( cb )
    mat::Matrix{Int64} = zeros( Int64 , N , N )
    tmp = deepcopy(cb[1])
    @inbounds for j=1:N
        mul!( tmp , ap , cb[j] )
        idx = get_bs_idx( tmp , cb )
        mat[idx,j] = 1
    end
    return MatrixPermutation(mat)
end
function MatrixPermutation!( tmp::MatrixPermutation , ap::AbstractPermutation , cb::T ; verbose=false ) where {T<:CanonicalBasis}
    N = length( cb )
    tmp.mat .= zero(tmp.mat[1,1])
    tmpbe = deepcopy(cb[1])
    @inbounds for j=1:N
        mul!( tmpbe , ap , cb[j] )
        idx = get_bs_idx( tmpbe , cb )
        if verbose 
            @show cb[j] 
            @show tmpbe 
            println( "$j => $idx" )
        end
        tmp.mat[idx,j] = 1
    end
end

# DISPLAY
show( io::IO , p::MatrixPermutation ) = print( io , p.mat )
show( io::IO, ::MIME"text/plain", p::MatrixPermutation ) = show( p ) 

# ***************
# SYMMETRIC GROUP
# ...............
struct SymmetricGroup 
    permutations::Vector{AbstractPermutation}
end

# CONSTRUCTOR
function SymmetricGroup( n::Int64 )
    return SymmetricGroup([AbstractPermutation(p) 
                           for p in permutations(1:n)])
end

# OPERATIONS
# getindex (array-like interface)
getindex( sg::SymmetricGroup , n::Int64 ) = sg.permutations[n]
# length
length( sg::SymmetricGroup ) = length( sg.permutations )

# DISPLAY
function show( io::IO , sg::SymmetricGroup )
    n = length( sg.permutations )
    print( io , "Symmetric group of $n element permutations" )
end

# ##############################################
# ROW AND COLUMN SYMMETRIES OF A TABLEAU
#
# - Using the structures above, we can find 
#   the row- and column-symmetries of a tableau.
#
# ##############################################
# row symmetries
function get_rowsyms( t::Tableau ; verbose=false ) 
    if verbose 
        println( "GETTING ROW SYMS" )
        pp(t)
    end
    rowsyms::Vector{AbstractPermutation} = []
    p = t.partition
    n::Int64 = sum(p.vec)
    sg = SymmetricGroup(n)
    for i_sg in 1:length(sg)
        perm = sg[i_sg]
        verbose && print( "perm = $perm => " )
        transt = perm*t 
        sym = true 
        for i=1:size(t,1) 
            rowt  = sort(t[i,1:p[i]])
            rowtt = sort(transt[i,1:p[i]])
            if !(rowt==rowtt)
                sym = false 
                verbose && println( "$sym" )
                break 
            end
        end
        (verbose && sym) && println( "$sym" )
        sym && push!( rowsyms , perm )
    end
    return rowsyms
end
# column symmetries
function get_colsyms( t::Tableau ) 
    colsyms::Vector{AbstractPermutation} = []
    p = t.partition
    n::Int64 = sum(p.vec)
    sg = SymmetricGroup(n)
    for i_sg in 1:length(sg)
        perm = sg[i_sg]
        transt = perm*t 
        sym = true 
        for j=1:size(t,2) 
            colsize = maximum(collect( i for i in 1:length(p) 
                                      if p[i]>=j ))
            colt  = sort(t[1:colsize,j])
            coltt = sort(transt[1:colsize,j])
            if !(colt==coltt)
                sym = false 
                break 
            end
        end
        sym && push!( colsyms , perm )
    end
    return colsyms
end

# #######################################
# YOUNG SYMMETRIZERS
#
# - Symmetrizers from Young's theory.
#       
#       Symmetrizer 
#       
#       YoungSymmetrizer  <: Symmetrizer 
#       RowSymmetrizer    <: Symmetrizer
#       ColumnSymmetrizer <: Symmetrizer
#
# #######################################

# ***********
# SYMMETRIZER
# ...........
abstract type Symmetrizer <: Operator end

# DISPLAY
show( io::IO , sym::Symmetrizer ) = print( io , sym.mat )
show( io::IO, ::MIME"text/plain", sym::Symmetrizer ) = show( sym ) 

# ***************************
# ROW AND COLUMN SYMMETRIZERS
# ...........................
struct RowSymmetrizer <: Symmetrizer
    mat::Matrix{Int64} 
end
struct ColumnSymmetrizer <: Symmetrizer
    mat::Matrix{Int64}
end

# *****************
# YOUNG SYMMETRIZER
# .................
struct YoungSymmetrizer <: Symmetrizer
    mat::Matrix{Int64}
end

# CONSTRUCTORS
function RowSymmetrizer( t::Tableau , cb::T ) where {T<:CanonicalBasis}
    rowsyms = get_rowsyms( t::Tableau ) 
    N = length(cb)
    mat = zeros( Int64 , N , N )
    for perm in rowsyms 
        mat .+= MatrixPermutation( perm , cb ).mat
    end
    return RowSymmetrizer( mat )
end
function ColumnSymmetrizer( t::Tableau , cb::T ) where {T<:CanonicalBasis}
    colsyms = get_colsyms( t::Tableau ) 
    N = length(cb)
    mat = zeros( Int64 , N , N )
    tmp = MatrixPermutation( zeros( Int64 , N , N ) )
    for perm in colsyms 
        MatrixPermutation!( tmp , perm , cb )
        mat .+= parity(perm) .* tmp.mat
    end
    return ColumnSymmetrizer( mat )
end
function YoungSymmetrizer( t::Tableau , cb::T ) where {T<:CanonicalBasis}
    return RowSymmetrizer(t,cb)*ColumnSymmetrizer(t,cb)
end
function YoungSymmetrizer!( ys::YoungSymmetrizer , t::Tableau , cb::T ) where {T<:CanonicalBasis}
    mul!( ys , RowSymmetrizer(t,cb) , ColumnSymmetrizer(t,cb) )
end

# ALGEBRA
function *( s1::RowSymmetrizer , s2::ColumnSymmetrizer ) 
    return YoungSymmetrizer( s1.mat*s2.mat )
end
function mul!( tmp::YoungSymmetrizer , 
               s1::RowSymmetrizer , 
               s2::ColumnSymmetrizer ) 
    mul!( tmp.mat , s1.mat , s2.mat )
end
function *( s1::ColumnSymmetrizer , s2::RowSymmetrizer ) 
    return YoungSymmetrizer( s1.mat*s2.mat )
end
function mul!( tmp::YoungSymmetrizer ,
               s1::ColumnSymmetrizer ,
               s2::RowSymmetrizer ) 
    mul!( tmp.mat , s1.mat , s2.mat )
end
