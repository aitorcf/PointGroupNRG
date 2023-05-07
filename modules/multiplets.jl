#%%
using Combinatorics 
using LinearAlgebra
import Base.:show
import Base.:*
import Base.:+
import Base.:-
import Base.:(==)
import Base.:getindex
import Base.:size
import Base.:display
import Base.:length
import Base.:sort
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

# **********************
# COMPOUND BASIS ELEMENT
# ......................
struct CompoundBasisElement <: BasisElement 
    orbital::OrbitalBasisElement
    spin::SpinBasisElement
end

# OPERATIONS
# equality
(==)( be1::CompoundBasisElement , be2::CompoundBasisElement ) =
    ( be1.orbital==be2.orbital && be1.spin==be2.spin )
(==)( be_c::CompoundBasisElement , be_os::Tuple{OrbitalBasisElement,SpinBasisElement} ) =
    ( be_c.orbital==be_os[1] && be_c.spin==be_os[2] )
# order and sorting (see above)
sort( be::CompoundBasisElement ) = CompoundBasisElement( sort(be.orbital) , sort(be.spin) )
isorted( be::CompoundBasisElement ) = be.orbital==sort(be.orbital) || be.spin==sorted(be.spin)

# DISPLAY
function string( be::CompoundBasisElement )
    opart = "| " * reduce(*,map(x->"$x ",be.orbital.occupations))*")" 
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


# ************************
# COMPOUND CANONICAL BASIS
# ........................
struct CompoundCanonicalBasis <: CanonicalBasis 
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
# getindex (array interface)
getindex( s::State , i::Int64 ) = s.vec[i]
# cleaning: keep only vector components with norm above tolerance (1e-6)
function clean_state!( s::State )
    s.vec = ComplexF64[ isapprox(LinearAlgebra.norm(x),0.0,atol=1e-6) ? 0.0im : x for x in s.vec ] 
end
function clean_symstates!( symstates::Dict{Tuple{Tuple, String, Int64, Int64, Int64}, State} )
    for k in keys(symstates) 
        clean_state!( symstates[k] ) 
    end
end
function clean_symstates!( symstates::Dict{Tuple{String, Float64, Int64, Float64, Int64}, State} )
    for k in keys(symstates) 
        clean_state!( symstates[k] ) 
    end
end

# tensor multiplication of an orbital state and a spin state
function tensormult( s_o::State , 
                     s_s::State , 
                     basis_c::CompoundCanonicalBasis )

    s_c = State( zeros(ComplexF64,length(basis_c)) , basis_c )

    basis_o::OrbitalCanonicalBasis = s_o.basis
    basis_s::SpinCanonicalBasis = s_s.basis

    idx::Int64 = 0

    for (n_o,e_o) in enumerate(s_o.vec),
        (n_s,e_s) in enumerate(s_s.vec)

        reftuple = (basis_o[n_o],basis_s[n_s])

        #idx = [i for i in 1:length(basis_c)
        #         if basis_c[i]==(basis_o[n_o],basis_s[n_s])][1]
        @inbounds for i in 1:length(basis_c) 
            if basis_c[i]==reftuple
                s_c.vec[i] += e_o*e_s
                break
            end
            #basis_c[i]==(basis_o[n_o],basis_s[n_s]) && (idx=i)
        end
        #idx = findfirst( i->basis_c[i]==(basis_o[n_o],basis_s[n_s]) ,
        #                 1:length(basis_c) )

        #s_c += e_o*e_s * State( basis_c[idx] , basis_c )
        #s_c.vec[idx] += e_o*e_s
    end
    return s_c
end
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
    show( s )
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
        return [ Partition(p) for p in pp]
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
function get_partitions2tableaux( 
            partitions::Vector{Partition} )
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
function *( p::AbstractPermutation , bs::CompoundBasisElement )
    return CompoundBasisElement( p*bs.orbital , p*bs.spin )
end
function mul!( tmp::CompoundBasisElement ,
               p::AbstractPermutation ,
               bs::CompoundBasisElement )
    mul!( tmp.orbital , p , bs.orbital )
    mul!( tmp.spin    , p , bs.spin    )
end
# state permutation 
function *( p::AbstractPermutation , s::State ) 
    return State( MatrixPermutation(p,s.basis)*s.vector , s.basis )
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
        tmp.vec[i] = p[t[i]]::Int64 
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
        #for i=1:N 
        #    #mat[i,j] = dot( cb[i] , ap*cb[j] )
        #    mat[i,j] = dot( cb[i] , tmp )
        #end
    end
    return MatrixPermutation(mat)
end
function MatrixPermutation!( tmp::MatrixPermutation , ap::AbstractPermutation , cb::T ; verbose=false ) where {T<:CanonicalBasis}
    N = length( cb )
    mat::Matrix{Int64} = zeros( Int64 , N , N )
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
        #for i=1:N 
        #    #mat[i,j] = dot( cb[i] , ap*cb[j] )
        #    mat[i,j] = dot( cb[i] , tmp )
        #end
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
struct RowSymmetrizer <:Symmetrizer
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
        #mat .+= parity(perm) .* MatrixPermutation( perm , cb ).mat
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



# ############################################
# OBTENTION OF CLEBSCH-GORDAN COEFFICIENTS
# 
# - Functions for reading CG files and 
#   constructing the appropriate dictionaries.
#
# - Also some related utilities.
#
# ############################################

# get orbital dimension M
function get_M( I , cg_path ) 
    cgo = cg_orbital( I , I , cg_path )
    return maximum([k[2] for k in keys(cgo)])
end

# cg dictionary for a combination of orbitals I_1 and I_2
function cg_orbital( I_1 , I_2 , path ; verbose=false )
    # STRING version. 
    # Given two orbital irreps I_1 and I_2, it searches in path 
    # for the file containing CG information and returns it in 
    # the form of a dictionary:
    #
    #           cg[I_1,i_1,I_2,i_2,I_3,i_3] = ( I_1 , i_1 ; I_2 , i_2 | I_3 , i_3 )
    #
    cg::Dict{ Tuple{String,Int64,String,Int64,String,Int64} , ComplexF64 } = 
        Dict{ Tuple{String,Int64,String,Int64,String,Int64} , ComplexF64 }()
    file = [ x for x in readdir("$(path)/") 
               if (occursin("$(I_1)x$(I_2)",x) || occursin("$(I_2)x$(I_1)",x)) ][1]
    verbose && @show file 

    inverted = false
    occursin("$(I_2)x$(I_1)",file) && (inverted=true)

    I_3::String = "a"
    for line in readlines( "$(path)/$(file)" ) 
        line=="" && continue
        sline = split(strip(line)," ")
        I_3 = length(sline)==2 ? sline[2] : I_3
        length(sline)==2 && continue
        sline = [sline[2:3]...,sline[5],reduce(*,sline[8:end])]
        sline[end] = reduce( * , replace.( sline[end] , "I"=>"im" ) )
        sline = map( x -> eval(Meta.parse(x)) , sline )
        if ! inverted
            push!( cg , (I_1,sline[1]::Int64,I_2,sline[2]::Int64,I_3,sline[3]::Int64)=>sline[4] )
        else
            push!( cg , (I_1,sline[2]::Int64,I_2,sline[1]::Int64,I_3,sline[3]::Int64)=>sline[4] )
        end
    end
    return cg    
end

# cg dictionary for all possible combinations of orbitals that may appear in the 
# calculation
function get_cg_o_fulldict( oirreps , cg_path )
    # Given a collection of orbital irreps 'oirreps', it searches in cg_path 
    # for CG information and returns the coefficients for every possible
    # combination (I_1,I_2) for I_1 and I_2 in oirreps:
    #
    #           cg[I_1,i_1,I_2,i_2,I_3,i_3] = ( I_1 , i_1 ; I_2 , i_2 | I_3 , i_3 )
    #
    cg_o_full = Dict{ Tuple{String,Int64,String,Int64,String,Int64} , ComplexF64 }()
    for I1 in oirreps, I2 in oirreps 
        merge!( cg_o_full , cg_orbital( I1 , I2 , cg_path ) )
    end
    return cg_o_full
end

# given some initial orbital irreps, compute all orbital irreps that 
# may appear by recursive combination.
function cg_shortcircuit( CG_PATH , oirreps... ; verbose=false )
    verbose && println( "recursion call" )
    seeds::Vector{String} = collect( oirreps )
    verbose && @show seeds
    produced::Vector{String} = []
    for seed_pair in with_replacement_combinations( seeds , 2 ) 
        cg_1 = cg_orbital( seed_pair[1] , seed_pair[2] , CG_PATH ; verbose=verbose )
        cg_2 = cg_orbital( seed_pair[2] , seed_pair[1] , CG_PATH ; verbose=verbose )
        append!( produced , [k[5] for k in keys(cg_1)] )
        append!( produced , [k[5] for k in keys(cg_2)] )
        append!( produced , seed_pair )
        append!( produced , reverse(seed_pair) )
    end
    produced = collect(Set(produced))
    verbose && @show produced 
    verbose && println()
    if Set(produced)==Set(seeds)
        return produced 
    else 
        return cg_shortcircuit( CG_PATH , produced... )
    end
end

# dictionary I => dim(I)
get_oirreps2dimensions( cg_o ) = Dict( 
            k[1]=>maximum(K[2] for K in keys(cg_o) if K[1]==k[1]) 
            for k in keys(cg_o) )


# ##########################################
# SYMSTATES 
#
# - Calculation of symstates-like 
#   dictionaries. 
#
# - These functions are the main steps in
#   the calculation of antisymmetric states.
#
# ##########################################

# ********
# Ssa part
# ........
function get_dSsa2states( basis_s::SpinCanonicalBasis , 
                          bipartitions::Vector{Partition} ;
                          verbose=false )

    if verbose 
        println( "==========================================================" )
        println( "CALCULATION OF PERMUTATION- AND SPIN-SYMMETRIC SPIN STATES" )
        println( "==========================================================" )
    end

    dSsa2states = Dict{ Tuple{Tuple,Float64,Float64,Int64} , 
                        State }()
    #space = Matrix{ComplexF64}(undef,length(basis_s),0)
    space::Dict{ Tuple , Matrix{ComplexF64} } = Dict()

    for i in 1:length(basis_s)

        seed = basis_s[i]
        sortedspins = (sort(seed.occupations)...,)
        if sortedspins ∉ keys(space) 
            space[sortedspins] = Matrix{ComplexF64}(undef,length(basis_s),0) 
        end
        if verbose 
            println( "**********************" )
            println( "STATES FROM SEED $seed" )
            println( "......................" )
        end

        s = spinz( seed )
        if verbose 
            println( "S_z = $s" )
            println()
        end
        seedstate = State( seed , basis_s )

        verbose && println( "BIPARTITIONS" )
        for (bipart,tabs) in get_partitions2tableaux(bipartitions)

            S = bipart2spin( bipart )
            d = ( bipart.vec... ,) 
            verbose && @show d, S

            for (a,t) in enumerate(tabs)
                verbose && pp( t )
                ys = YoungSymmetrizer( t , basis_s )
                ts = ys*seedstate 
                if norm(ts)==zero(ComplexF64)
                    verbose && println( "state annihilated by symmetrizer" )
                    continue
                end
                if is_dependent( ts.vec , space[sortedspins] )
                    verbose && println( "linearly dependent combination" )
                    continue
                end
                space[sortedspins] = hcat( space[sortedspins] , ts.vec )
                normalize!(ts)
                sym = (d,S,s,a)
                dSsa2states[sym] = ts
                verbose && println( "$sym ==> $ts" )
            end
            verbose && println()
        end
    end
    return dSsa2states
end

# ********
# dar part
# ........
function get_dar2states( basis_o::OrbitalCanonicalBasis , 
                         partitions_o::Vector{Partition} ;
                         verbose=false )

    if verbose 
        println( "==============================================" )
        println( "COMPUTING PERMUTATION-SYMMETRIC ORBITAL STATES" )
        println( "==============================================" )
    end

    dar2states_o = Dict{ Tuple{Tuple,Int64,Int64} , State }()
    space::Dict{ Tuple , Matrix{ComplexF64} } = Dict()
    #space = Matrix{ComplexF64}(undef,length(basis_o),0)

    r = 1
    for i in 1:length(basis_o)

        seed = basis_o[i]

        # NEW
        sortednums = (sort(seed.occupations)...,) 
        if sortednums ∉ keys(space)
            space[sortednums] = Matrix{ComplexF64}(undef,length(basis_o),0)
        end

        seedstate = State( seed , basis_o )
        if verbose 
            println( "*************************" )
            println( "STATES FROM SEED $seed" )
            println( ".........................." )
        end

        nonzero = false
        saturated = false
        for (part,tabs) in get_partitions2tableaux(partitions_o)
            
            d = ( part.vec... ,) 
            verbose && @show d

            for (a,t) in enumerate(tabs)
                verbose && pp( t )
                ys = YoungSymmetrizer( t , basis_o )
                ts = ys*seedstate 
                if norm(ts)==zero(ComplexF64)
                    verbose && println( "state annihilated by symmetrizer" )
                    continue
                end
                nonzero = true
                if is_dependent( ts.vec , space[sortednums] )
                    verbose && println( "linearly dependent" )
                    continue
                end
                #space = hcat( space , ts.vec )
                space[sortednums] = hcat( space[sortednums] , 
                                          ts.vec )
                normalize!(ts)
                sym = (d,a,r)
                dar2states_o[sym] = ts
                verbose && println( "$sym ==> $ts" )
                #if size(space,2)==length(basis_o)
                #    println( "space saturated. returning." )
                #    return dar2states_o
                #end
            end
            verbose && println()
        end
        nonzero && (r+=1)
    end
    return dar2states_o
end

# ***************
# Iir2states part
# ...............
function get_Iir2states( 
            N::Int64 , 
            M::Int64 , 
            orbital::String , 
            cg_path::String ;
            verbose=false )

    basis_o = OrbitalCanonicalBasis(N,M)
    oirreps = cg_shortcircuit( cg_path , orbital )
    cg_o = get_cg_o_fulldict( oirreps , cg_path )
    oirreps2dimensions = get_oirreps2dimensions( cg_o ) 
    symstates_block::Dict{ Tuple{String,Int64,Int64} , Vector{Tuple{ComplexF64,Vector{Int64}}} }=
        Dict( (orbital,i,1)=>[(1.0im,[i])] for i=1:oirreps2dimensions[orbital] ) 
    symstates_add::Dict{ Tuple{String,Int64,Int64} , Vector{Tuple{ComplexF64,Vector{Int64}}} }=
        symstates_block

    for n in 2:N
        symstates_new = combine_symstates( symstates_add ,
                                           symstates_block ,
                                           cg_o ,
                                           oirreps2dimensions ;
                                           verbose=verbose )
        symstates_block = symstates_new
        if verbose 
            for (k,vs) in symstates_block
                ss = reduce( + , [s[1]*State(OrbitalBasisElement(s[2]),OrbitalCanonicalBasis(n,M)) for s in vs] )
                verbose && @show ss
            end
        end
    end

    Iir2states::Dict{ Tuple{String,Int64,Int64} , State } = Dict()
    for (k,vs) in symstates_block
        Iir2states[k] = reduce( + , [s[1]*State(OrbitalBasisElement(s[2]),basis_o) for s in vs] )
    end

    return Iir2states
end

function combine_symstates( symstates_block::Dict{ Tuple{String,Int64,Int64} , Vector{Tuple{ComplexF64,Vector{Int64}}} },
                            symstates_add::Dict{ Tuple{String,Int64,Int64} , Vector{Tuple{ComplexF64,Vector{Int64}}} },
                            cg_o::Dict{ Tuple{String,Int64,String,Int64,String,Int64} , ComplexF64 },
                            oirreps2dimensions;
                            verbose=false )

    verbose && println( "COMBINING SYMSTATES..." )

    symnew = Dict{ Tuple{String,Int64,Int64} , Vector{Tuple{ComplexF64,Vector{Int64}}} }()

    multiplets_a = Set((k[1],k[3]) for k in keys(symstates_add))
    multiplets_b = Set((k[1],k[3]) for k in keys(symstates_block))

    for m_a in multiplets_a, m_b in multiplets_b 

        if verbose 
            println("decomposition of multiplets m_a=$m_a and m_b=$m_b") 
        end

        I_a, I_b = m_a[1], m_b[1]
        r_a, r_b = m_a[2], m_b[2]
        II_c = Set(k[5] for k in keys(cg_o) 
                        if (k[1],k[3])==(I_a,I_b))

        for I_c in II_c 

            r = get_r( symnew , I_c )
            verbose && @show r
            m_c = (I_c,r)

            verbose && println( "m_c=$m_c" )

            for i_c in 1:oirreps2dimensions[I_c],
                i_a in 1:oirreps2dimensions[I_a],
                i_b in 1:oirreps2dimensions[I_b]

                
                cg_key = (I_a,i_a,I_b,i_b,I_c,i_c) 
                coeff = get( cg_o , cg_key ,
                             zero(ComplexF64) )
                coeff==zero(ComplexF64) && continue

                verbose && println("$cg_key => $coeff") 

                vs_a = symstates_add[I_a,i_a,r_a]
                vs_b = symstates_block[I_b,i_b,r_b]
                vs_c = merge_vecstates(vs_a,vs_b,coeff)
                
                if (I_c,i_c,r) in keys(symnew)
                    append!( symnew[I_c,i_c,r] , vs_c )
                else 
                    symnew[I_c,i_c,r] = vs_c
                end
            end
        end
    end

    symnew = Dict( k=>fuse_states(v) 
                   for (k,v) in symnew )
    if verbose
        println( "SYMSTATES:" )
        for (k,v) in symnew 
            @show k 
            @show v 
            println()
        end
    end
    return symnew
end

function fuse_states( vs::Vector{Tuple{ComplexF64,Vector{Int64}}} ) 
    basisels = Set( x[2] for x in vs ) 
    fused::Vector{Tuple{ComplexF64,Vector{Int64}}} = []
    for be in basisels 
        coeff = reduce( + , [x[1] for x in vs if x[2]==be] )::ComplexF64
        push!( fused , (coeff,be) )
    end
    return fused
end

@inline function get_r( symstates_new::Dict{ Tuple{String,Int64,Int64} , Vector{Tuple{ComplexF64,Vector{Int64}}} }, 
                        I::String )
    if I in Set( k[1] for k in keys(symstates_new) )
        r::Int64 = maximum(Set(kk[3] for kk in keys(symstates_new) if kk[1]==I))
        r += 1
        return r
    else
        return 1::Int64
    end
end

function merge_vecstates( vs_add::Vector{Tuple{ComplexF64,Vector{Int64}}}  , 
                          vs_block::Vector{Tuple{ComplexF64,Vector{Int64}}}  , 
                          coeff::ComplexF64=1.0 )
    vs_new::Vector{Tuple{ComplexF64,Vector{Int64}}} = []
    for e_add in vs_add, e_block in vs_block 
        c = coeff*e_add[1]*e_block[1] 
        s = vcat(e_add[2],e_block[2])
        push!( vs_new , (c,s) ) 
    end
    return vs_new
end

# *********
# dIair part
# .........
function get_dIair2states( dar2states::Dict{ Tuple{Tuple,Int64,Int64} , State },
                           Iir2states::Dict{ Tuple{String,Int64,Int64} , State } ,
                           basis_o::OrbitalCanonicalBasis ;
                           verbose=false )

    dIair2states::Dict{ Tuple{Tuple,String,Int64,Int64,Int64} , State } = Dict()

    for (d,a) in Set(k[1:2] for k in keys(dar2states)),
        (I,i) in Set(k[1:2] for k in keys(Iir2states))

        S_dar = hcat([ s.vec for (k,s) in dar2states 
                             if k[1:2]==(d,a) ]...)
        S_Iir = hcat([ s.vec for (k,s) in Iir2states 
                             if k[1:2]==(I,i) ]...)

        int = subspace_intersection( S_dar , S_Iir )
        size(int,2)==0 && continue

        for r in 1:size(int,2) 
            dIair2states[d,I,a,i,r] = State(int[:,r],basis_o)
        end
        verbose && println()
    end

    return dIair2states
end

# ************
# dSsdIar part
# ............
function get_dSsadIair2states( dIair2states_o::Dict{ Tuple{Tuple,String,Int64,Int64,Int64} , State } ,
                               dSsa2states_s::Dict{ Tuple{Tuple,Float64,Float64,Int64} , State } ,
                               basis_c::CompoundCanonicalBasis ;
                               verbose=false )

    verbose && println( "GETTING dSsadIair2states" )
    dSsadIair2states = Dict{ Tuple{Tuple,Float64,Float64,Int64,Tuple,String,Int64,Int64,Int64} , State }()

    for (k_o,s_o) in dIair2states_o,
        (k_s,s_s) in dSsa2states_s 

        p_o = Partition([k_o[1]...])
        p_s = Partition([k_s[1]...])
        if verbose
            pp(p_o)
            pp(p_s)
        end
        p_o==complementary(p_s) || continue

        #if verbose
        #    pp(p_o)
        #    pp(p_s)
        #end

        k_c = (k_s...,k_o...)
        s_c = tensormult( s_o , s_s , basis_c )
        dSsadIair2states[k_c] = s_c

        if verbose
            println( "INGREDIENTS AND RESULT" )
            println( "$k_o => $s_o" )
            println( "$k_s => $s_s" )
            @show s_c
            println()
        end

    end
    return dSsadIair2states
end


# ******************
# antisymmetric part 
# ..................
function get_M( basis_c::CompoundCanonicalBasis ) 
    return basis_c[length(basis_c)].orbital.occupations[1]*2
end
function get_N( basis_c::CompoundCanonicalBasis )
    return length(basis_c[1].orbital.occupations) 
end
function get_asymdim( basis_c::CompoundCanonicalBasis ) 
    M = get_M( basis_c ) 
    N = get_N( basis_c )
    res::Int64 = factorial(M)/factorial(M-N)/factorial(N)
    return res
end
function get_asymsubspace( N::Int64 , 
                           basis_c::CompoundCanonicalBasis ;
                           verbose=false )
    if verbose 
        println( "================================" )
        println( "COMPUTING ANTISYMMETRIC SUBSPACE" )
        println( "================================" )
    end
    asymdim = get_asymdim(basis_c)
    ap = Partition([1 for _ in 1:N])
    at = Tableau( ap , collect(1:N) )
    if verbose 
        @show asymdim 
        pp(at)
    end
    ay = YoungSymmetrizer( at , basis_c )
    ts = State( basis_c[1] , basis_c )
    asymsubspace = Matrix{ComplexF64}(undef,length(basis_c),asymdim)
    #if (asymdim==1 && N>0)  # particle saturation
    #    if verbose 
    #        println( "********************" )
    #        println( "shortcut calculation" )
    #        println( "--------------------" )
    #    end
    #    be = BasisElement( [1,2,3,1,2,3] , ['↑','↓','↑','↓','↑','↓'] )
    #    verbose && @show bs
    #    s = State( be , basis_c )
    #    mul!( ts , ys , s ) 
    #    verbose && @show ts
    #    normalize!( ts )
    #    asymsubspace[:,1] .= ts.vec 
    #    return asymsubspace 
    #end
    idx = 1
    for i in 1:length(basis_c)
        s = State( basis_c[i] , basis_c ) 
        #ts = ay*s 
        mul!( ts , ay , s )
        if verbose
            @show s 
            @show ts 
        end
        if ts==0 
            verbose && println()
            continue 
        end
        if idx==1
            verbose && println( "first state in" )
            verbose && println()
            normalize!(ts)
            #asymsubspace = hcat( asymsubspace , ts.vec )
            asymsubspace[:,idx] .= ts.vec
            idx += 1
            verbose && @show idx 
            verbose && println()
            idx>asymdim && break
        else
            if is_dependent( ts.vec , asymsubspace[:,1:(idx-1)] )
                if verbose 
                    println( "dependent state" )
                    println()
                end
                continue
            end
            verbose && println( "new state is independent. including..." )
            normalize!(ts)
            #asymsubspace = hcat( asymsubspace , ts.vec )
            asymsubspace[:,idx] .= ts.vec
            idx += 1
            verbose && @show idx
            verbose && println()
            idx>asymdim && break
        end
    end
    if verbose 
        println( "result:" )
        pp(asymsubspace)
        verbose && println()
    end
    return asymsubspace
end
#function get_asymsubspace( N::Int64 , 
#                           basis_c::CompoundCanonicalBasis ;
#                           verbose=false )
#    if verbose 
#        println( "================================" )
#        println( "COMPUTING ANTISYMMETRIC SUBSPACE" )
#        println( "================================" )
#    end
#    ap = Partition([1 for _ in 1:N])
#    at = Tableau( ap , [i for i in 1:N] )
#    ay = YoungSymmetrizer( at , basis_c )
#    asymsubspace = Matrix{ComplexF64}(undef,length(basis_c),0)
#    for i in 1:length(basis_c)
#        s = State( basis_c[i] , basis_c ) 
#        ts = ay*s 
#        if verbose
#            @show s 
#            @show ts 
#        end
#        if ts==0 
#            verbose && println()
#            continue 
#        end
#        if size(asymsubspace,2)==0
#            verbose && println( "first state in" )
#            verbose && println()
#            normalize!(ts)
#            asymsubspace = hcat( asymsubspace , ts.vec )
#            continue
#        end
#        is_dependent( ts.vec , asymsubspace ) && continue
#        verbose && println( "new state is independent. including..." )
#        normalize!(ts)
#        asymsubspace = hcat( asymsubspace , ts.vec )
#        verbose && println()
#    end
#    if verbose 
#        println( "result:" )
#        pp(asymsubspace)
#        verbose && println()
#    end
#    @show size(asymsubspace)
#    return asymsubspace
#end

# **********
# ISisr part
# ..........
function get_asym_ISisr( dSsadIair2states::Dict{ Tuple{Tuple,Float64,Float64,Int64,Tuple,String,Int64,Int64,Int64} , State } ,
                         asymsubspace::Matrix{ComplexF64} ,
                         basis_c::CompoundCanonicalBasis ;
                         verbose=false )

    if verbose
        println( "===========================" )
        println( "ANTISYMMETRIC INTERSECTIONS" )
        println( "===========================" )
    end

    ISisr = Dict{ Tuple{String,Float64,Int64,Float64,Int64} , State }()
    for sym in Set( (k[1:3]...,k[5],k[6],k[8]) 
                    for k in keys(dSsadIair2states) )

        if verbose
            @show sym
            println( "***********************"^2 )
        end

        (_,S,s,_,I,i) = sym

        sub = hcat([s.vec for (k,s) in dSsadIair2states 
                    if (k[1:3]==sym[1:3] && (k[5],k[6],k[8])==sym[4:end])]...)
        
        int = subspace_intersection( asymsubspace , sub ;
                                     verbose=false)

        asymsym = [State(int[:,i],basis_c) for i in 1:size(int,2)]
        size(int,2)==0 && continue

        if verbose
            for as in asymsym 
                @show as 
            end
            println()
        end

        for (r,as) in enumerate(asymsym)
            ISisr[I,S,i,s,r] = as 
        end
    end
    return ISisr
end

# #############################################
#
# LINEAR INDEPENDENCE AND SUBSPACE INTERSECTION 
#
# #############################################
function is_dependent( v::Vector{T} , S::Matrix{T} ; atol::Float64=1e-6 ) where {T<:Number}
    return any([ isapprox(0.0,s,atol=atol) for s in svd(hcat(S,v)).S ])
end

function gram_schmidt!( sub::Matrix{T} ) where {T<:Number}
    temp = similar(sub[:,1])
    @inbounds for j in 1:size(sub,2) 
        temp = sub[:,j]
        j==1 || (temp .-= sum([ (sub[:,j]'*sub[:,jj])*sub[:,jj] for jj=1:(j-1) ]))
        sub[:,j] .= temp ./ norm(temp)
    end
    return sub
end

function subspace_intersection( S1::Matrix{T} , 
                                S2::Matrix{T} , 
                                atol::Float64=1e-6 ;
                                verbose=false ) where {T<:Number}
    int = S1 * nullspace(hcat(S1,-S2),atol=atol)[1:size(S1,2),:] 
    size( int , 2 )==0 && return int
    gram_schmidt!(int)
    return int
end

function pp( m::Matrix ) 
    for i in 1:size(m,1)
        println( m[i,:] )
    end
end


# ##############################################
# AUTOMATIZATION 
#
# - Functions that compute all the antisymmetric 
#   states with minimal input.
#
# ##############################################
function compute_asymstates_allN( orbital::String ,
                                  cg_path::String ,
                                  multiplets_path::String ;
                                  verbose::Bool=false )

    # add convention name to multiplet folder
    multiplets_path *= "/$(orbital)_julia"

    # create asym dir if it does not exist
    isdir(multiplets_path) || mkdir( multiplets_path )

    # compute multiplet states
    for n in 2:2*get_M(orbital,cg_path)
        compute_asymstates_N( orbital , n , cg_path , multiplets_path ; verbose=verbose )
    end
end

function compute_asymstates_N( 
            orbital::String , 
            N::Int64 ,
            cg_path::String ,
            asym_path::String ;
            verbose=false ,
            identityrep="" )

    # orbital dimensions
    M = get_M( orbital , cg_path )


    # compound basis 
    basis_c = CompoundCanonicalBasis( N , M )

    # write down antisymmetric basis 
    asymsubspace = get_asymsubspace( N , basis_c , verbose=verbose )
    open( "$(asym_path)/N$(N)basis.txt" , "w" ) do io
        for i in 1:size(asymsubspace,2) 

            v = asymsubspace[:,i]
            idx = findfirst( x->(!isapprox(abs(x),0.0)) , v )
            basis_element = basis_c.states[idx]

            o_part = basis_element.orbital.occupations
            s_part = basis_element.spin.occupations

            o_string = reduce( (x,y)->"$x $y" , o_part )
            s_string = reduce( (x,y)->"$x $y" , s_part )

            write( io , "$o_string $s_string \n" )

        end
    end

    if (N==M*2 && identityrep!=="")

        if verbose 
            println( "************************************************" )
            println( "SHORTCUT COMPUTATION FOR THE SATURATED CASE N=2M" )
            println( "------------------------------------------------" )
        end

        s = State( asymsubspace[:,1] ) 
        sym = (identityrep,0.0,1,0.0,1)
        ISisr = Dict( sym => s )

        filename = "$(asym_path)/N$N.txt" 
        cmd = `touch $filename` 
        run(cmd)
        open( "$(asym_path)/N$N.txt" , "w" ) do io 
            for (k,s) in ISisr 
                symstring = reduce( * , map(x->"$x ",[k...]) )
                canrep = s.vec 
                canstring = reduce( * , map(x->"($x)  ",canrep) ) 
                asymrep = collect(flatten(nullspace(hcat(asymsubspace,-canrep),atol=1e-6)[1:size(asymsubspace,2),:]))
                normalize!(asymrep)
                asymstring = reduce( * , map(x->"($x)  ",asymrep) ) 
                toprint = reduce( (x,y)->x*"| "*y , [symstring,canstring,asymstring] )
                verbose && println( toprint )
                println( io , toprint )
            end
        end
        return ISisr
    end



    # *********
    # SPIN PART
    # .........

    if verbose
        println( "# · ######### · #" )
        println( "# | --------- | #" )
        println( "# | SPIN PART | #" )
        println( "# | --------- | #" )
        println( "# · ######### · #" )
        println()
    end

    #%% BASIS
    basis_s = SpinCanonicalBasis(N)
    if verbose 
        println( "==========" )
        println( "SPIN BASIS" )
        println( "==========" )
        println( basis_s )
        println()
    end

    #%% BIPARTITIONS
    partitions_s = bipartitions(N)
    if verbose 
        println( "=======================================" )
        println( "SPIN VALUES AND ASSOCIATED BIPARTITIONS" )
        println( "=======================================" )
        for p in partitions_s
            println( "S = $(bipart2spin(p)) ==> partition: $(p.vec)" )
            pp( p )
            println()
        end
    end

    #%% spin symstates
    dSsa2states_s = get_dSsa2states( basis_s , partitions_s ;
                                     verbose=verbose )
    if verbose
        println( "===========================================" )
        println( "SPIN- AND PERMUTATION-SYMMETRIC SPIN STATES" )
        println( "===========================================" )
        for (k,s) in dSsa2states_s 
            @show k 
            @show s
            println()
        end
    end


    # ************
    # ORBITAL PART 
    # ............
    if verbose
        println( "# · ############ · #" )
        println( "# | ------------ | #" )
        println( "# | ORBITAL PART | #" )
        println( "# | ------------ | #" )
        println( "# · ############ · #" )
        println()
    end

    #%% orbital basis
    basis_o = OrbitalCanonicalBasis( N , M )
    if verbose 
        println( "=============" )
        println( "ORBITAL BASIS" )
        println( "=============" )
        show( basis_o )
        println()
    end

    #%% partitions
    partitions_so = Dict( p=>complementary(p) for p in partitions_s )
    partitions_o = [complementary(p) for p in partitions_s]
    if verbose 
        println( "==========================================" )
        println( "ORBITAL PARTITIONS (complementary to spin)" )
        println( "==========================================" )
        for (sp,op) in partitions_so 
            println( sp , " ==> " , op )
            pp( op )
            println()
        end
    end

    #%% permutation-symmetric orbital states 
    dar2states_o = get_dar2states( basis_o , 
                                   partitions_o ; 
                                   verbose=verbose )
    if verbose 
        println( "====================================" )
        println( "PERMUTATION-SYMMETRIC ORBITAL STATES" )
        println( "====================================" )
        for (k,v) in dar2states_o 
            println( "$k ==> $v" )
            println()
        end
    end


    #%% orbital-symmetric orbital states
    Iir2states_o = get_Iir2states( N , M , orbital , cg_path ; verbose=verbose )
    if verbose 
        println( "================================" )
        println( "ORBITAL-SYMMETRIC ORBITAL STATES" )
        println( "================================" )
        for (k,v) in Iir2states_o
            println( "$k ==> $v" )
            println()
        end
    end


    #%% permutation and orbital symmetric states
    dIair2states_o = get_dIair2states( dar2states_o , 
                                       Iir2states_o ,
                                       basis_o ;
                                       verbose=verbose )
    if verbose
        println( "=========================================" )
        println( "ORBITAL- AND PERMUTATION-SYMMETRIC STATES" )
        println( "=========================================" )
        clean_symstates!( dIair2states_o )
        for (k,v) in dIair2states_o
            println( "$k ==> $v" ) 
        end
        println()
    end

    
    # *************
    # COMBINED PART
    # .............
    if verbose 
        println( "# · ############# · #" )
        println( "# | ------------- | #" )
        println( "# | COMBINED PART | #" )
        println( "# | -------.----- | #" )
        println( "# · ############# · #" )
        println()
    end

    #%% compound basis 
    basis_c = CompoundCanonicalBasis( N , M )
    if verbose 
        println( "==============" )
        println( "COMPOUND BASIS" )
        println( "==============" )
        show(basis_c)
    end

    #%% independently symmetric states
    dSsadIair2states = get_dSsadIair2states( dIair2states_o , 
                                             dSsa2states_s ,
                                             basis_c ;
                                             verbose=verbose )
    if verbose
        println( "================" )
        println( "dSsadIair STATES" )
        println( "================" )
        for (k,v) in dSsadIair2states 
            @show k 
            @show v 
            println()
        end
    end

    # ******************
    # ANTISYMMETRIC PART
    # ..................
    #%% completely antisymmetric states
    asymsubspace = get_asymsubspace( N , basis_c , verbose=verbose )

    # *************************************
    # SYMMETRY-ADAPTED ANTISYMMETRIC STATES
    # .....................................
    #%%
    ISisr = get_asym_ISisr( dSsadIair2states , 
                            asymsubspace ,
                            basis_c ;
                            verbose=verbose )
    clean_symstates!( ISisr )

    # *******
    # WRITING
    # .......
    filename = "$(asym_path)/N$N.txt" 
    cmd = `touch $filename` 
    run(cmd)
    open( "$(asym_path)/N$N.txt" , "w" ) do io 
        for (k,s) in ISisr 
            symstring = reduce( * , map(x->"$x ",[k...]) )
            canrep = s.vec 
            canstring = reduce( * , map(x->"($x)  ",canrep) ) 
            asymrep = collect(Base.Iterators.flatten(nullspace(hcat(asymsubspace,-canrep),atol=1e-6)[1:size(asymsubspace,2),:]))
            for e in asymrep 
                if abs(e)!==0.0 
                    if real(e)<0.0
                        asymrep = -asymrep
                    end
                    break 
                end
            end
            normalize!(asymrep)
            asymstring = reduce( * , map(x->"($x)  ",asymrep) ) 
            toprint = reduce( (x,y)->x*"| "*y , [symstring,canstring,asymstring] )
            verbose && println( toprint )
            println( io , toprint )
        end
    end
          
    return ISisr
end
