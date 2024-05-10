using PartialWaveFunctions

# **********************************************
#
# This module includes functions that deal with 
# the symmetry-related part of the calculations.
# This includes generation of symmetry-adapted 
# states, manipulation of dictionaries classi-
# fied according to symmetry, and Clebsh-Gordan 
# coefficients, among others.
#
# **********************************************

const ClearIrrep   = Tuple{Int64,String,Float64}
const ClearPartner = Tuple{Int64,Float64}
const ClearQNums   = Tuple{Int64,String,Float64,Int64,Float64,Int64}
const IntIrrep   = NTuple{3,Int64}
const IntPartner = NTuple{2,Int64}
const IntQNums   = NTuple{6,Int64}

const AbstractIrrep = Union{ClearIrrep,IntIrrep}
const AbstractPartner = Union{ClearPartner,IntPartner}
const AbstractQNums = Union{ClearQNums,IntQNums}

const ClearTripleQ = NTuple{3,ClearQNums}
const IntTripleQ   = NTuple{3,IntQNums}
const AbstractTripleQ = Union{ClearTripleQ,IntTripleQ}
const ClearTripleG = NTuple{3,ClearIrrep}
const IntTripleG   = NTuple{3,IntIrrep}
const AbstractTripleG = Union{ClearTripleG,IntTripleG}

const ClearSymstateDict = Dict{ ClearQNums , State }
const IntSymstateDict   = Dict{ IntQNums , State }
const AbstractSymstateDict = Union{ClearSymstateDict,IntSymstateDict}

const ClearQPCG = Dict{ ClearTripleQ , ComplexF64 }
const IntQPCG = Dict{ IntTripleQ , ComplexF64 }
const AbstractQPCG = Union{ClearQPCG,IntQPCG}

const ClearIrrepPCG = Dict{ ClearTripleG , Array{ComplexF64,3} }
const IntIrrepPCG   = Dict{ IntTripleG ,   Array{ComplexF64,3} }
const IntIrrepPCGNS = Dict{ IntTripleG ,   Array{ComplexF64,4} }
const AbstractIrrepPCG = Union{ClearIrrepPCG,IntIrrepPCG}

const ClearMultiplet = Tuple{Int64,String,Float64,Int64}
const IntMultiplet   = NTuple{4,Int64}
const AbstractMultiplet = Union{ClearMultiplet,IntMultiplet}

const ClearMultipletSet = Set{ClearMultiplet}
const IntMultipletSet = Set{IntMultiplet}
const AbstractMultipletSet = Union{ClearMultipletSet,IntMultipletSet}
const ClearMultipletVector = Vector{ClearMultiplet}
const IntMultipletVector = Vector{IntMultiplet}
const AbstractMultipletVector = Union{ClearMultipletVector,IntMultipletVector}

const ClearIrrepSet = Set{ClearIrrep}
const IntIrrepSet = Set{IntIrrep}
const AbstractIrrepSet = Union{ClearIrrepSet,IntIrrepSet}
const ClearIrrepVector = Vector{ClearIrrep}
const IntIrrepVector = Vector{IntIrrep}
const AbstractIrrepVector = Union{ClearIrrepVector,IntIrrepVector}

const IntCGKey = NTuple{3,Int64}
const IntCG = Dict{ IntCGKey , Array{ComplexF64,3} }
const IntCGNS = Dict{ IntCGKey , Array{ComplexF64,4} } # NS=NonSimple

# ##################
# ORBITAL GENERATORS 
# ..................

function get_atomic_states( 
            external_label::Int64 ,
            ostring::String , 
            oirreps2dimensions::Dict{String,Int64} )::Vector{Tuple{Int64,String,Int64,String}}

    root = ( external_label , ostring )
    statetuples = Tuple{Int64,String,Int64,String}[]

    for orbital in 1:oirreps2dimensions[ostring], 
        spin in ["u","d"]

        push!( statetuples , (root...,orbital,spin) )

    end

    return statetuples
end

function Eg_states( shell::Int64 )::Vector{Tuple{Int64,String,Int64,String}}
    # input:
    # - shell: number of the shell for which to 
    #          compute the states (usually=0)
    # output:
    # - statuples: a vector of tuples, each of 
    #              characterizes a state
    root = ( shell , "e" )
    statetuples::Vector{Tuple{Int64,String,Int64,String}} = []
    for orbital in 1:2, spin in ["u","d"]
        push!( statetuples , (root...,orbital,spin) )
    end
    return statetuples
end

function T2g_states( shell::Int64 )::Vector{Tuple{Int64,String,Int64,String}}
    # input:
    # - shell: number of the shell for which to 
    #          compute the states (usually=0)
    # output:
    # - statuples: a vector of tuples, each of 
    #              characterizes a state
    root = ( shell , "t" )
    statetuples::Vector{Tuple{Int64,String,Int64,String}} = []
    for orbital in 1:3, spin in ["u","d"]
        push!( statetuples , (root...,orbital,spin) )
    end
    return statetuples
end

# C4v
function A1B1_states( shell::Int64 )::Vector{Tuple{Int64,String,Int64,String}} 
    # input:
    # - shell: number of the shell for which to 
    #          compute the states (usually=0)
    # output:
    # - statuples: a vector of tuples, each of 
    #              characterizes a state
    roots = [( shell , "a" ),( shell , "b" )]
    statuples::Vector{Tuple{Int64,String,Int64,String}} = []
    for root in roots, spin in ["u","d"] 
        push!( statuples , (root...,1,spin) ) 
    end
    return statuples 
end
function A1_states( shell::Int64 )::Vector{Tuple{Int64,String,Int64,String}}
    # input:
    # - shell: number of the shell for which to 
    #          compute the states (usually=0)
    # output:
    # - statuples: a vector of tuples, each of 
    #              characterizes a state
    root = ( shell , "a" )
    statetuples::Vector{Tuple{Int64,String,Int64,String}} = []
    for spin in ["u","d"]
        push!( statetuples , (root...,1,spin) )
    end
    return statetuples
end
function B1_states( shell::Int64 )::Vector{Tuple{Int64,String,Int64,String}}
    # input:
    # - shell: number of the shell for which to 
    #          compute the states (usually=0)
    # output:
    # - statuples: a vector of tuples, each of 
    #              characterizes a state
    root = ( shell , "b" )
    statetuples::Vector{Tuple{Int64,String,Int64,String}} = []
    for spin in ["u","d"]
        push!( statetuples , (root...,1,spin) )
    end
    return statetuples
end

function shell_hilbertstates( orbital::String , r::Int64 , dimension::Int64 )::Vector{Tuple{Int64,String,Int64,String}}
    root = ( r , "orbital" )
    statetuples::Vector{Tuple{Int64,String,Int64,String}} = [] 
    for spin in ["u","d"], gamma in 1:dimension 
        push!( statetuples , (root...,gamma,spin) )
    end
    return statetuples
end


# ##########
# ASYMSTATES
# ..........

function read_asymfile_fockbasis( filename::String ; verbose=false )::Dict{ Vector{String} , Vector{ComplexF64} }
    # read N3.txt kind of file and return a dictionary
    #       qnums => vector
    d::Dict{ Vector{String} , Vector{ComplexF64} } = Dict()
    for line in readlines(filename)
        verbose && @show line
        line[1]=='#' && continue
        sym,can,fock = (strip.(split(line,"|"))...,)
        verbose && @show sym,can,fock
        k = split( sym , " " )
        fockels = split( fock , "  " ) 
        verbose && @show fockels
        fockv = map( x -> eval(Meta.parse(x)) , fockels )
        get!( d , k , fockv )
    end
    return d
end

# OLD, for sage-generated files
function read_asymfile_fockbasis_sage( filename::String )
    # read N3.txt kind of file and return a dictionary
    #       qnums => vector
    d = Dict()
    for line in readlines(filename)
        line[1]=='#' && continue
        get!( d , split(line," ")[1:5] , 
                  map( x -> eval(Meta.parse(x)) , split(strip(split(line,"|")[3])," ") ))
    end
    return d
end



# #################################################
# HELPER FUNCTIONS FOR DICTIONARIES AND SETS
#
# - Obtain information about "symstates" type
#   variable or any other dictionary with structure
#           ( qnums => anything )
# .................................................

function print_dict( dict )
    for (k,v) in dict 
        println( k )
        println( v )
        println()
    end
end

function get_multiplets( symstates::Dict{T,S} ) where {T<:Tuple,S<:State}
    return Set( (k[1],k[2],k[3],k[end]) for k in keys(symstates) )
end

function get_irreps( symstates::Dict{T,S} ; multiplicity=false )::ClearIrrepSet where {T<:Tuple,S<:State}
    # gets irreps (with multiplicities, in tuple form) 
    # of "symstates" or other Dict(q=>any) object.
    irreps::ClearIrrepSet = Set( k[1:3] for k in keys(symstates) ) 
    if multiplicity
        return Set( (irr,get_multiplicity(symstates,irr)) for irr in irreps )
    else
        return irreps
    end
end
function get_irreps( multiplets::Set{T} ; multiplicity=false ) where {T<:Tuple}
    # gets irreps (with multiplicities, in tuple form) 
    # of "symstates" or other Dict(q=>any) object.
    irreps = Set( m[1:3] for m in multiplets ) 
    if multiplicity
        return Set( (irr,get_multiplicity(multiplets,irr)) for irr in irreps )
    else
        return irreps
    end
end

function get_irrepstates( symstates::Dict{TS,S} , irrep::TI ) where {TS<:Tuple,S<:State,TI<:Tuple}
    # returns symstates belonging the specified irrep
    return Dict{ TS , S }( q=>s for (q,s) in symstates if q[1:3]==irrep )
end

function get_irrepstates_onepartner( symstates::Dict{TS,S} , irrep::TI ; partner="o1_smax" ) where {TS<:Tuple,S<:State,TI<:Tuple}
    # return symstates corresponding to a chosen partner 
    # of the specified irrep 
    if partner=="o1_smax"
        return Dict{ TS , S }( q=>s for (q,s) in get_irrepstates(symstates,irrep) if (q[4],q[5])==(1,q[3]) )
    end
end

function get_partners( symstates::Dict{TS,S} ; irreps=false , multiplicity=false ) where {TS<:Tuple,S<:State}
    # gets partners (with multiplicities and/or irreps, in tuple form) 
    # of "symstates" or other Dict(q=>any) object.
    partners = Set( k[1:5] for k in keys(symstates) )
    if ( multiplicity && !irreps )
        return Set( (p,get_multiplicity(symstates,p[1:3])) for p in partners )
    elseif ( irreps && !multiplicity )
        return Set( (p,p[1:3]) for p in partners )
    elseif ( irreps && multiplicity )
        return Set( (p,p[1:3],get_multiplicity(symstates,p[1:3]))
                    for p in partners )
    else
        return partners
    end
end

function get_partners_of_irrep( symstates::Dict{TS,S} , irrep::TI ) where {TS<:Tuple,S<:State,TI<:Tuple} 
    return Set( k[1:5] for k in keys(symstates) if k[1:3]==irrep )
end

function get_multiplicity( symstates::Dict{TS,S} , irrep::TI ) where {TS<:Tuple,S<:State,TI<:Tuple} 
    return length(collect(m for m in get_multiplets(symstates) if m[1:3]==irrep))
end
function get_multiplicity( multiplets::Set{TM} , irrep::TI )::Int64 where {TM<:Tuple,TI<:Tuple}
    return sum(map( m->( (m[1],m[2],m[3])==irrep ? 1 : 0 ) , collect(multiplets) ))
end


# ###########################
# CLEBSCH GORDAN COEFFICIENTS
# ...........................

function cg_orbital( I_1::String , 
                     I_2::String , 
                     path::String ; 
                     verbose::Bool=false )
    # STRING version. 
    # Given two orbital irreps I_1 and I_2, it searches in path 
    # for the file containing CG information and returns it in 
    # the form of a dictionary:
    #
    #           cg[I_1,i_1,I_2,i_2,I_3,i_3] = ( I_1 , i_1 ; I_2 , i_2 | I_3 , i_3 )
    #
    cg::Dict{ Tuple{String,Int64,String,Int64,String,Int64} , ComplexF64 } = Dict()

    file = [ x for x in readdir( "$(path)/" ) 
             if (occursin("_$(I_1)x$(I_2).",x) || occursin("_$(I_2)x$(I_1).",x)) ][1]
    verbose && @show file 

    inverted = I_1!==I_2 && occursin("$(I_2)x$(I_1)",file)

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

function get_cg_o_fulldict( 
            oirreps::Vector{String} , 
            cg_path::String )::Dict{Tuple{String,Int64,String,Int64,String,Int64}}
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

function get_cg_o_fullmat( cg_o_fulldict::Dict{Tuple{String,Int64,String,Int64,String,Int64},ComplexF64} )::Dict{Tuple{String,String,String},Array{ComplexF64,3}}
    # Given a CG coefficient dictionary in the form 
    #
    #           cg[I_1,i_1,I_2,i_2,I_3,i_3] = ( I_1 , i_1 ; I_2 , i_2 | I_3 , i_3 ),
    #
    # it transforms it into a dictionary of matrices 
    #
    #           cg[I_1,I_2,I_3] = M(I_1,I_2,I_3), where
    #           M(I_1,I_2;I_3)_{i_1,i_2,i_3}=(I_1,i_1;I_2,i_2|I_3,i_3),
    #
    oirreps = Set( k[1] for k in keys(cg_o_fulldict) )
    combs = Set( (k[1],k[3],k[5]) for k in keys(cg_o_fulldict) )
    cg_o_fullmat::Dict{Tuple{String,String,String},Array{ComplexF64,3}} = Dict() 
    for (I1,I2,I3) in combs
        D1 = maximum(Set( k[2] for k in keys(cg_o_fulldict) if k[1]==I1 ))
        D2 = maximum(Set( k[4] for k in keys(cg_o_fulldict) if k[3]==I2 ))
        D3 = maximum(Set( k[6] for k in keys(cg_o_fulldict) if k[5]==I3 ))
        push!( cg_o_fullmat , 
               (I1,I2,I3) => 
               ComplexF64[ get(cg_o_full,(I1,i1,I2,i2,I3,i3),zero(ComplexF64))
                                         for i1 in 1:D1, i2 in 1:D2, i3 in 1:D3 ]
         )
    end
    return cg_o_fullmat 
end

function get_cg_o_fullmatint( 
            cg_o_fulldict::Dict{Tuple{String,Int64,String,Int64,String,Int64},ComplexF64} , 
            oirreps::Vector{String} 
            )::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} } 
    # 
    # Given a CG coefficient dictionary in the form 
    #
    #           cg[I_1,i_1,I_2,i_2,I_3,i_3] = ( I_1 , i_1 ; I_2 , i_2 | I_3 , i_3 ),
    #
    # where I_1, I_2 and I_3 are STRINGS, it transforms it into a dictionary of matrices 
    #
    #           cg[I_1,I_2,I_3] = M(I_1,I_2,I_3), where
    #           M(I_1,I_2;I_3)_{i_1,i_2,i_3}=(I_1,i_1;I_2,i_2|I_3,i_3),
    #
    # where all indices are INTEGERS.
    #
    oirreps2indices = Dict( o=>i for (i,o) in enumerate(oirreps) )
    combs = Set( (k[1],k[3],k[5]) for k in keys(cg_o_fulldict) )
    cg_o_fullmat::Dict{Tuple{Int64,Int64,Int64},Array{ComplexF64,3}} = Dict() 
    for (Is1,Is2,Is3) in combs
        D1 = maximum(Set( k[2] for k in keys(cg_o_fulldict) if k[1]==Is1 ))
        D2 = maximum(Set( k[4] for k in keys(cg_o_fulldict) if k[3]==Is2 ))
        D3 = maximum(Set( k[6] for k in keys(cg_o_fulldict) if k[5]==Is3 ))
        I1 = oirreps2indices[Is1]
        I2 = oirreps2indices[Is2]
        I3 = oirreps2indices[Is3]
        push!( cg_o_fullmat , 
               (I1,I2,I3) => 
               ComplexF64[ get(cg_o_fulldict,(Is1,i1,Is2,i2,Is3,i3),zero(ComplexF64))
                                         for i1=1:D1, i2=1:D2, i3=1:D3 ]
         )
    end
    return cg_o_fullmat 
end

function cg_spin( 
            S_1::Float64 , 
            S_2::Float64 
            )::Dict{ NTuple{6,Float64} , ComplexF64 }

    cg_s::Dict{ NTuple{6,Float64} , ComplexF64 } = Dict()

    cp::NTuple{6,Float64} = (0,0,0,0,0,0)

    S2_1::Int64 = convert(Int64,2*S_1)
    S2_2::Int64 = convert(Int64,2*S_2)

    @inbounds for S2_3::Int64=abs(S2_1-S2_2):2:(S2_1+S2_2)
        @inbounds for s2_1::Int64=(-S2_1):2:S2_1, 
                      s2_2::Int64=(-S2_2):2:S2_2, 
                      s2_3::Int64=(-S2_3):2:S2_3
                    cp = map( x->x/2.0 , (S2_1,s2_1,S2_2,s2_2,S2_3,s2_3) )
                    push!( cg_s , 
                        cp=>PartialWaveFunctions.clebschgordan_doublearg(S2_1,s2_1,S2_2,s2_2,S2_3,s2_3) )
        end
    end

    return cg_s
end


function cg_shortcircuit( CG_PATH::String , 
                          oirreps::Vector{String} ; 
                          verbose=false )::Vector{String}

    verbose && println( "CG SHORTCIRCUIT" )
    seeds::Vector{String} = collect( oirreps )
    verbose && @show seeds
    produced::Vector{String} = []
    for seed_pair::Vector{String} in with_replacement_combinations( seeds , 2 ) 
        cg_1 = cg_orbital( seed_pair[1] , seed_pair[2] , CG_PATH ; verbose=verbose )
        cg_2 = cg_orbital( seed_pair[2] , seed_pair[1] , CG_PATH ; verbose=verbose )
        append!( produced , [k[5] for k in keys(cg_1)] )
        append!( produced , [k[5] for k in keys(cg_2)] )
        append!( produced , seed_pair )
        append!( produced , reverse(seed_pair) )
    end
    produced = collect(Set(produced))
    verbose && @show produced 
    println()
    if Set(produced)==Set(seeds)
        return produced 
    else 
        return cg_shortcircuit( CG_PATH , produced )
    end
end

function cg_combination( reduced_combination::NTuple{2,Tuple{Int64,String,Float64}} , 
                         oh_path::String )
    # Used by cg_reduce_product_states().
    # Given a reduces_combination=( (N_1,I_1,S_1) , (N_2,I_2,S_2) ) 
    # in STRING format, it returns a tuple (cg_o,cg_s) with a 
    # CG dictionary (see above) in STRING format containing 
    # the information in a format that can be used by 
    # cg_reduces_product_states().
    (I_0::String,S_0::Float64) = reduced_combination[1][2:end]
    (I_1::String,S_1::Float64) = reduced_combination[2][2:end]
    cg_o::Dict{ Tuple{String,Int64,String,Int64,String,Int64} , ComplexF64 }  = cg_orbital( I_0::String , I_1::String , oh_path::String )
    cg_s::Dict{ NTuple{6,Float64} , ComplexF64 } = cg_spin( S_0::Float64 , S_1::Float64 )
    return tuple(cg_o,cg_s)::Tuple{Dict{ Tuple{String,Int64,String,Int64,String,Int64} , ComplexF64 },Dict{ NTuple{6,Float64} , ComplexF64 }}   
end

# ###################################
# SYMSTATE CONSTRUCTION AND UTILITIES
# ...................................

function symstates_n_oneirrep( 
            basis::CanonicalBasis{SFS} , 
            n::Int64 , 
            hiztegia::D,
            symstates::Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , S }=Dict() ;
            filename::String="" ,
            identity::String="A1g" )::Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , S } where {SFS<:SymbolFockState,D<:Dict,S<:State}
    # input:
    # - canonical basis for N=n subspace 
    # - number of particles N=n 
    # - dictionary "hiztegia" of: state label => standard symmetry label 
    # - dictionary "symstates" of labelled symmetric states
    # output: # - symstates actualized
    if n==0  
        get!( symstates , 
              (0,identity,0.0,1,0.0,1) ,
              GenericBasis(basis).states[1] ) 
        return symstates
    elseif n==1 
        for sfs::SFS in basis.states
            (g::Int64,I::String,mu::Int64,m::String) = sfs.hilbert.states[ findfirst(sfs.occ) ]
            get!( symstates , 
                  (1,hiztegia[I],m=="-" ? 0.0 : 0.5,mu,hiztegia[m]::Float64,g) , 
                  State(sfs,basis) )
        end
        return symstates
    end
    fileinfo::Dict{ Vector{String} , Vector{ComplexF64} } = read_asymfile_fockbasis( filename )
    for (sym,v) in fileinfo
        s::State{CanonicalBasis{SFS}} = State( v , basis )
        ssym::Tuple{String,Float64,Int64,Float64,Int64} = ( sym[1] , map( x->eval(Meta.parse(x)) , sym[2:end] )... )
        get!( symstates , (n,ssym...) , s )
    end
    return symstates
end

#function cg_reduce_product_states( symstates_1 , symstates_2 , oh_path )
#    # take symmetric states symstates_1 and symstates_2, which 
#    # usually belong to two different hilbert spaces and thus
#    # have two different bases, and generate a new set of symmetric 
#    # states defined on a combined basis using clebsch-gordan 
#    # coefficients. 
#    # input
#    # - symstates_1 : symmetric states (will appear on the left) 
#    # - symstates_2 : symmetric states (will appear on the right)
#    # - cg : for each irrep combination to be reduced (redcomb), 
#    #        it contains cg coefficients for the orbital and 
#    #        spin parts --> cg[redcomb][1] is orbital,
#    #                       cg[redcomb][2] is spin.
#    # output
#    # - symstates_new : symmetric states containing states from 
#    #                   both subspaces.
#
#    irrep_combinations = Iterators.product( get_multiplets(symstates_1) ,
#                                            get_multiplets(symstates_2) )    
#
#    symstates = Dict()
#    for ircomb in irrep_combinations
#
#        (N_1,I_1,S_1,r_1) = ircomb[1]
#        (N_2,I_2,S_2,r_2) = ircomb[2]
#
#        redcomb = ( (N_1,I_1,S_1) , (N_2,I_2,S_2) )
#        N_3 = N_1 + N_2
#        r_3 = (r_1,r_2)
#
#        cg_comb = cg_combination( redcomb , oh_path )
#
#        cg_o = cg_comb[1]
#        cg_s = cg_comb[2]
#
#        for (q_o,c_o) in cg_o, (q_s,c_s) in cg_s 
#            I_3 = q_o[end-1]
#            S_3 = q_s[end-1]
#            mu_3 = q_o[end]
#            m_3 = q_s[end]
#
#            q_3 = ((N_1,N_2,N_3),(I_1,I_2,I_3),(S_1,S_2,S_3),mu_3,m_3,r_3)
#
#            s1 = symstates_1[( N_1 , I_1 , S_1 , q_o[2] , q_s[2] , r_1 )]
#            s2 = symstates_2[( N_2 , I_2 , S_2 , q_o[4] , q_s[4] , r_2 )]
#            merge!( + , symstates , Dict(q_3 => c_o*c_s*x(s1,s2)) )
#        end
#    end
#    print_symstates_ordered( symstates )
#
#    symstates_new = Dict()
#    principals = Dict()
#    for (qnums,state) in symstates
#        if all(state.vector.==0)
#            println( qnums )
#        end
#        N_3 = qnums[1][3]
#        I_3 = qnums[2][3]
#        S_3 = qnums[3][3]
#        mu_3 = qnums[4]
#        m_3 = qnums[5]
#
#        p = ( N_3 , I_3 , S_3 , mu_3 , m_3 )
#        merge!( + , principals , Dict(p => 1) )
#        q_3 = ( p... , principals[p] )
#
#        get!( symstates_new , q_3 , state )
#    end
#    return symstates_new
#end
function cg_reduce_product_states( symstates_1::Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , S } , 
                                   symstates_2::Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , S } ,
                                   oh_path::String )::Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , S } where {S<:State}
    # take symmetric states symstates_1 and symstates_2, which 
    # usually belong to two different hilbert spaces and thus
    # have two different bases, and generate a new set of symmetric 
    # states defined on a combined basis using clebsch-gordan 
    # coefficients. 
    # input
    # - symstates_1 : symmetric states (will appear on the left) 
    # - symstates_2 : symmetric states (will appear on the right)
    # - cg : for each irrep combination to be reduced (redcomb), 
    #        it contains cg coefficients for the orbital and 
    #        spin parts --> cg[redcomb][1] is orbital,
    #                       cg[redcomb][2] is spin.
    # output
    # - symstates_new : symmetric states containing states from 
    #                   both subspaces.

    irrep_combinations = Iterators.product( get_multiplets(symstates_1) ,
                                            get_multiplets(symstates_2) )    

    symstates = Dict()
    for ircomb in irrep_combinations

        (N_1::Int64,I_1::String,S_1::Float64,r_1::Int64) = ircomb[1]
        (N_2::Int64,I_2::String,S_2::Float64,r_2::Int64) = ircomb[2]

        redcomb::NTuple{2,Tuple{Int64,String,Float64}} = ( (N_1,I_1,S_1) , (N_2,I_2,S_2) )
        N_3::Int64 = N_1 + N_2
        r_3::NTuple{2,Int64} = (r_1,r_2)

        cg_comb::Tuple{Dict{ Tuple{String,Int64,String,Int64,String,Int64} , ComplexF64 },Dict{ NTuple{6,Float64} , ComplexF64 }} = 
                cg_combination( redcomb , oh_path )

        cg_o::Dict{ Tuple{String,Int64,String,Int64,String,Int64} , ComplexF64 } = cg_comb[1]
        cg_s::Dict{ NTuple{6,Float64} , ComplexF64 } = cg_comb[2]

        for (q_o::Tuple{String,Int64,String,Int64,String,Int64},c_o::ComplexF64) in cg_o, 
            (q_s::NTuple{6,Float64},c_s::ComplexF64) in cg_s 

            I_3::String = q_o[end-1]
            S_3::Float64 = q_s[end-1]
            mu_3::Int64 = q_o[end]
            m_3::Float64 = q_s[end]

            q_3::Tuple{ NTuple{3,Int64} , NTuple{3,String} , NTuple{3,Float64} , Int64 , Float64 , NTuple{2,Int64} } = 
                ((N_1,N_2,N_3),(I_1,I_2,I_3),(S_1,S_2,S_3),mu_3,m_3,r_3)

            s1::S = symstates_1[( N_1 , I_1 , S_1 , q_o[2] , q_s[2] , r_1 )]
            s2::S = symstates_2[( N_2 , I_2 , S_2 , q_o[4] , q_s[4] , r_2 )]
            merge!( + , symstates , Dict(q_3 => c_o*c_s*x(s1,s2)) )
        end
    end

    symstates_new::Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , S } = Dict()
    principals::Dict{ Tuple{Int64,String,Float64,Int64,Float64} , Int64 } = Dict()
    sortedkeys = collect(keys(symstates))
    sort!( sortedkeys , by=x->x[1] )
    sort!( sortedkeys , by=x->x[2] )
    sort!( sortedkeys , by=x->x[3] )
    sort!( sortedkeys , by=x->x[end] )

    for qnums in sortedkeys

        state::State = symstates[qnums]

        N_3::Int64   = qnums[1][3]
        I_3::String  = qnums[2][3]
        S_3::Float64 = qnums[3][3]
        mu_3::Int64  = qnums[4]
        m_3::Float64 = qnums[5]

        p::Tuple{Int64,String,Float64,Int64,Float64} = ( N_3 , I_3 , S_3 , mu_3 , m_3 )
        merge!( + , principals , Dict(p => 1) )

        q_3::Tuple{Int64,String,Float64,Int64,Float64,Int64} = ( p... , principals[p] )

        get!( symstates_new , q_3 , state )
    end

    return symstates_new
end

#function read_basis( hilbert , irrepath )
#
#    N = length(basis.hilbert.states)
#    r = hilbert.states[1][1]
#    I = hilbert.states[1][2]
#
#    hizt = Dict( 
#            "↑"=>"u",
#            "↓"=>"d" 
#    )
#
#    for i in 0:1 
#        basis_number = [CanonicalBasis(hilbert,i) for i=0:length(hilbert.states)]
#    end
#    for i in 2:N
#        open( "$(irrepath)N$(i)basis.txt" , "r" ) do io 
#            for line in io
#                o_s = strip.(split(line," "))
#                halfline = (length(o_s)/2)
#                o = o_s[1:halfline]
#                s = [hizt[x] for x in o_s[(halfline+1):(2*halfline)]]
#                for i in 1:halfline
#                    state = (r,I,o[i],s[i])
#                    stateindex = findfirst( x->x==state , hilbert.states )
#                end
#
#
#            end
#        end
#    end
#
#end
function oneirrep_symstates( hilbert::HS , 
                             hiztegia::D ,
                             identityrep::String , 
                             irrepath::String ) where {HS<:HilbertSpace,D<:Dict}

    basis = CanonicalBasis( hilbert )
    basis_number = [CanonicalBasis(hilbert,i) for i=0:length(hilbert.states)]

    symstates::Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , State } = Dict()
    for (i::Int64,b) in enumerate(basis_number)
        n = i-1
        symstates = symstates_n_oneirrep( b , 
                                          n , 
                                          hiztegia , 
                                          symstates ;
                                          filename=irrepath*"N$n.txt" , 
                                          identity=identityrep )
    end
    for (k,v) in symstates 
        symstates[k] = extend_basis( v , basis )
    end

    return symstates
end


function symstates_irreps_and_partners( symstates )
    irr_part = Set( q[1:(end-1)] for q in keys(symstates) ) 
    irrs     = Set( q[1:3] for q in irr_part )
    return (irr_part, irrs)
end

function irrpart_from_irr( irr , irrs_parts )
    for irr_part in irrs_parts 
        if irr_part[1:3]==irr
            return irr_part 
        end
    end
end


# ###############################
# SYMMETRIC OPERATOR CONSTRUCTORS
# ...............................

function pickrandom( d )
    for (k,v) in d 
        return (k,v)
    end
end

function symstates_n( symstates , n )
    syms = Dict()
    for (k,v) in symstates
        k[1]!==n && continue
        merge!( syms , Dict(k[2:end]=>v) )
    end
    return syms 
end

function sfs2coperator( sfs::SFS , 
                        basis::CB ) where {SFS<:SymbolFockState,CB<:CanonicalBasis}
    # returns the creation operator that creates the state 
    # represented by the symbolic state sfs when applied to 
    # the vacuum state |)
    snames = [sfs.hilbert.states[i] for i=1:length(sfs.hilbert) 
              if sfs.occ[i]==1 ]
    coperator = Operator( 1 , basis )
    for sname in snames
        sop = SymbolCreationOperator( sname )
        coperator *= Operator( sop , basis )
    end
    return coperator
end

function basis2coperators( basis ; n=-1 )
    coperators = []
    for (i,sfs) in enumerate(basis.states)
        if n==-1
            push!( coperators , sfs2coperator(sfs,basis) )
        else
            if sum(sfs.occ)==n
                push!( coperators , sfs2coperator(sfs,basis) )
            end
        end
    end
    return coperators 
end


function check_atomic( s::State ; verbose=false )
    basis = s.basis 
    sfs_components = [ basis.states[i] for i=1:length(basis.states)
                       if ! isapprox(s.vector[i],0) ]
    if verbose
        for sfs in sfs_components 
            println( sfs )
        end
    end
    sfs_occupations= [ [sfs.hilbert.states[i] for i=1:length(sfs.hilbert.states) 
                        if !isapprox(sfs.occ[i],0)]
                        for sfs in sfs_components ]
    sfs_occupations= collect(Iterators.flatten( sfs_occupations ))
    verbose && println( sfs_occupations )
    any([ o[1]!==0 for o in sfs_occupations ]) && (return false)
    return true 
end


function epsilon_sym( symstates , symparams::Dict{Tuple{String,Int64},ComplexF64} ; verbose=false )
    symstates_n1 = symstates_n( symstates , 1 )
    basis = pickrandom( symstates_n1 )[2].basis
    n = length( basis.states )
    coperators = basis2coperators( basis )
    epsop = Operator( 0 , basis )
    for (q,s) in symstates_n1
        verbose && println( s )
        # construct the corresponding part of the operator
        for (i,component) in enumerate(s.vector)
            component==0 && continue
            cop = coperators[i]
            epsop += symparams[(q[1],q[5])]*component*cop*adjoint(cop)
        end
    end
    return epsop
end
function epsilon_sym( symstates , symparams::Dict{String,Vector{ComplexF64}} ; verbose=false )
    # transform to previous convention (above)
    symparams = Dict( (k,r)=>v for (k,V) in symparams for (r,v) in enumerate(V) )
    symstates_n1 = symstates_n( symstates , 1 )
    basis = pickrandom( symstates_n1 )[2].basis
    n = length( basis.states )
    coperators = basis2coperators( basis )
    epsop = Operator( 0 , basis )
    for (q,s) in symstates_n1
        verbose && println( s )
        # construct the corresponding part of the operator
        for (i,component) in enumerate(s.vector)
            component==0 && continue
            cop = coperators[i]
            epsop += symparams[(q[1],q[5])]*component*cop*adjoint(cop)
        end
    end
    return epsop
end

function electron_counter_sym( symstates , 
                               multiplet::Tuple{String,Int64} ; 
                               verbose=false )

    # gather one-particle symstates
    symstates_n1 = symstates_n( symstates , 1 )

    # basis
    basis = pickrandom( symstates_n1 )[2].basis
    n = length( basis.states )

    # creation operators for all the one-electron states
    coperators = basis2coperators( basis )

    # initialize counting operator
    counter = Operator( 0 , basis )

    # select symstate from selected multiplet
    symstates_chosen = [ s for (q,s) in symstates if (q[1],q[2],q[end])==(1,multiplet...) ]

    # find the non-zero (=1) component of the symstates in the
    # and create the number operator for the corresponding site.
    for symstate in symstates_chosen,
        (i,component) in enumerate(symstate.vector)

        # filter components
        isapprox(component,zero(component)) && continue
        @assert isapprox(component,1.0) "Problem with one-electron state occupation: $component."

        # add contribution to counter operator
        cop = coperators[i]
        counter += cop*adjoint(cop)

    end

    return counter
end

function u_sym( symstates , symparams ; verbose=false )
    if verbose 
        println( "CONSTRUCTING U..." )
        println()
    end
    symstates_n2 = symstates_n( symstates , 2 )
    basis = pickrandom( symstates_n2 )[2].basis 
    coperators = basis2coperators( basis )
    u = Operator( 0 , basis )
    # u = sum_q1q2 create_q1 annihilate_q2
    for (q1,s1) in symstates_n2, (q2,s2) in symstates_n2 
        # last index is multiplicity 
        (q1[1:(end-1)]!==q2[1:(end-1)]) && continue 
        if verbose
            println( "*******" )
            println( "symmetric states" )
            print( q1 , " ==> " )
            println( s1 )
            print( q2 , " ==> " )
            println( s2 )
            println()
        end
        for (i1,c1) in enumerate(s1.vector), (i2,c2) in enumerate(s2.vector)
            (c1==0 || c2==0) && continue 
            if verbose 
                println( "coefficients: c1=$c1, c2=$c2" )
                println( "states" )
                println( basis.states[i1] )
                println( basis.states[i2] )
            end
            creop = coperators[i1]
            annop = adjoint( coperators[i2] )
            c = symparams[q1[1:2]][q1[end],q2[end]]*c1*conj(c2)
            u += c * creop*annop
            verbose && println( "( sfs1 | u | sfs2 ) = $c" )
        end
        verbose && println()
    end
    return u
end

function gamma_sym( symstates , symparams , inshell )
    symstates_n1 = symstates_n( symstates , 1 )
    basis = pickrandom( symstates_n1 )[2].basis
    n = length( basis.states )
    coperators = basis2coperators( basis )
    epsop = Operator( 0 , basis )
    for (q,s) in symstates_n1
        # check whether only atomic states 
        # are occupied
        !check_atomic( s ) && continue
        # construct the corresponding part of the operator
        for (i,component) in enumerate(s.vector)
            component==0 && continue
            cop = coperators[i]
            epsop += symparams[q[1]]*component*cop*adjoint(cop)
        end
    end
    return epsop
end

function irrep_and_energy( gb::GenericBasis , symstates , energies ;
                             verbose=false )
    symsubs = Set([ q[1:3] for q in keys(symstates) ])
    subvecs = Dict()
    symenergies = []
    for sub in symsubs
        vecs = hcat([s.vector for (q,s) in symstates 
                     if q[1:3]==sub]...)
        merge!( subvecs , Dict(sub=>vecs) )
    end
    for (i,s) in enumerate(gb.states)
        verbose && println( s )
        for (sub,vecs) in subvecs 
            all_vecs = hcat( vecs , s.vector )
            U,S,V = svd( all_vecs )
            if any( isapprox.(S,0;atol=1e-4) ) 
                push!( symenergies , [sub,energies[i]] )
                if verbose
                    println( "$sub ==> $(energies[i])" )
                    println()
                end
            end
        end
    end
    return symenergies
end

function test_parameters( eps_range , eps_N , u_range , u_N ,
                          symstates ; target=0 )
    Deps = eps_range[2] - eps_range[1]
    Du   = u_range[2] - u_range[1]
    eps_vals = [ (eps_range[1] + Deps/eps_N*i) for i=0:eps_N ]
    u_vals   = [ (u_range[1] + Du/u_N*i) for i=0:u_N ]

    symstates_n1 = symstates_n( symstates , 1 )
    symstates_n2 = symstates_n( symstates , 2 )

    basis = pickrandom( symstates )[2].basis
    coperators = basis2coperators( basis )

    for eps_e in eps_vals, eps_b in eps_vals
        for u_a1 in u_vals, u_b1 in u_vals, u_a2 in u_vals,
            u_e1 in u_vals, u_b2 in u_vals, u_e0 in u_vals 
            epsilon_symparams = Dict(
                "E" => eps_e,
                "B1"=> eps_b
            )
            epsilon = epsilon_sym(symstates_n1,epsilon_symparams,coperators)
            u_symparams = Dict(
                ("A1",0) => [u_a1 0; 0 u_a1],
                ("A2",1) => [u_a2][:,:],
                ("B1",0) => [u_b1][:,:],
                ("B2",0) => [u_b2][:,:],
                ("E",0)  => [u_e0][:,:],
                ("E",1)  => [u_e1][:,:]
            )
            u = u_sym(symstates,u_symparams;verbose=false)
            H = epsilon + u; H = 0.5*(H + H')
            (e,u) = diagonalize( H )
            eigenbasis = GenericBasis( basis , u )
            symenergy = energy_by_symmetry( eigenbasis , symstates , e )
            if (symenergy[1][1]==target || target==0)
                println( "eps_e = $eps_e ; eps_b = $eps_b" )
                println( "u_a1 = $u_a1 ; u_a2 = $u_a2 ; u_b1 = $u_b1" )
                println( "u_b2 = $u_b2 ; u_e0 = $u_e0 ; u_e1 = $u_e1" )
                println( "******************************************" )
                for (i,syme) in enumerate(symenergy)
                    i>5 && break 
                    println( syme[1] , " ==> " , syme[2] )
                end
                println()
                println()
            end
        end
    end
end


# ###############################################
# SYMMETRIC DIAGONALIZATION OF ATOMIC HAMILTONIAN
# ...............................................

function symdiag( irreps_0::ClearIrrepSet , symstates_0::SSD , H::O ) where {SSD<:ClearSymstateDict,O<:Operator}
    irrEU_imp = Dict{ ClearIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} }()
    for G in irreps_0 
        #size = Int64( oirreps2dimensions[G[2]] * (2*G[3]+1) )
        Gstates = get_irrepstates_onepartner( symstates_0 , G ; partner="o1_smax" )
        N = length( Gstates )
        hmat = Matrix{ComplexF64}( undef , N , N )
        for r1=1:N, r2=1:N
            s1 = Gstates[( G... , 1 , G[3] , r1 )]
            s2 = Gstates[( G... , 1 , G[3] , r2 )]
            hmat[r1,r2] = s1 * H * s2
        end
        hmat = 0.5*(hmat + conj(hmat'))
        (e,u) = eigen( hmat )
        e = real.(e)
        push!( irrEU_imp , G=>(e,u) )
    end
    return irrEU_imp
end

get_irrEU_clean( identityrrep ) = 
    Dict{ Tuple{Int64,String,Float64} , Tuple{Vector{Float64},Matrix{ComplexF64}} }( 
        (0,identityrrep,0)=>( Float64[0.0] , ComplexF64[1.0im][:,:] ),
    )


# ###################
# SPIN CLEBSCH-GORDAN
# ...................
function get_cg_s_fullmatint( max_spin2::Int64 )

    cg_s_fullmatint = Dict{ NTuple{3,Int64} , Array{ComplexF64,3} }()

    for S1=0:1:max_spin2,
        S2=0:1:max_spin2,
        S3=abs(S1-S2):2:(S1+S2)

        dS1 = S1 + 1
        dS2 = S2 + 1
        dS3 = S3 + 1

        cg_s_fullmatint[S1,S2,S3] = zeros( ComplexF64 , dS1 , dS2 , dS3 )
        cgsview = @view cg_s_fullmatint[S1,S2,S3][:,:,:]

        for (i1,s1) in enumerate((-S1):2:S1),
            (i2,s2) in enumerate((-S2):2:S2),
            (i3,s3) in enumerate((-S3):2:S3)

            cgsview[i1,i2,i3] = CG_doublearg(S1,s1,S2,s2,S3,s3) 

        end
    end
    
    return cg_s_fullmatint

end

# #############
# MISCELLANEOUS
# .............

get_oirreps2indices( oirreps ) = Dict( o=>i for (i,o) in enumerate(oirreps) )

function print_multiplets_Nordered( multiplets ) 
    mm = collect(multiplets) 
    sort!( mm , by=x->(1000*x[1]+10*x[3]+x[4]) ) 
    for m in mm
        println( m )
    end
    println()
end

function print_symstates_Nordered( symstates ) 
    ss = collect(symstates) 
    sort!( ss , by=x->(100_000*x[1][1]
                      +1000*x[1][3]
                      +10*x[1][5]
                      +x[1][6]) ) 
    for (q,s) in ss
        @show q 
        println( s ) 
        println()
    end
end


