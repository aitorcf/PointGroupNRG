# ========= #
# UTILITIES #
# ========= #

# get symmetry type
ispointspin(symmetry::String) = symmetry=="PS" || symmetry=="pointspin"
isdoublegroup(symmetry::String) = symmetry=="D" || symmetry=="doublegroup"
istotalangularmomentum(symmetry::String) = symmetry=="J" || symmetry=="totalangularmomentum" || symmetry=="S" || symmetry=="spin"
isorbital(symmetry::String) = ispointspin(symmetry) || isdoublegroup(symmetry)

# total angular momentum
doublej(J::Float64) = Int64(2J)
jdim(J::Float64) = Int64(2J+1)
jdim(J::Int64) = J+1

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
               if (occursin("_$(I_1)x$(I_2).",x) || occursin("_$(I_2)x$(I_1).",x)) ][1]
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
    for k in keys(cg_o) 
)

function cg_orbital_nonsimple( I_1 , I_2 , path ; verbose=false )

    # STRING version. 
    # Given two orbital irreps I_1 and I_2, it searches in path 
    # for the file containing CG information and returns it in 
    # the form of a dictionary:
    #
    #           cg[I_1,i_1,I_2,i_2,I_3,i_3] = ( I_1 , i_1 ; I_2 , i_2 | I_3 , i_3 )
    #
    #
    cg::Dict{ Tuple{String,Int64,String,Int64,String,Int64,Int64} , ComplexF64 } = 
        Dict{ Tuple{String,Int64,String,Int64,String,Int64,Int64} , ComplexF64 }()
    file = [ x for x in readdir("$(path)/") 
               if (occursin("_$(I_1)x$(I_2).",x) || occursin("_$(I_2)x$(I_1).",x)) ][1]
    verbose && @show file 

    inverted = (I_2!==I_1 && occursin("$(I_2)x$(I_1)",file))

    I_3::String = "a"
    r_3::Int64 = 1

    sline1::String = ""
    sline2::String = ""
    sline3::String = ""
    sline4::String = ""

    i1::Int64 = 0
    i2::Int64 = 0
    i3::Int64 = 0
    c::ComplexF64 = 0.0

    # iterate over lines in file
    for sline::Vector{String} in split.(strip.(eachline("$(path)/$(file)"),' '))

        # skip empty lines
        length(sline)==0 && continue

        # header line --> assign irrep and multiplicity, then continue
        I_3,r_3 = length(sline)==3 ? (sline[2]::String,parse(Int64,sline[3])) : (I_3,r_3)
        length(sline)==3 && continue

        sline1 = sline[2]
        sline2 = sline[3]
        sline3 = sline[5]
        @views sline4 = reduce(*,replace.(sline[8:end],"I"=>"im"))
        # sline4 = reduce( * , replace.( sline4 , "I"=>"im" ) )
        # sline = [sline[2:3]...,sline[5],reduce(*,sline[8:end])]
        # sline[end] = reduce( * , replace.( sline[end] , "I"=>"im" ) )

        i1 = eval(parse(Int64,sline1))
        i2 = eval(parse(Int64,sline2))
        i3 = eval(parse(Int64,sline3))
        c  = ComplexF64(eval(Meta.parse(sline4)))
        # sline = map( x -> eval(Meta.parse(x)) , sline )
        if ! inverted
            push!( cg , (I_1,i1,I_2,i2,I_3,i3,r_3)=>c )
            # push!( cg , (I_1,sline[1]::Int64,I_2,sline[2]::Int64,I_3,sline[3]::Int64,r_3)=>sline[4] )
        else
            push!( cg , (I_1,i2,I_2,i1,I_3,i3,r_3)=>c )
            # push!( cg , (I_1,sline[2]::Int64,I_2,sline[1]::Int64,I_3,sline[3]::Int64,r_3)=>sline[4] )
        end
    end
    return cg
end

# orbital part multiplicity
function get_M_nonsimple( I , cg_path ; symmetry="D") 
    if isorbital(symmetry)

        cgo = cg_orbital_nonsimple( I , I , cg_path )
        return maximum([k[2] for k in keys(cgo)])

    else

        return 1

    end
end
# number of particles in irrep subspace
function get_N_nonsimple( 
            IJ::SF ,
            cg_path::String ,
            symmetry::String 
    ) where {SF<:Union{String,Float64}}

    if ispointspin(symmetry)

        return 2*get_M_nonsimple(IJ,cg_path;symmetry=symmetry)

    elseif isdoublegroup(symmetry)

        return get_M_nonsimple(IJ,cg_path;symmetry=symmetry)

    else # totalangularmomentum

        return jdim(IJ)

    end
end

function cg_shortcircuit_nonsimple( CG_PATH , oirreps... ; verbose=false )
    verbose && println( "recursion call" )
    seeds::Vector{String} = collect( oirreps )
    verbose && @show seeds
    produced::Vector{String} = []
    for seed_pair in with_replacement_combinations( seeds , 2 ) 
        cg_1 = cg_orbital_nonsimple( seed_pair[1] , seed_pair[2] , CG_PATH ; verbose=verbose )
        cg_2 = cg_orbital_nonsimple( seed_pair[2] , seed_pair[1] , CG_PATH ; verbose=verbose )
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
        return cg_shortcircuit_nonsimple( CG_PATH , produced... )
    end
end

function get_cg_o_fulldict_nonsimple( oirreps , cg_path )
    # Given a collection of orbital irreps 'oirreps', it searches in cg_path 
    # for CG information and returns the coefficients for every possible
    # combination (I_1,I_2) for I_1 and I_2 in oirreps:
    #
    #           cg[I_1,i_1,I_2,i_2,I_3,i_3] = ( I_1 , i_1 ; I_2 , i_2 | I_3 , i_3 )
    #
    cg_o_full = Dict{ Tuple{String,Int64,String,Int64,String,Int64,Int64} , ComplexF64 }()
    for I1 in oirreps, I2 in oirreps 
        merge!( cg_o_full , cg_orbital_nonsimple( I1 , I2 , cg_path ) )
    end
    return cg_o_full
end

