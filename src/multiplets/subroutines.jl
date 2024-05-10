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
function get_dSsa2states( basis_s::SpinZeroCanonicalBasis , 
                          bipartitions::Vector{Partition} ;
                          verbose=false )

    if verbose 
        println( "==========================================================" )
        println( "CALCULATION OF PERMUTATION- AND SPIN-SYMMETRIC SPIN STATES" )
        println( "==========================================================" )
    end

    dSsa2states = Dict{ Tuple{Tuple,Float64,Float64,Int64} , 
                        State }()
    space::Dict{ Tuple , Matrix{ComplexF64} } = Dict()

    for i in 1:length(basis_s)

        seed = basis_s[i]
        sortedspins = (seed.occupations...,)
        if sortedspins ∉ keys(space) 
            space[sortedspins] = Matrix{ComplexF64}(undef,length(basis_s),0) 
        end
        if verbose 
            println( "**********************" )
            println( "STATES FROM SEED $seed" )
            println( "......................" )
        end

        s = 0.0

        seedstate = State( seed , basis_s )

        verbose && println( "BIPARTITIONS" )
        for (bipart,tabs) in get_partitions2tableaux(bipartitions)

            S = 0.0
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
function get_dSsa2states( basis_s::JCanonicalBasis , 
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

        s = jz( seed )
        if verbose 
            println( "J_z = $s" )
            println()
        end
        seedstate = State( seed , basis_s )

        verbose && println( "PARTITIONS" )
        for (part,tabs) in get_partitions2tableaux(bipartitions)

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
function get_dar2states( basis_o::JCanonicalBasis , 
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

using PartialWaveFunctions
function get_Jjr2states( 
            N::Int64 ,
            J::Float64 ;
            verbose=false )

    basis_j = JCanonicalBasis(N,J)
    symstates_block::Dict{ Tuple{Float64,Float64,Int64} , Vector{Tuple{ComplexF64,Vector{Float64}}} } =
        Dict( (J,j,1)=>[(1.0im,[j])] for j=-J:J ) 
    symstates_add::Dict{ Tuple{Float64,Float64,Int64} , Vector{Tuple{ComplexF64,Vector{Float64}}} }=
        symstates_block

    for n in 2:N
        symstates_new = combine_symstates_J( symstates_add ,
                                             symstates_block ;
                                             verbose=verbose )
        symstates_block = symstates_new
        if verbose 
            for (k,vs) in symstates_block
                ss = reduce( + , [s[1]*State(JBasisElement(s[2]),JCanonicalBasis(n,J)) for s in vs] )
                verbose && @show ss
            end
        end
    end

    Jjr2states::Dict{ Tuple{Float64,Float64,Int64} , State } = Dict()
    for (k,vs) in symstates_block
        Jjr2states[k] = reduce( + , [s[1]*State(JBasisElement(s[2]),basis_j) for s in vs] )
    end

    return Jjr2states
end

function combine_symstates( symstates_block::Dict{ Tuple{String,Int64,Int64} , Vector{Tuple{ComplexF64,Vector{N}}} },
                            symstates_add::Dict{ Tuple{String,Int64,Int64} , Vector{Tuple{ComplexF64,Vector{N}}} },
                            cg_o::Dict{ Tuple{String,Int64,String,Int64,String,Int64} , ComplexF64 },
                            oirreps2dimensions;
                            verbose=false ) where {N<:Number}

    verbose && println( "COMBINING SYMSTATES..." )

    symnew = Dict{ Tuple{String,Int64,Int64} , Vector{Tuple{ComplexF64,Vector{N}}} }()

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
function combine_symstates_J( 
            symstates_block::Dict{ Tuple{Float64,Float64,Int64} , Vector{Tuple{ComplexF64,Vector{Float64}}} } ,
            symstates_add::Dict{ Tuple{Float64,Float64,Int64} , Vector{Tuple{ComplexF64,Vector{Float64}}} } ;
            verbose=false 
    )

    verbose && println( "COMBINING SYMSTATES..." )

    # combined symmetric states
    symnew = Dict{ Tuple{Float64,Float64,Int64} , Vector{Tuple{ComplexF64,Vector{Float64}}} }()

    # multiplets of the added and block sectors
    multiplets_a = Set((k[1],k[3]) for k in keys(symstates_add))
    multiplets_b = Set((k[1],k[3]) for k in keys(symstates_block))

    # iterate over multiplet combinations m_a,m_b
    for m_a in multiplets_a,
        m_b in multiplets_b

        if verbose 
            println("decomposition of multiplets m_a=$m_a and m_b=$m_b") 
        end

        # multiplet quantum numbers
        J_a, J_b = m_a[1], m_b[1]
        r_a, r_b = m_a[2], m_b[2]

        # iterate over decomposition of angular momentum product
        for J_c in abs(J_a-J_b):(J_a+J_b)

            # outer multiplicity
            r = get_r( symnew , J_c )
            verbose && @show r
            m_c = (J_c,r)

            verbose && println( "m_c=$m_c" )

            # iterate over irreps
            for j_c in -J_c:J_c,
                j_a in -J_a:J_a,
                j_b in -J_b:J_b

                # obtain clebsch-gordan coefficient
                cg_key = (J_a,j_a,J_b,j_b,J_c,j_c)
                coeff = ComplexF64(clebschgordan_doublearg(Int64.(2 .* cg_key)...))
                coeff==zero(ComplexF64) && continue

                verbose && println("$cg_key => $coeff") 

                vs_a = symstates_add[J_a,j_a,r_a]
                vs_b = symstates_block[J_b,j_b,r_b]
                vs_c = merge_vecstates(vs_a,vs_b,coeff)

                if (J_c,j_c,r) in keys(symnew)
                    append!( symnew[J_c,j_c,r] , vs_c )
                else
                    symnew[J_c,j_c,r] = vs_c
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

function fuse_states( vs::Vector{Tuple{ComplexF64,Vector{N}}} ) where {N<:Number}

    basisels = Set( x[2] for x in vs ) 
    fused::Vector{Tuple{ComplexF64,Vector{N}}} = []

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
@inline function get_r( symstates_new::Dict{ Tuple{Float64,Float64,Int64} , Vector{Tuple{ComplexF64,Vector{Float64}}} }, 
                        J::Float64 )::Int64

    if J in Set( k[1] for k in keys(symstates_new) )

        r::Int64 = maximum(Set(kk[3] for kk in keys(symstates_new) if kk[1]==J))
        r += 1
        return r

    end

    return 1
end

function merge_vecstates( vs_add::Vector{Tuple{ComplexF64,Vector{N}}}  , 
                          vs_block::Vector{Tuple{ComplexF64,Vector{N}}}  , 
                          coeff::ComplexF64=1.0 ) where {N<:Number}

    vs_new::Vector{Tuple{ComplexF64,Vector{N}}} = []

    for e_add in vs_add, 
        e_block in vs_block 

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


    if verbose 
        println()
        println("Calculating dIair2states.")
        println()
    end


    dIair2states::Dict{ Tuple{Tuple,String,Int64,Int64,Int64} , State } = Dict()

    # gather phases
    if verbose
        println("Storing phases")
        println("***")
    end
    phases = Dict{String,Vector{Tuple{Int64,Float64}}}()
    for I in Set(k[1] for k in keys(Iir2states))

        phases[I] = []

        dim_I = maximum([k[2] for k in keys(Iir2states) if k[1]==I])
        for i in 1:dim_I
            q = (I,i,1)
            state = Iir2states[q]
            idx = findfirst(!iszero,state.vec)
            c = angle(state.vec[idx])
            push!( phases[I] , (idx,c) )

            if verbose
                @show state
                @show c
                println("---")
            end
        end
    end

    if verbose
        println("Calculating intersections")
        println("***")
    end
    for (d,a) in Set(k[1:2] for k in keys(dar2states)),
        (I,i) in Set(k[1:2] for k in keys(Iir2states))

        S_dar = hcat([ s.vec for (k,s) in dar2states 
                             if k[1:2]==(d,a) ]...)
        S_Iir = hcat([ s.vec for (k,s) in Iir2states 
                             if k[1:2]==(I,i) ]...)

        int = subspace_intersection( S_dar , S_Iir )
        size(int,2)==0 && continue

        idx,phase = phases[I][i]
        if verbose
            @show d,a,I,i
            @show phase
        end
        for r in 1:size(int,2) 
            verbose && println( "pre-phase-fix:")
            state = State(int[:,r],basis_o)
            if verbose
                println(state)
                @show state.vec[idx],angle(state.vec[idx])
            end
            state.vec .*= cis(phase-angle(state.vec[idx]))
            if verbose
                println("post-phase-fix:")
                println(state)
                println("---")
            end
            dIair2states[d,I,a,i,r] = state
        end
    end

    return dIair2states
end
function get_dJajr2states( dar2states::Dict{ Tuple{Tuple,Int64,Int64} , State },
                           Jjr2states::Dict{ Tuple{Float64,Float64,Int64} , State } ,
                           basis_j::JCanonicalBasis ;
                           verbose=false )

    if verbose
        println()
        println("Calculating dJajr2states")
        println()
    end

    dJajr2states::Dict{ Tuple{Tuple,Float64,Int64,Float64,Int64} , State } = Dict()

    # gather phases
    #
    if verbose
        println("Storing phases")
        println("---")
    end
    phases = Dict{Float64,Vector{Tuple{Float64,Float64}}}()
    for J in Set(k[1] for k in keys(Jjr2states))

        phases[J] = []

        for j in -J:J
            q = (J,j,1)
            state = Jjr2states[q]
            idx = findfirst(!iszero,state.vec)
            c = angle(state.vec[idx])
            push!( phases[J] , (idx,c) )

            if verbose
                @show state
                @show c
                println("---")
            end
        end

    end

    if verbose
        println("Calculating intersection.")
        println("---")
    end
    for (d,a) in Set(k[1:2] for k in keys(dar2states)),
        (J,j) in Set(k[1:2] for k in keys(Jjr2states))

        S_dar = hcat([ s.vec for (k,s) in dar2states 
                             if k[1:2]==(d,a) ]...)
        S_Jjr = hcat([ s.vec for (k,s) in Jjr2states 
                             if k[1:2]==(J,j) ]...)

        int = subspace_intersection( S_dar , S_Jjr )
        size(int,2)==0 && continue

        idx = Int64(2J+1-(J-j))
        _,phase = phases[J][idx]

        if verbose
            @show d,a,J,j
            @show phase
        end

        for r in 1:size(int,2) 
            verbose && println( "pre-phase-fix:")
            state = State(int[:,r],basis_j)
            if verbose
                println(state)
                @show state.vec[idx],angle(state.vec[idx])
            end
            state.vec .*= cis(phase-angle(state.vec[idx]))
            if verbose
                println("post-phase-fix:")
                println(state)
            end
            dJajr2states[d,J,a,j,r] = state
        end
        verbose && println()
    end

    return dJajr2states
end

# ************
# dSsdIar part
# ............
function get_dSsadIair2states( dIair2states_o::Dict{ Tuple{Tuple,String,Int64,Int64,Int64} , State } ,
                               dSsa2states_s::Dict{ Tuple{Tuple,Float64,Float64,Int64} , State } ,
                               basis_c::CCBG ;
                               verbose=false ) where {CCBG<:CompoundCanonicalBasisGeneral}

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
# J version
function get_dJajdIair2states( dIair2states_o::Dict{ Tuple{Tuple,String, Int64,Int64,  Int64} , State } ,
                               dJajr2states_j::Dict{ Tuple{Tuple,Float64,Int64,Float64,Int64} , State } ,
                               basis_c::CCBG ;
                               verbose=false ) where {CCBG<:CompoundCanonicalBasisGeneral}

    verbose && println( "GETTING dJajdIair2states" )
    dJajdIair2states = Dict{ Tuple{Tuple,Float64,Int64,Float64,Tuple,String,Int64,Int64,Int64} , State }()

    multiplicities = Dict{ Tuple{Tuple{Tuple,Float64,Int64,Float64},Tuple{Tuple,String,Int64,Int64}} , Int64 }()

    for (k_o,s_o) in dIair2states_o,
        (k_j,s_j) in dJajr2states_j

        p_o = Partition([k_o[1]...])
        p_j = Partition([k_j[1]...])
        if verbose
            pp(p_o)
            pp(p_j)
        end
        p_o==complementary(p_j) || continue

        k_irreps_partners::Tuple{Tuple{Tuple,Float64,Int64,Float64},Tuple{Tuple,String,Int64,Int64}} = (k_j[1:end-1],k_o[1:end-1])
        if haskey( multiplicities , k_irreps_partners )
            multiplicities[k_irreps_partners] += 1
        else
            multiplicities[k_irreps_partners] = 1
        end
        r_c = multiplicities[k_irreps_partners]
        k_c = (k_j[1:end-1]...,k_o[1:end-1]...,r_c)
        s_c = tensormult( s_o , s_j , basis_c )
        dJajdIair2states[k_c] = s_c

        if verbose
            println( "INGREDIENTS AND RESULT" )
            println( "$k_o => $s_o" )
            println( "$k_j => $s_j" )
            @show s_c
            println()
        end

    end
    return dJajdIair2states
end


# ******************
# antisymmetric part 
# ..................
function get_M( basis_c::CompoundCanonicalBasis ) 
    return basis_c[length(basis_c)].orbital.occupations[1]*2
end
function get_M( basis_c::DoubleGroupCompoundCanonicalBasis ) 
    return basis_c[length(basis_c)].orbital.occupations[1]
end
function get_M( basis_c::JCompoundCanonicalBasis )
    return 2*basis_c[length(basis_c)].spin.occupations[1]+1
end
function get_N( basis_c::CCBG ) where {CCBG<:CompoundCanonicalBasisGeneral}
    return length(basis_c[1].orbital.occupations) 
end
function get_asymdim( basis_c::CCBG ) where {CCBG<:CompoundCanonicalBasisGeneral}
    M::Int64 = get_M( basis_c ) 
    N::Int64 = get_N( basis_c )
    return Int64(factorial(M)/factorial(M-N)/factorial(N))
end
function get_asymsubspace( N::Int64 ,
                           basis_c::CCBG ;
                           verbose=false ) where {CCBG<:CompoundCanonicalBasisGeneral}
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
    idx = 1
    for i in 1:length(basis_c)
        s = State( basis_c[i] , basis_c ) 
        #ts = ay*s 
        mul!( ts , ay , s )
        if verbose
            @show i
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
                         basis_c::CCBG ;
                         verbose=false ) where {CCBG<:CompoundCanonicalBasisGeneral}

    if verbose
        println( "===========================" )
        println( "ANTISYMMETRIC INTERSECTIONS" )
        println( "===========================" )
    end

    # get phases
    if verbose
        println()
        println("Storing phases:")
        println("***")
    end
    phases = Dict()
    for ((dS,S,s,aS,dI,I,aI,i,r),state) in dSsadIair2states

        r==1 || continue
        aS==1 || continue
        aI==1 || continue

        idx   = findfirst(x->!isapprox(x,zero(x);atol=1e-6),state.vec)
        phase = angle(state.vec[idx])


        phases[dS,S,s,dI,I,i] =  (idx,phase)
        if verbose
            @show state
            @show phase
            println("---")
        end
    end

    if verbose
        println("Calculating intersections")
        println("***")
    end
    ISisr = Dict{ Tuple{String,Float64,Int64,Float64,Int64} , State }()
    for sym in Set( (k[1:3]...,k[5],k[6],k[8]) 
                    for k in keys(dSsadIair2states) )

        if verbose
            @show sym
            println( "***********************"^2 )
        end

        (_,S,s,_,I,i) = sym

        idx,phase = phases[sym]
        verbose && (@show phase)

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

        if verbose
            println("Phase correction")
            println("***")
        end
        for (r,as) in enumerate(asymsym)
            as_phasecorrected = as*cis(phase-angle(as.vec[idx])) 
            ISisr[I,S,i,s,r] = as_phasecorrected
            if verbose
                @show as
                @show as_phasecorrected
                println("")
            end
        end
    end
    return ISisr
end
function get_asym_IJijr( dJajdIair2states::Dict{ Tuple{Tuple,Float64,Int64,Float64,Tuple,String,Int64,Int64,Int64} , State } ,
                         asymsubspace::Matrix{ComplexF64} ,
                         basis_c::CCBG ;
                         verbose=false ) where {CCBG<:CompoundCanonicalBasisGeneral}

    if verbose
        println( "===========================" )
        println( "ANTISYMMETRIC INTERSECTIONS" )
        println( "===========================" )
    end

    # get phases
    if verbose
        println()
        println("Storing phases")
        println("***")
    end
    phases = Dict()
    for ((dJ,J,aJ,j,dI,I,aI,i,r),state) in dJajdIair2states

        r==1  || continue
        aJ==1 || continue
        aI==1 || continue

        idx   = findfirst(x->!isapprox(x,zero(x);atol=1e-6),state.vec)
        phase = angle(state.vec[idx])

        phases[dJ,J,j,dI,I,i] =  (idx,phase)
        if verbose
            @show state
            @show phase
            println("---")
        end
    end

    IJijr = Dict{ Tuple{String,Float64,Int64,Float64,Int64} , State }()
    for sym in Set( (k[1],k[2],k[4],k[5],k[6],k[8]) for k in keys(dJajdIair2states) )

        if verbose
            @show sym
            println( "***********************"^2 )
        end

        (_,J,J,_,I,i) = sym

        idx,phase = phases[sym]
        verbose && (@show phase)

        sub = hcat([s.vec for (k,s) in dJajdIair2states 
                          if ((k[1],k[2],k[4])==sym[1:3] && (k[5],k[6],k[8])==sym[4:end])]...)

        int = subspace_intersection( asymsubspace , 
                                     sub ;
                                     verbose=false)

        asymsym = [State(int[:,i],basis_c) for i in 1:size(int,2)]
        size(int,2)==0 && continue

        if verbose
            for as in asymsym 
                @show as 
            end
            println()
        end

        verbose && println("\nPhase correction in antisymmetric intersections")
        for (r,as) in enumerate(asymsym)
            as_phasecorrected = as*cis(phase-angle(as.vec[idx])) 
            verbose && (@show as_phasecorrected)
            IJijr[I,J,i,j,r] = as_phasecorrected
        end
    end
    return IJijr
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


function combine_symstates_nonsimple( symstates_block::Dict{ Tuple{String,Int64,Int64} , Vector{Tuple{ComplexF64,Vector{Int64}}} },
                                      symstates_add::Dict{ Tuple{String,Int64,Int64} , Vector{Tuple{ComplexF64,Vector{Int64}}} },
                                      cg_o::Dict{ Tuple{String,Int64,String,Int64,String,Int64,Int64} , ComplexF64 },
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

        I_a,I_b = m_a[1], m_b[1]
        r_a,r_b = m_a[2], m_b[2]
        IIrr_c = Set((k[5],k[7]) for k in keys(cg_o) 
                            if (k[1],k[3])==(I_a,I_b))

        for (I_c,r_c) in IIrr_c

            r = get_r( symnew , I_c )
            verbose && @show r
            m_c = (I_c,r)

            verbose && println( "m_c=$m_c" )

            for i_c in 1:oirreps2dimensions[I_c],
                i_a in 1:oirreps2dimensions[I_a],
                i_b in 1:oirreps2dimensions[I_b]

                
                cg_key = (I_a,i_a,I_b,i_b,I_c,i_c,r_c) 
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

function get_Iir2states_nonsimple( 
            N::Int64 , 
            M::Int64 , 
            orbital::String , 
            cg_path::String ;
            verbose=false )

    basis_o = OrbitalCanonicalBasis(N,M)
    oirreps = cg_shortcircuit_nonsimple( cg_path , orbital )
    cg_o = get_cg_o_fulldict_nonsimple( oirreps , cg_path )
    oirreps2dimensions = get_oirreps2dimensions( cg_o ) 
    symstates_block::Dict{ Tuple{String,Int64,Int64} , Vector{Tuple{ComplexF64,Vector{Int64}}} }=
        Dict( (orbital,i,1)=>[(1.0,[i])] for i=1:oirreps2dimensions[orbital] ) 
    symstates_add::Dict{ Tuple{String,Int64,Int64} , Vector{Tuple{ComplexF64,Vector{Int64}}} }=
        symstates_block

    for n in 2:N
        symstates_new = combine_symstates_nonsimple( symstates_add ,
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

