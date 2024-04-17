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
    for line in readlines( "$(path)/$(file)" ) 
        line=="" && continue
        sline = split(strip(line)," ")
        I_3,r_3 = length(sline)==3 ? (sline[2],parse(Int64,sline[3])) : (I_3,r_3)
        length(sline)==3 && continue
        sline = [sline[2:3]...,sline[5],reduce(*,sline[8:end])]
        sline[end] = reduce( * , replace.( sline[end] , "I"=>"im" ) )
        sline = map( x -> eval(Meta.parse(x)) , sline )
        if ! inverted
            push!( cg , (I_1,sline[1]::Int64,I_2,sline[2]::Int64,I_3,sline[3]::Int64,r_3)=>sline[4] )
        else
            push!( cg , (I_1,sline[2]::Int64,I_2,sline[1]::Int64,I_3,sline[3]::Int64,r_3)=>sline[4] )
        end
    end
    return cg
end

function get_M_nonsimple( I , cg_path ) 
    cgo = cg_orbital_nonsimple( I , I , cg_path )
    @show cgo
    return maximum([k[2] for k in keys(cgo)])
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

function compute_asymstates_N_doublegroup_nonsimple( 
            orbital::String , 
            N::Int64 ,
            cg_path::String ,
            asym_path::String ;
            verbose::Bool=false ,
            identityrep::String="" )

    # orbital dimensions
    M = get_M_nonsimple( orbital , cg_path )

    # basis
    basis_c = DoubleGroupCompoundCanonicalBasis( N , M )

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
            s_string = reduce( (x,y)->"- -" , s_part )

            write( io , "$o_string $s_string \n" )

        end
    end

    if (N==M && identityrep!=="")

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
    basis_s = SpinZeroCanonicalBasis(N)
    if verbose 
        println( "==========" )
        println( "SPIN BASIS" )
        println( "==========" )
        println( basis_s )
        println()
    end

    #%% BIPARTITIONS
    partitions_s = [Partition([N])]
    if verbose 
        println( "===============" )
        println( "SPIN PARTITIONS" )
        println( "===============" )
        for p in partitions_s
            println( "partition: $(p.vec)" )
            pp( p )
            println()
        end
    end

    #%% spin symstates
    dSsa2states_s = get_dSsa2states( basis_s , 
                                     partitions_s ;
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
    Iir2states_o = get_Iir2states_nonsimple( N , M , orbital , cg_path ; verbose=verbose )
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
    basis_c = DoubleGroupCompoundCanonicalBasis( N , M )
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
            phase_canrep = angle(filter(!iszero,canrep)[1])
            canstring = reduce( * , map(x->"($x)  ",canrep) ) 
            asymrep = collect(Base.Iterators.flatten(nullspace(hcat(asymsubspace,-canrep),atol=1e-6)[1:size(asymsubspace,2),:]))
            phase_asymrep = angle(filter(!iszero,asymrep)[1])
            asymrep *= exp(1.0im*(phase_canrep-phase_asymrep))
            #for e in asymrep 
            #    if abs(e)!==0.0 
            #        if real(e)<0.0
            #            asymrep = -asymrep
            #        end
            #        break 
            #    end
            #end
            normalize!(asymrep)
            asymstring = reduce( * , map(x->"($x)  ",asymrep) ) 
            toprint = reduce( (x,y)->x*"| "*y , [symstring,canstring,asymstring] )
            verbose && println( toprint )
            println( io , toprint )
        end
    end

    return ISisr
end
