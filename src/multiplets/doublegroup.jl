function compute_asymstates_N_doublegroup( 
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
