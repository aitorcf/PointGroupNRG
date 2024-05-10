function compute_asymstates_N_totalangularmomentum( 
            J::Float64 , # total angular momentum J of the electrons
            N::Int64 , # number of electrons
            asym_path::String ;
            verbose::Bool=false )

    # orbital dimensions
    M = 2J+1

    # basis
    basis_c = JCompoundCanonicalBasis( N , J )

    # write down antisymmetric basis 
    asymsubspace = get_asymsubspace( N , basis_c , verbose=verbose )
    open( "$(asym_path)/N$(N)basis.txt" , "w" ) do io
        for i in 1:size(asymsubspace,2) 

            v = asymsubspace[:,i]
            idx = findfirst( x->(!isapprox(abs(x),0.0)) , v )
            basis_element = basis_c.states[idx]

            o_part = basis_element.orbital.occupations
            s_part = basis_element.spin.occupations

            o_string = reduce( (x,y)->"- -" , o_part )
            s_string = reduce( (x,y)->"$x $y" , s_part )

            write( io , "$o_string $s_string \n" )

        end
    end

    #if (N==M && identityrep!=="")

    #    if verbose 
    #        println( "************************************************" )
    #        println( "SHORTCUT COMPUTATION FOR THE SATURATED CASE N=2M" )
    #        println( "------------------------------------------------" )
    #    end

    #    s = State( asymsubspace[:,1] ) 
    #    sym = (identityrep,0.0,1,0.0,1)
    #    ISisr = Dict( sym => s )

    #    filename = "$(asym_path)/N$N.txt" 
    #    cmd = `touch $filename` 
    #    run(cmd)
    #    open( "$(asym_path)/N$N.txt" , "w" ) do io 
    #        for (k,s) in ISisr 
    #            symstring = reduce( * , map(x->"$x ",[k...]) )
    #            canrep = s.vec 
    #            canstring = reduce( * , map(x->"($x)  ",canrep) ) 
    #            asymrep = collect(flatten(nullspace(hcat(asymsubspace,-canrep),atol=1e-6)[1:size(asymsubspace,2),:]))
    #            normalize!(asymrep)
    #            asymstring = reduce( * , map(x->"($x)  ",asymrep) ) 
    #            toprint = reduce( (x,y)->x*"| "*y , [symstring,canstring,asymstring] )
    #            verbose && println( toprint )
    #            println( io , toprint )
    #        end
    #    end
    #    return ISisr
    #end


    # ******
    # J PART
    # ......

    if verbose
        println( "# · ###### · #" )
        println( "# | ------ | #" )
        println( "# | J PART | #" )
        println( "# | ------ | #" )
        println( "# · ###### · #" )
        println()
    end

    #%% basis
    basis_j = JCanonicalBasis(N,J)
    if verbose 
        println( "=======" )
        println( "J BASIS" )
        println( "=======" )
        println( basis_j )
        println()
    end

    #%% partitions
    partitions_j = [Partition([1 for _ in 1:N])]
    if verbose 
        println( "============" )
        println( "J PARTITIONS" )
        println( "============" )
        for p in partitions_j
            println( "partition: $(p.vec)" )
            pp( p )
            println()
        end
    end

    #%% permutation-symmetric j states
    dar2states_j = get_dar2states( basis_j ,
                                   partitions_j ;
                                   verbose=verbose )
    if verbose 
        println( "==============================" )
        println( "PERMUTATION-SYMMETRIC J STATES" )
        println( "==============================" )
        for (k,v) in dar2states_j
            println( "$k ==> $v" )
            println()
        end
    end

    #%% j-symmetric j states
    Jjr2states_j = get_Jjr2states( N , J ; verbose=verbose )
    if verbose 
        println( "====================" )
        println( "J-SYMMETRIC J STATES" )
        println( "====================" )
        for (k,v) in Jjr2states_j
            println( "$k ==> $v" )
            println()
        end
    end

    #%% permutation- and j-symmetric states
    dJajr2states_j = get_dJajr2states( dar2states_j , 
                                       Jjr2states_j ,
                                       basis_j ;
                                       verbose=verbose )
    if verbose
        println( "===================================" )
        println( "J- AND PERMUTATION-SYMMETRIC STATES" )
        println( "===================================" )
        clean_symstates!( dJajr2states_j )
        for (k,v) in dJajr2states_j
            println( "$k ==> $v" ) 
        end
        println()
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
    basis_o = OrbitalNullCanonicalBasis( N )
    if verbose 
        println( "=============" )
        println( "ORBITAL BASIS" )
        println( "=============" )
        show( basis_o )
        println()
    end

    #%% partitions
    partitions_so = Dict( p=>complementary(p) for p in partitions_j )
    partitions_o = [complementary(p) for p in partitions_j]
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

    # permutation- and orbital-symmetric orbital states
    dIair2states_o = Dict{ Tuple{Tuple,String,Int64,Int64,Int64} , State }(
        ((N,),"A",1,1,1) => State(basis_o.states[1],basis_o)
    )
    if verbose
        println( "=================================================" )
        println( "ORBITAL- AND PERMUTATION-SYMMETRIC ORBITAL STATES" )
        println( "=================================================" )
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
    basis_c = JCompoundCanonicalBasis( N , J )
    if verbose 
        println( "==============" )
        println( "COMPOUND BASIS" )
        println( "==============" )
        show(basis_c)
    end

    #%% independently symmetric states
    dJajdIair2states = get_dJajdIair2states( dIair2states_o , 
                                             dJajr2states_j ,
                                             basis_c ;
                                             verbose=verbose )
    if verbose
        println( "================" )
        println( "dJajdIair STATES" )
        println( "================" )
        for (k,v) in dJajdIair2states 
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
    IJijr = Dict{ Tuple{String,Float64,Int64,Float64,Int64} , State }(
        (I,J,i,j,r)=>state
        for ((_,J,_,j,_,I,_,i,r),state) in dJajdIair2states
    )
    if verbose
        println( "============" )
        println( "IJijr STATES" )
        println( "============" )
        for (k,v) in IJijr
            @show k 
            @show v 
            println()
        end
    end

    #ISisr = get_asym_ISisr( dSsadIair2states , 
    #                        asymsubspace ,
    #                        basis_c ;
    #                        verbose=verbose )
    #clean_symstates!( ISisr )

    # *******
    # WRITING
    # .......
    filename = "$(asym_path)/N$N.txt" 
    cmd = `touch $filename` 
    run(cmd)
    open( "$(asym_path)/N$N.txt" , "w" ) do io 
        for (k,s) in IJijr 
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

    return IJijr
end
