using Printf
function compute_asymstates_N_allsymmetries( 
            symmetry::String ,
            orbital::SF , 
            N::Int64 ,
            asym_path::String ;
            cg_path::String="" ,
            verbose::Bool=false ,
    ) where {SF<:Union{String,Float64}}

    symmetryname = begin
        if ispointspin(symmetry)
            "Orbital point group ⊗ spin SU(2)"
        elseif isdoublegroup(symmetry)
            "Double group"
        else
            "Spin-orbital SU(2)"
        end
    end
    println("| ==========================================")
    println("| Calculation of multiplet states.")
    println("|  - symmetry: $symmetryname.")
    println("|  - orbital: $orbital.")
    println("| ==========================================")
    println()
    println("* Notation:")
    println("*  P   : partition, permutation-symmetry class.")
    println("*  p   : member of permutation-symmetry class.")
    println("*  I   : orbital irrep.")
    println("*  i   : orbital partner.")
    println("*  S/J : total spin (S) / angular momentum (J).")
    println("*  s/j : projection of spin (s) / total angular momentum (j).")
    println("*  r   : outer multiplicity.")
    println()

    if istotalangularmomentum(symmetry)
        J::Float64 = orbital
    end

    # orbital dimensions
    M::Int64 = begin
        if isorbital(symmetry)
            get_M_nonsimple( orbital , cg_path )
        else
            2J+1
        end
    end

    # basis
    basis_c = begin
        if isdoublegroup(symmetry)

            DoubleGroupCompoundCanonicalBasis( N , M )

        elseif ispointspin(symmetry)

            CompoundCanonicalBasis( N , M )

        else # total angular momentum

            JCompoundCanonicalBasis( N , J )

        end
    end

    # write down antisymmetric basis 
    asymsubspace = get_asymsubspace( N , basis_c , verbose=verbose )
    open( "$(asym_path)/N$(N)basis.txt" , "w" ) do io
        for i in 1:size(asymsubspace,2) 

            v = asymsubspace[:,i]
            idx = findfirst( x->(!isapprox(abs(x),0.0)) , v )
            basis_element = basis_c.states[idx]

            o_part = basis_element.orbital.occupations
            s_part = basis_element.spin.occupations

            o_string = begin
                if isorbital(symmetry)
                    reduce( (x,y)->"$x $y" , o_part )
                else
                    reduce( (x,y)->"- -" , o_part )
                end
            end
            s_string = begin
                if !isdoublegroup(symmetry)
                    reduce( (x,y)->"$x $y" , s_part )
                else
                    reduce( (x,y)->"- -" , s_part )
                end
            end

            write( io , "$o_string $s_string \n" )

        end
    end

    # *********
    # SPIN PART
    # .........

    if verbose
        println( "# · ######## · #" )
        println( "# | -------- | #" )
        println( "# | S/J PART | #" )
        println( "# | -------- | #" )
        println( "# · ######## · #" )
        println()
    end

    #%% BASIS
    basis_s = begin
        if ispointspin(symmetry)

            SpinCanonicalBasis(N)

        elseif isdoublegroup(symmetry)

            SpinZeroCanonicalBasis(N)

        else # total angular momentum

            JCanonicalBasis(N,J)

        end
    end
    println( "=========" )
    println( "S/J BASIS" )
    println( "=========" )
    println( basis_s )
    println()

    #%% spin partitions
    partitions_s = begin
        if ispointspin(symmetry)

            bipartitions(N)

        elseif isdoublegroup(symmetry)

            [Partition([N])]

        else

            [Partition([1 for _ in 1:N])]

        end
    end
    println( "==============" )
    println( "S/J PARTITIONS" )
    println( "==============" )
    for p in partitions_s
        println( "partition: $(p.vec)" )
        pp( p )
        println()
    end

    #%% S/J symstates
    # TODO: Check whether dSsa2states needs r label
    states_s = begin
        if isorbital(symmetry)

            get_dSsa2states( 
                basis_s , 
                partitions_s ;
                verbose=verbose 
            )

        else

            Jjr2states = get_Jjr2states( N , J ; verbose=verbose )

            get_dJajr2states_alternative(
                partitions_s,
                Jjr2states,
                basis_s
            )

        end
    end
    println( "==========================================" )
    println( "S/J- AND PERMUTATION-SYMMETRIC STATES" )
    println( "==========================================" )
    println()
    @printf  "%-14s %-3s %-3s %-3s\n" "P" "p" "S/J" "s/j"
    @printf        "state\n"
    if istotalangularmomentum(symmetry)
        for ((d,J,a,j,r),state) in sort(collect(states_s),by=x->x[1])
            @printf "%-14s %-3i %-3.1f %-3.1f\n" d a J j
            println(state)
            println()
        end
    else
        for ((d,S,s,a),state) in sort(collect(states_s),by=x->x[1])
            @printf "%-14s %-3i %-3.1f %-3.1f\n" d a S s
            println(state)
            println()
        end
    end


    # ************
    # ORBITAL PART 
    # ............
    println( "# · ############ · #" )
    println( "# | ------------ | #" )
    println( "# | ORBITAL PART | #" )
    println( "# | ------------ | #" )
    println( "# · ############ · #" )
    println()

    #%% orbital basis
    basis_o = begin
        if isorbital(symmetry)

            OrbitalCanonicalBasis( N , M )

        else

            OrbitalNullCanonicalBasis(N)

        end
    end
    println( "=============" )
    println( "ORBITAL BASIS" )
    println( "=============" )
    println( basis_o )
    println()

    #%% partitions
    partitions_so = Dict( p=>complementary(p) for p in partitions_s )
    partitions_o = [complementary(p) for p in partitions_s]
    println( "==========================================" )
    println( "ORBITAL PARTITIONS (complementary to S/J)" )
    println( "==========================================" )
    println()
    @printf "%-20s <== %-20s\n" "orbital partition" "complementary spin partition"
    println()
    for (sp,op) in partitions_so 
        @printf "%-20s <== %-20s\n" op sp
        pp( op )
        println()
    end

    #%% orbital- and permutation-symmetric orbital states 
    states_o = begin
        if isorbital(symmetry)

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

            get_dIair2states_alternative( 
                partitions_o , 
                Iir2states_o ,
                basis_o ;
                verbose=verbose 
            )

        else

            Dict{ Tuple{Tuple,String,Int64,Int64,Int64} , State }(
                ((N,),"A",1,1,1) => State(basis_o.states[1],basis_o)
            )

        end
    end
    println( "=================================================" )
    println( "ORBITAL- AND PERMUTATION-SYMMETRIC ORBITAL STATES" )
    println( "=================================================" )
    clean_symstates!( states_o )
    @printf "%-14s %-3s %-5s %-3s %-3s\n" "P" "p" "I" "i" "r"
    @printf "state\n"
    println()
    for ((d,I,a,i,r),state) in states_o
        @printf "%-14s %-3s %-5s %-3i %-3i\n" d a I i r
        println(state)
        println()
    end
    println()

    # *************
    # COMBINED PART
    # .............
    println( "# · ############# · #" )
    println( "# | ------------- | #" )
    println( "# | COMBINED PART | #" )
    println( "# | -------.----- | #" )
    println( "# · ############# · #" )
    println()

    #%% compound basis 
    println( "==============" )
    println( "COMPOUND BASIS" )
    println( "==============" )
    show(basis_c)
    println()

    #%% independently symmetric states
    states_so = begin
        if isorbital(symmetry)

            get_dSsadIair2states( 
                states_o , 
                states_s ,
                basis_c ;
                verbose=verbose 
            )

        else

            get_dJajdIair2states(
                states_o , 
                states_s ,
                basis_c ;
                verbose=verbose 
            )

        end
    end
    println( "============================================" )
    println( " STATES WITH COUPLED ORBITAL AND S/J SECTORS" )
    println( "============================================" )
    println()
    @printf "%-14s %-4s %-3s %-3s " "P_sj" "p_sj" "S/J" "s/j"
    @printf "%-14s %-4s %-3s %-3s"  "P_o"  "p_o"  "I"   "i"
    @printf "%-3s\n" "r"
    @printf "state\n"
    for ((d_j,J,a_j,j,d_o,I,a_o,i,r),state) in states_so 
        @printf "%-14s %-4i %-3.1f %-3.1f " d_j a_j J j
        @printf "%-14s %-4i %-3s %-3i "     d_o a_o I i
        @printf "%-3s\n" r
        println(state)
        println()
    end

    # *************************************
    # SYMMETRY-ADAPTED ANTISYMMETRIC STATES
    # .....................................
    #%%
    ISisr = begin 
        if isorbital(symmetry)

            get_asym_ISisr_alternative( 
                states_so , 
                basis_c ;
                verbose=verbose 
            )

        else

            Dict{ Tuple{String,Float64,Int64,Float64,Int64} , State }(
                (I,J,i,j,r)=>state
                for ((_,J,_,j,_,I,_,i,r),state) in states_so
            )

        end
    end
    clean_symstates!( ISisr )
    println("====================================")
    println("ANTISYMMETRIC STATES (phase unfixed)")
    println("====================================")
    println()
    @printf "%-4s %-4s %-4s %-4s %-4s\n" "I" "i" "S/J" "s/j" "r"
    @printf "state\n"
    println()
    for ((I,S,i,s,r),state) in ISisr
        @printf "%-4s %-4.1f %-4i %-4.1f %-4i\n" I S i s r
        println(state)
        println()
    end
    println("Read final input in $(asym_path)/")

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
            normalize!(asymrep)
            asymstring = reduce( * , map(x->"($x)  ",asymrep) ) 
            toprint = reduce( (x,y)->x*"| "*y , [symstring,canstring,asymstring] )
            verbose && println( toprint )
            println( io , toprint )
        end
    end

    return ISisr
end
