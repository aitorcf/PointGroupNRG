const joirrep = "A"

# dimension of subspace with total angular momentum J
jdim( J::Float64 ) = 2J+1

# convert J to int by doubling
doublej( J::Float64 )::Int64 = Int64(2J)

# j1 ⊠ j2 = … ⊕ j3 ⊕ … ?
function isj3inj1timesj2(j1::Float64,j2::Float64,j3::Float64)
    (abs(j1-j2) ≤ j3 ≤ (j1+j2)) || return false
    isinteger(j1+j2) && (isinteger(j3) || return false)
    !isinteger(j1+j2) && (!isinteger(j3) || return false)
    return true
end
# integer version (j means 2×j)
function isj3inj1timesj2(j1::Int64,j2::Int64,j3::Int64)
    (abs(j1-j2) ≤ j3 ≤ (j1+j2)) || return false
    iseven(j1+j2) && (iseven(j3) || return false)
    isodd(j1+j2) && (isodd(j3) || return false)
    return true
end
 
# ------------------- #
# BASIS AND SYMSTATES #
# ------------------- #

# atomic states with the modified pointspin convention
#   | r I i m ⟩ → | r J j - ⟩
function get_atomic_states_totalangularmomentum( 
            external_label::Int64 ,
            J::Float64 )::Vector{Tuple{Int64,String,Int64,String}}

    root = ( external_label , string(doublej(J)) )
    statetuples = Tuple{Int64,String,Int64,String}[]

    for j in -J:J

        push!( statetuples , (root...,doublej(j),"-") )

    end

    return statetuples
end

function symstates_n_oneirrep_totalangularmomentum(
            basis::CanonicalBasis{SFS} , 
            n::Int64 , 
            symstates::Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , S } ;
            filename::String="" ,
            identity::String=joirrep )::Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , S } where {SFS<:SymbolFockState,S<:State}

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
            (g::Int64,J2string::String,j2::Int64,_) = sfs.hilbert.states[ findfirst(sfs.occ) ]
            get!( symstates ,
                  (1,joirrep,parse(Float64,J2string)/2.0,1,Float64(j2/2.0),g) , 
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

function oneirrep_symstates_totalangularmomentum( 
            hilbert::HS , 
            irrepath::String ) where {HS<:HilbertSpace}

    basis = CanonicalBasis(hilbert)
    basis_number = [CanonicalBasis(hilbert,i) for i=0:length(hilbert.states)]

    symstates::Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , State } = Dict()
    for (i::Int64,b) in enumerate(basis_number)
        n = i-1
        symstates = symstates_n_oneirrep_totalangularmomentum( 
            b ,
            n ,
            symstates ;
            filename=irrepath*"N$n.txt" , 
        )
    end
    for (k,v) in symstates 
        symstates[k] = extend_basis( v , basis )
    end

    return symstates
end

function cg_reduce_product_states_totalangularmomentum( 
            symstates_1::Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , S } , 
            symstates_2::Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , S } 
    )::Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , S } where {S<:State}

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

    I_1 = I_2 = joirrep

    symstates = Dict()
    for ircomb in irrep_combinations

        (N_1::Int64,_,J_1::Float64,r_1::Int64) = ircomb[1]
        (N_2::Int64,_,J_2::Float64,r_2::Int64) = ircomb[2]

        N_3::Int64 = N_1 + N_2
        I_3::String = "A"

        for J_3 in abs(J_1-J_2):(J_1+J_2)

            ro_3::Int64 = 1

            for j_1 in -J_1:J_1,
                j_2 in -J_2:J_2,
                j_3 in -J_3:J_3

                r_3 = (r_1,r_2,ro_3)

                c_j = clebschgordan(J_1,j_1,J_2,j_2,J_3,j_3)

                q_3::Tuple{ NTuple{3,Int64} , NTuple{3,String} , NTuple{3,Float64} , Int64 , Float64 , NTuple{3,Int64} } = 
                    ((N_1,N_2,N_3),(I_1,I_2,I_3),(J_1,J_2,J_3),1,j_3,r_3)

                s1::S = symstates_1[( N_1 , I_1 , J_1 , 1 , j_1 , r_1 )]
                s2::S = symstates_2[( N_2 , I_2 , J_2 , 1 , j_2 , r_2 )]
                merge!( + , symstates , Dict(q_3 => c_j*x(s1,s2)) )

            end
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
        J_3::Float64 = qnums[3][3]
        mu_3::Int64  = qnums[4]
        m_3::Float64 = qnums[5]

        p::Tuple{Int64,String,Float64,Int64,Float64} = ( N_3 , I_3 , J_3 , mu_3 , m_3 )
        merge!( + , principals , Dict(p => 1) )

        q_3::Tuple{Int64,String,Float64,Int64,Float64,Int64} = ( p... , principals[p] )

        get!( symstates_new , q_3 , state )
    end

    return symstates_new
end

function get_symstates_basis_multiplets_totalangularmomentum( 
            atom_config::Dict{Float64,Int64},
            asym_dir::String ;
            verbose::Bool=true )

    # symstates
    separate_symstates = []
    for (J,multiplicity) in atom_config
        for m in 1:multiplicity
            onemult_states = get_atomic_states_totalangularmomentum(m,J)
            onemult_hilbert = HilbertSpace( onemult_states )
            onemult_symstates = oneirrep_symstates_totalangularmomentum( 
                onemult_hilbert ,
                "$(asym_dir)/$(J)/" 
            )
            push!( separate_symstates , onemult_symstates )
        end
    end
    symstates = separate_symstates[1]
    for symstates_new in separate_symstates[2:end]
        symstates = cg_reduce_product_states_totalangularmomentum(
                        symstates , 
                        symstates_new )
    end

    # basis 
    basis::CanonicalBasis = collect(values(symstates))[1].basis 

    # multiplets 
    multiplets::Set{Tuple{Int64,String,Float64,Int64}} = get_multiplets( symstates )
    multiplets_a::Set{Tuple{Int64,String,Float64,Int64}} = Set( m for m in multiplets if m[1]==1 )

    # printing
    if verbose 
        println( "BASIS" )
        println( basis )
        println()
        println( "SYMSTATES" )
        print_symstates_ordered( symstates )
        println()
        println( "MULTIPLETS" )
        print_multiplets_Nordered( multiplets )
        println()
        println( "ONE-ELECTRON MULTIPLETS" )
        print_multiplets_Nordered( multiplets_a )
        println()
    end

    return (symstates,basis,multiplets,multiplets_a)
end

# ------------------ #
# LEHMANN AMPLITUDES #
# ------------------ #

function shell_coperators_totalangularmomentum( 
            shell_basis::CB 
    )::Dict{ ClearQNums , Operator{CB} } where {CB<:CanonicalBasis,D<:Dict}

    # same as basis2coperators (above) but stores the creation
    # operators as a dictionary.
    #
    # input
    #
    # - shell_basis : canonical basis of shell.
    # - hiztegia : irrep pseudoname => irrep name
    #
    # output
    #
    # - coperator : sym => creation operatorcoperator : sym => creation operator
    #

    coperators = Dict{ ClearQNums , Operator{CB} }()
    @inbounds for (i,sfs) in enumerate(shell_basis.states)

        sum(sfs.occ)!==1 && continue

        idx = findfirst( sfs.occ )
        tup = sfs.hilbert.states[idx] # Int, String, Int, Int

        # sym = ( N , I , J , i , j , r )
        sym::ClearQNums = (
            1 , 
            joirrep , 
            parse(Float64,tup[2])/2.0 ,
            1 ,
            tup[3]/2.0 , 
            tup[1]
        )
        cop::Operator{CB} = sfs2coperator(sfs,shell_basis)
        sym2cop::Dict{ ClearQNums , Operator{CB} } = Dict{ClearQNums,Operator{CB}}( sym => cop )
        merge!( coperators , sym2cop )

    end

    return coperators 
end

# ClearX to IntX conversion
function clear2int_totalangularmomentum( multiplet::ClearMultiplet )::IntMultiplet
    ( N , _ , J , r ) = multiplet
    return ( N , 1 , doublej(J) , r )
end
function clear2int_totalangularmomentum( irrep::ClearIrrep )::IntIrrep
    ( N , _ , J ) = irrep 
    return ( N , 1 , doublej(J) )
end
function clear2int_totalangularmomentum( q::ClearQNums )::IntQNums
    ( N , _ , J , _ , j , r ) = q
    return (N,1,doublej(J),1,doublej(j),r)
end

function get_lehmann_totalangularmomentum( 
                symstates_shell::Dict{ClearQNums,S} , 
                basis_shell::CB ;
                verbose=false )::IntQPCG where {CB<:CanonicalBasis,S<:State,D<:Dict}
    # compute the pseudo-CG coefficients for the shells.
    #
    # input
    #
    # - symstates_shell : all symstates of a single shell (any)
    # - basis_shell : basis for only that shell 
    # - hiztegia : dict( irrep name => standard notation )
    # - [format : "int" or "standard"]
    #
    # output
    #
    # - pseudo-CG coefficients as dict( (q_nu,q_a,q_mu)=>coeff )
    
    pseudoCG_clear::ClearQPCG = ClearQPCG()

    # shell creation operators
    shell_cops::Dict{ ClearQNums , Operator{CB} } = shell_coperators_totalangularmomentum( basis_shell )

    for (q_nu::ClearQNums,s_nu::S) in symstates_shell, 
        (q_mu::ClearQNums,s_mu::S) in symstates_shell

        for (q_a::ClearQNums,c_a::Operator{CB}) in shell_cops 
            q::ClearTripleQ = ( q_nu , q_a , q_mu ) 
            # ⟨ q_nu | c†_{q_a} | q_mu ⟩
            # notice that it is not the complex conjugate!!!
            coeff::ComplexF64 = (s_nu * c_a * s_mu)::ComplexF64
            isapprox( abs(coeff) , 0.0 , atol=1e-6 ) || push!( pseudoCG_clear , q=>coeff ) 
        end

    end

    pseudoCG_int::IntQPCG = IntQPCG( (clear2int_totalangularmomentum(k[1]),
                                      clear2int_totalangularmomentum(k[2]),
                                      clear2int_totalangularmomentum(k[3]))=>v 
                                      for (k::ClearTripleQ,v::ComplexF64) in pseudoCG_clear )

    if verbose 
        for (q,coeff) in pseudoCG_int
            println( "q = $q" )
            println( "coeff = $coeff" )
            println()
        end
    end

    return pseudoCG_int::IntQPCG
end

function get_redmat_nonsimple_totalangularmomentum( 
            dictmat , 
            multiplets_atom , 
            multiplets_operator ;
            verbose=false )

    if verbose 
        println( "~~~~~~~~~~~~~~~~~~~~~~~~~~" )
        println( "REDUCED MATRIX COMPUTATION" )
        println( "~~~~~~~~~~~~~~~~~~~~~~~~~~" )
        println()
        @show multiplets_atom 
        @show multiplets_operator
        println()
    end

    # irrep combinations in dictmat
    Gcombs = Set( (q_i[1:3],q_a[1:3],q_j[1:3]) 
                  for (q_i,q_a,q_j) in keys(dictmat) )
    irrmult_ij = get_irreps( Set(multiplets_atom) ; multiplicity=true )
    irrmult_a  = get_irreps( Set(multiplets_operator) ; multiplicity=true )

    if verbose
        @show irrmult_ij 
        @show irrmult_a
        @show Gcombs
        println()
    end

    # reduced matrix with zeros
    #       < m_i || f^+_{m_a} || m_j > = 0
    irr2mult_ij = Dict( G=>R for (G,R) in irrmult_ij )
    irr2mult_a  = Dict( G=>R for (G,R) in irrmult_a  )
    redmat::Dict{NTuple{3,NTuple{3,Int64}},Array{ComplexF64,4}} = 
        Dict( (G1,G2,G3)=>zeros(ComplexF64,1,irr2mult_ij[G1],irr2mult_a[G2],irr2mult_ij[G3]) 
              for (G1,G2,G3) in Gcombs )

    # iterate over irreps (with multiplicities)
    @inbounds for (G_i,R_i) in irrmult_ij,
                  (G_a,R_a) in irrmult_a,
                  (G_j,R_j) in irrmult_ij

        ((G_i,G_a,G_j) in Gcombs) || continue

        # irrep quantum numbers
        (N_i,_,J2_i) = G_i
        (N_j,_,J2_j) = G_j
        (N_a,_,J2_a) = G_a

        # filter irrep combination
        N_i==(N_j+N_a)                  || continue 
        isj3inj1timesj2(J2_a,J2_j,J2_i) || continue

        if verbose
            println( "IRREPS" )
            @show G_i,G_a,G_j
            println()
        end

        # iterate over irrep multiplets 
        @inbounds for r_i in 1:R_i,
                      r_a in 1:R_a,
                      r_j in 1:R_j

            verbose && (@show r_i,r_a,r_j)

            # SU(2) is simply reducible
            α = 1

            # Store  < m_i || f^+_{m_a} || m_j > = 0
            # without having to search in array
            matel = zero(ComplexF64)

            # arbitrary
            j2_i = J2_i 
            for j2_j in -J2_j:2:J2_j,
                j2_a in -J2_a:2:J2_a

                # full quantum numbers
                q_i = (G_i...,1,j2_i,r_i)
                q_a = (G_a...,1,j2_a,r_a)
                q_j = (G_j...,1,j2_j,r_j)

                # lehmann amplitude
                lehmann = get(dictmat,(q_i,q_a,q_j),zero(ComplexF64))
                iszero(lehmann) && continue

                # clebsch-gordan
                cg = clebschgordan_doublearg(J2_a,j2_a,J2_j,j2_j,J2_i,j2_i)

                # compute reduced matrix element
                matel += cg*lehmann

                if verbose
                    @show q_i,q_a,q_j
                    @show lehmann
                    @show cg
                    @show matel
                    println()
                end

            end

            redmat[(G_i,G_a,G_j)][α,r_i,r_a,r_j] = matel

        end
    end
    return redmat
end

function pcg_nonsimple_sanity_check_totalangularmomentum( 
            pcgdict::IntQPCG ,
            pcgred::IntIrrepPCGNS ;
            verbose::Bool=false )

    if verbose
        println()
        println("NON REDUCED")
        print_dict(pcgdict)
        println()
        println("REDUCED")
        print_dict(pcgred)
        println()
    end

    for ((q_u,q_a,q_v),matel) in pcgdict

        Nu,Iu,J2u,iu,j2u,ru = q_u
        Na,Ia,J2a,ia,j2a,ra = q_a
        Nv,Iv,J2v,iv,j2v,rv = q_v

        Gu = (Nu,Iu,J2u)
        Ga = (Na,Ia,J2a)
        Gv = (Nv,Iv,J2v)

        cg = clebschgordan_doublearg(J2a,j2a,J2v,j2v,J2u,j2u)

        if verbose
            println("Dict matrix element:")
            @show q_u,q_a,q_v
            @show matel
            println()
            println( "------------" )
            println("Orbital Clebsch-Gordan coefficients")
            @show Gu,Ga,Gv
            @show Ia,Iu,Iv

            α = 1
            println( "outer multipliticty α=$α" )

            for j2a in -J2a:2:J2a,
                j2u in -J2u:2:J2u,
                j2v in -J2v:2:J2v

                c = clebschgordan_doublearg(J2a,j2a,J2v,j2v,J2u,j2u)
                println( "( $j2a $j2v | $j2u ) = $c" )

            end
            println()
            println( "------------" )
        end

        wignereckart = 0.0
        α = 1

        red = pcgred[Gu,Ga,Gv][α,ru,ra,rv]
        verbose && println( "red = ",red," || ",j2a,j2v,j2u," → " , cg )
        wignereckart += red*conj(cg)

        if verbose
            println("Comparison:")
            println("  - dictionary element:       $matel")
            println("  - reduced × clebsch-gordan: $wignereckart")
            println()
        end
        if !isapprox(matel,wignereckart;atol=1e-6)
            error("Error in reduced matrix element.")
        end

    end
    verbose && println("Reduced Lehmann amplitudes correct!")
end

function get_lehmann_reduced_totalangularmomentum(
            basis::CB ,
            symstates_noint::ClearSymstateDict ,
            multiplets::IntMultipletSet ;
            verbose::Bool=false )::IntIrrepPCGNS where {CB<:CanonicalBasis,D<:Dict}

    # pcg in dict format
    pcg::IntQPCG = get_lehmann_totalangularmomentum(
        symstates_noint ,
        basis 
    )

    # pcg in reduced matrix format
    multiplets_a::IntMultipletVector = collect(filter( x->x[1]==1 , multiplets ))
    pcgred = get_redmat_nonsimple_totalangularmomentum(
        pcg,
        multiplets ,
        multiplets_a ;
        verbose=verbose 
    )
    if verbose
        println( "PCGRED" )
        print_dict( pcgred )
        println()
    end

    pcg_nonsimple_sanity_check_totalangularmomentum( pcg , pcgred ; verbose=verbose )

    return pcgred
end

# ================= #
# HAMILTONIAN TERMS #
# ================= #

# occupation energies
function epsilon_sym_totalangularmomentum( symstates , symparams::Dict{Float64,Vector{ComplexF64}} ; verbose=false )

    # one-particle symstates
    symstates_n1 = symstates_n( symstates , 1 )

    # basis
    basis = pickrandom(symstates_n1)[2].basis
    n = length( basis.states )

    # creation operators
    coperators = basis2coperators(basis)

    # initialize operator
    epsop = Operator(0,basis)

    # iterate through one-particle symstates
    for (q,s) in symstates_n1

        verbose && println(s)

        # construct the corresponding part of the operator
        for (i,component) in enumerate(s.vector)
            component==0 && continue
            cop = coperators[i]
            epsop += symparams[q[2]][q[end]]*component*cop*adjoint(cop)
        end
    end

    return epsop
end
# coulomb repulsion
function u_sym_totalangularmomentum( symstates , symparams ; verbose=false )

    if verbose 
        println( "CONSTRUCTING U..." )
        println()
    end

    # two-particle symstates
    symstates_n2 = symstates_n( symstates , 2 )

    # basis
    basis = pickrandom( symstates_n2 )[2].basis 

    # creation operators
    coperators = basis2coperators(basis)

    # initialize operator
    u = Operator( 0 , basis )

    # u = ∑_{q1q2} c†_q1 c_q2
    for (q1,s1) in symstates_n2, 
        (q2,s2) in symstates_n2 

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

        # iterate through slater expansions
        for (i1,c1) in enumerate(s1.vector), 
            (i2,c2) in enumerate(s2.vector)

            (c1==0 || c2==0) && continue 

            if verbose 
                println( "coefficients: c1=$c1, c2=$c2" )
                println( "states" )
                println( basis.states[i1] )
                println( basis.states[i2] )
            end

            creop = coperators[i1]
            annop = adjoint(coperators[i2])
            c = symparams[q1[1:2]][q1[end],q2[end]]*c1*conj(c2)
            u += c * creop*annop
            verbose && println( "( sfs1 | u | sfs2 ) = $c" )
        end
        verbose && println()
    end

    return u
end
# number operator
function electron_counter_sym_totalangularmomentum(
            symstates::Dict{ ClearQNums , State } , # noint
            multiplet::Tuple{Float64,Int64} ;
            verbose=false
    )

    # gather one-particle symstates
    symstates_n1 = symstates_n(symstates,1)

    # basis
    basis = pickrandom(symstates_n1)[2].basis

    # creation operators for all the one-electron states
    coperators = basis2coperators( basis )

    # initialize counting operator
    counter = Operator( 0 , basis )

    # select symstate from selected multiplet
    symstates_chosen = [ s for (q,s) in symstates if (q[1],q[3],q[end])==(1,multiplet...) ]

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

# =========================== #
# HAMILTONIAN DIAGONALIZATION #
# =========================== #

# multiplet combinations
function get_combinations_Gu_muiualpha_totalangularmomentum( 
            multiplets_block::Set{NTuple{4,Int64}} , 
            multiplets_shell::Set{NTuple{4,Int64}} ;
            verbose=false )

    # construct dictionary
    combinations_Gu_muiualpha::Dict{IntIrrep,Vector{Tuple{IntMultiplet,IntMultiplet,IntMultiplet,Int64}}} = Dict()

    G2R::Dict{ NTuple{3,Int64} , Int64 } = Dict{ NTuple{3,Int64} , Int64 }()

    m_u::IntMultiplet = (0,0,0,0)

    # iterate over combinations |i⟩⊗|μ⟩
    for m_i::IntMultiplet  in multiplets_block, 
        m_mu::IntMultiplet in multiplets_shell 

        if verbose 
            println( "m_mu = $m_mu ; m_i = $m_i" )
            println( "============================" )
        end

        (N_i::Int64,I_i::Int64,S_i::Int64,r_i) = m_i
        (N_mu::Int64,I_mu::Int64,S_mu::Int64,r_mu) = m_mu 

        N_u::Int64 = N_i + N_mu
        IIAA_u::Vector{NTuple{2,Int64}} = [(1,1)]
        SS_u::Vector{Int64} = collect(abs(S_mu-S_i):2:(S_mu+S_i))

        for (I_u,A_u) in IIAA_u, 
            α_u in 1:A_u,
            S_u::Int64 in SS_u

            G_u::IntIrrep = (N_u,I_u,S_u)
            if G_u in keys(G2R)
                G2R[G_u] += 1 
                m_u = ( G_u... , G2R[G_u] )
                push!( combinations_Gu_muiualpha[G_u] , (m_mu,m_i,m_u,α_u) )
            else 
                G2R[G_u] = 1
                m_u = ( G_u... , G2R[G_u] )
                combinations_Gu_muiualpha[G_u] = [(m_mu,m_i,m_u,α_u)]
            end
            verbose && println( m_u )
        end
        verbose && println()
    end
    return combinations_Gu_muiualpha
end

# construct and diagonalize step hamiltonian ⟨u|H|v⟩
function construct_and_diagonalize_uHv_totalangularmomentum( 
        multiplets_block::Set{NTuple{4,Int64}} , 
        multiplets_shell::Set{NTuple{4,Int64}} ,
        irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
        hop_symparams::Dict{ Int64 , Matrix{ComplexF64} } ,
        dsum::DSum ,
        lehmann_iaj::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,4} } , # ⟨i||f†_a||j⟩^[n-1]
        lehmann_nuamu::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,4} } , # ⟨ν||f†_a||μ⟩^[n]
        multiplets_a_block::Vector{NTuple{4,Int64}} , 
        multiplets_a_shell::Vector{NTuple{4,Int64}} ;
        conduction_diagonals::Dict{IntIrrep,Vector{Float64}}=Dict{IntIrrep,Vector{Float64}}() ,
        impinfo::Bool=false ,
        verbose::Bool=false )

    # Compute matrix elements 
    #
    #   ⟨u,αu||H||v,αv⟩
    #
    # from the block multiplets 
    #
    #   |i⟩,|j⟩
    #
    # and the local multiplets
    #
    #   |μ⟩,|ν⟩

    verbose && println( "CALCULATION OF <u||H||v> MATRIX\n" )

    ## block-shell combination multiplets 
    combinations_Gu_muiualpha_new = get_combinations_Gu_muiualpha_totalangularmomentum( 
        multiplets_block ,
        multiplets_shell
    )

    # irrep and multiplet combinations 
    #
    # Dict(
    #   G_uv => Dict(
    #       (G_i,G_j) => [(r_i,r_j),...]
    #   )
    # )
    Gu2GmuGi2rmuiualpha::Dict{IntIrrep,Dict{NTuple{2,IntIrrep},Vector{NTuple{4,Int64}}}} = Dict(
        G_u => Dict(
        (G_mu,G_i) => [(m_mu[4],m_i[4],m_u[4],α_u) for (m_mu,m_i,m_u,α_u) in multiplet_combinations_Gu if (m_i[1:3]==G_i && m_mu[1:3]==G_mu)]
            for (G_mu,G_i) in Set( (m_mu[1:3],m_i[1:3]) for (m_mu,m_i,_,_) in multiplet_combinations_Gu )
        )
        for (G_u,multiplet_combinations_Gu) in combinations_Gu_muiualpha_new
    )

    multiplets_a_combs::Vector{NTuple{2,NTuple{4,Int64}}} = [
        (m_a_block,m_a_shell) for m_a_block in multiplets_a_block 
                              for m_a_shell in multiplets_a_shell 
                              if m_a_block[2]==m_a_shell[2] 
    ]
    Ga2amultcombs::Dict{IntIrrep,Vector{NTuple{2,Int64}}} = Dict(
        G_a => [ (m_a_block[4],m_a_shell[4]) for m_a_block in multiplets_a_block 
                                             for m_a_shell in multiplets_a_shell 
                                             if m_a_block[1:3]==G_a ]
        for G_a in Set( m_a_block[1:3] for (m_a_block,_) in multiplets_a_combs )
    )

    # full-sized hblock matrix
    R_uv_max::Int64 = maximum(
        mapreduce( length , + , values(GmuGi2combinations) )
        for (G_u,GmuGi2combinations) in Gu2GmuGi2rmuiualpha
    )
    hblock_full::Array{ComplexF64,3} = zeros(ComplexF64,R_uv_max,R_uv_max,2)

    # construct new irrEU
    irrEU_new::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } = Dict()
    for (G_uv,GmuGi2rmuiualpha) in Gu2GmuGi2rmuiualpha

        # irrep quantum numbers
        N_uv,I_uv,S_uv = G_uv

        # hblock
        R_uv::Int64 = length( combinations_Gu_muiualpha_new[G_uv] )
        @views hblock = hblock_full[1:R_uv,1:R_uv,:]
        hblock .= zero(ComplexF64)

        # iterate through u decompositions only, for diagonal part
        for ((G_mu,G_i),rmuiualphas) in GmuGi2rmuiualpha

            # diagonal part 
            @inbounds for (r_mu,r_i,r_u,α_u) in rmuiualphas
                hblock[r_u,r_u,1] = irrEU[G_i][1][r_i] + conduction_diagonals[G_mu][r_mu]
            end

            # hopping part
            for ((G_nu,G_j),rnujvalphas) in GmuGi2rmuiualpha

                # irrep quantum numbers
                N_i,I_i,S_i = G_i
                N_j,I_j,S_j = G_j
                N_mu,I_mu,S_mu = G_mu
                N_nu,I_nu,S_nu = G_nu

                # hopping allowed?
                hopallowed = (N_nu==(N_mu+1) && N_i==(N_j+1)) 
                if hopallowed

                    # iterate through hopper irreps
                    for (G_a,amultcombs) in Ga2amultcombs

                        # hopping irrep quantum numbers
                        N_a,I_a,S_a = G_a

                        # early discard
                        haskey( lehmann_iaj   , (G_i,G_a,G_j)   ) || continue
                        haskey( lehmann_nuamu , (G_nu,G_a,G_mu) ) || continue

                        # outer multiplicity arrays
                        lehmann_array_iaj   = lehmann_iaj[G_i,G_a,G_j]
                        lehmann_array_nuamu = lehmann_nuamu[G_nu,G_a,G_mu]

                        # clebsch gordan sum D
                        d_spin_and_sign,d_orbital_array = dsum[(G_uv,G_a,G_i,G_j,G_mu,G_nu)]
                        (iszero(length(d_orbital_array)) || iszero(d_spin_and_sign)) && continue

                        # multiplet iteration
                        for (r_mu,r_i,r_u,αu) in rmuiualphas,
                            (r_nu,r_j,r_v,αv) in rnujvalphas

                            # iterate over α,β lehmann outer multiplicities
                            for α in axes(lehmann_array_iaj,1),
                                β in axes(lehmann_array_nuamu,1)

                                d = d_spin_and_sign * d_orbital_array[αu,αv,α,β]
                                iszero(d) && continue

                                # hopping parameter matrix h(Γa)
                                hoparam_matrix = hop_symparams[G_a[3]]

                                # hopping contribution as matrix operation
                                @views hblock[r_u,r_v,2] += d * dot(
                                    lehmann_array_nuamu[β,r_nu,:,r_mu], # automatically complex-conjugated
                                    hoparam_matrix,
                                    lehmann_array_iaj[α,r_i,:,r_j]
                                )

                            end # end of α,β iteration
                        end # end of u,i,mu,v,j,nu multiplet iteration
                    end # end of G_a iteration
                end # end of hopping part
            end # end of v decomposition
        end # end of u decomposition

        # add hopping part with hermitian conjugate
        @inbounds for r_v::Int64 in 1:R_uv, 
                      r_u::Int64 in 1:R_uv
            hblock[r_u,r_v,1] += hblock[r_u,r_v,2] + conj(hblock[r_v,r_u,2])
        end

        # diagonalize
        @views F = eigen( hblock[:,:,1] )
        e, u = real.(F.values), F.vectors

        # insert in irrEU
        irrEU_new[G_uv] = (e,u)

    end 

    minE = minimum([e for (E,U) in values(irrEU_new) for e in E])
    irrEU_new = Dict( G=>(E.-minE,U) for (G,(E,U)) in irrEU_new )

    return ( irrEU_new , combinations_Gu_muiualpha_new )
end

# ==================== #
# FULL NRG CALCULATION #
# ==================== #

# nrg steps
function NRG_totalangularmomentum( 
              label::String ,
              calculation::String ,
              iterations::Int64, 
              cutoff_type::String, 
              cutoff_magnitude::Number,
              L::Float64,
              irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
              multiplets_shell::Set{NTuple{4,Int64}}, 
              cg_o_comb2A::Dict{ NTuple{3,Int64} , Int64 } ,
              dsum::DSum ,
              ksum::KSum ,
              lehmann_muanu::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,4} },
              multiplets_a::Vector{NTuple{4,Int64}} , 
              combinations_Gu_muiualpha::Dict{NTuple{3,Int64}, Vector{Tuple{IntMultiplet,IntMultiplet,IntMultiplet,Int64}}},
              betabar::Float64 ,
              oindex2dimensions::Vector{Int64} ,
              channels_codiagonals::Vector{Dict{Int64,Vector{Float64}}},
              max_spin2::Int64;
              verbose::Bool=false ,
              minmult::Int64=0 , 
              mine::Float64=0.0 ,
              z::Float64=0.0 ,
              spectral::Bool=false ,
              spectral_functions::Dict{String,Dict{IntMultiplet,Matrix{Float64}}}=Dict{String,Dict{IntMultiplet,Matrix{Float64}}}() ,
              K_factor::Float64=2.0 ,
              spectral_method::String="sakai1989",
              spectral_broadening::Float64=1.0 ,
              broadening_distribution::String="gauss",
              dmnrg::Bool=false ,
              dmnrg_run::Int64=1 ,
              shell_dimension::Int64=0 ,
              density_matrices::Vector{Dict{IntIrrep,Matrix{Float64}}}=Dict{IntIrrep,Matrix{Float64}}[] ,
              half_weight_idx::Int64=0 ,
              half_weight_energy::Float64=0.0 ,
              M::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} }=Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} }() ,
              AA::Vector{T}=[] ,
              fsum::FSum=EmptyClebschGordanSum(FSum) ,
              multiplets_atomhop::Vector{NTuple{4,Int64}}=NTuple{4,Int64}[] ,
              scale::Float64=1.0 ,
              compute_impmults::Bool=false ,
              mult2index::Dict{ClearMultiplet,Int64}=Dict{ClearMultiplet,Int64}() ,
              orbital_multiplets::Vector{ClearMultiplet}=ClearMultiplet[] ,
              mm_i::Dict{NTuple{4,Int64},Vector{Float64}}=Dict{NTuple{4,Int64},Vector{Float64}}() ,
              write_spectrum::Bool=false ,
              channels_diagonals::Vector{Dict{IntIrrep,Vector{Float64}}}=[] ,
              compute_selfenergy::Bool=false ,
              Mred_se::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} }=Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} }() ,
              AA_se=[] ,
              impurity_operators::Dict{String,Dict{IntTripleG,Array{ComplexF64,3}}}=Dict{String,Dict{IntTripleG,Array{ComplexF64,3}}}() ,
              spectral_temperature::Float64=0.0 ,
              extra_iterations::Int64=0 ) where {T}

    println( "=============" )
    println( "NRG PROCEDURE" )
    println( "=============" )
    println()

    # thermodynamic information from even and odd iterations
    thermo_even = zeros( Float64 , 0 , 8 )
    thermo_odd  = zeros( Float64 , 0 , 8 )
    thermo_evenodd  = zeros( Float64 , 0 , 8 )

    impspins       = []
    impnums        = []
    impmults       = []

    # spectrum after diagonalization in each step
    spectrum_even = Vector{Dict{IntIrrep,Vector{Float64}}}()
    spectrum_odd  = Vector{Dict{IntIrrep,Vector{Float64}}}()

    # performance at each step
    diagonalization_performances = []
    if spectral
        spectral_performances = []
    end

    # maximum energies and 2S (twice total spin)
    cutoff_eigenenergies = Float64[]
    cutoff_eigenenergy_all_iterations = 0.0
    maximum_spin2s = Int64[]
    maximum_spin2_all_iterations = 0

    # kept state ratios
    kept_discarded_ratios = Float64[]

    # create spectral directory
    isdir("spectral") || mkdir("spectral")

    # dmnrg
    diagonalizers::Vector{Dict{IntIrrep,Matrix{ComplexF64}}} = Dict{IntIrrep,Matrix{ComplexF64}}[]
    combinations_uprima_train = []
    multiplets_kept_train::Vector{Set{IntMultiplet}} = Set{IntMultiplet}[]
    multiplets_disc_train::Vector{Set{IntMultiplet}} = Set{IntMultiplet}[]
    spectrum_train = Dict{IntIrrep,Vector{Float64}}[]
    partition_dmnrg::Float64 = 0.0

    # NRG iterations
    nrg_performance = @timed for n in 2:iterations

        println( "-"^17 )
        @printf "ITERATION n = %3i\n" n
        println( "-"^17 )

        # cutoff
        #
        # apply cutoff
        irrEU_uncut = copy(irrEU)
        (multiplets_block, discarded) = cut_off!( irrEU ; 
                                                  type=cutoff_type , 
                                                  cutoff=cutoff_magnitude , 
                                                  safeguard=true ,
                                                  minmult=minmult ,
                                                  mine=mine ,
                                                  verbose=false ,
                                                  M=!iszero(length(impurity_operators)) ? impurity_operators["particle"] : Dict{IntTripleG,Array{ComplexF64,3}}() )

        # print info
        number_kept_multiplets = length(multiplets_block)
        number_kept_states = sum( (oindex2dimensions[m[2]]*(m[3]+1)) for m in multiplets_block )
        number_discarded_multiplets = length(discarded)
        number_total_multiplets = number_kept_multiplets+number_discarded_multiplets
        kept_discarded_ratio = number_kept_multiplets/number_total_multiplets
        cutoff_eigenenergy = maximum(collect( e for (G,(E,U)) in irrEU for e in E ))
        cutoff_eigenenergy_all_iterations = maximum([ cutoff_eigenenergy , cutoff_eigenenergy_all_iterations ])
        push!( cutoff_eigenenergies , cutoff_eigenenergy )
        println( "CUTOFF" )
        println( "  kept:      $number_kept_multiplets multiplets ($number_kept_states states)" )
        println( "  discarded: $number_discarded_multiplets multiplets" )
        println( "  ratio kept/total: $kept_discarded_ratio, $(number_kept_multiplets)/$(number_total_multiplets) " )
        println( "  cutoff eigenenergy: $cutoff_eigenenergy" )
        push!( kept_discarded_ratios , kept_discarded_ratio )

        if dmnrg && dmnrg_run==1 && !iszero(spectral_temperature) && !iszero(length(discarded))
            nn = n-1
            partition_dmnrg += sum( 
            oindex2dimensions[I]*(S+1)*exp(-irrEU_uncut[N,I,S][1][r]*iterscale(scale,L,nn)/spectral_temperature)*shell_dimension^(iterations-nn)
                for (N,I,S,r) in discarded
            )
        end

        # renormalize by √Λ
        for (G,(E,U)) in irrEU 
            irrEU[G] = ( E.*sqrt(L) , U )
        end

        # hopping parameter
        hop_symparams = Dict( 
            J_o => diagm(ComplexF64.(channels_codiagonals[n-1][J_o])) 
            for J_o in keys(channels_codiagonals[1])
        )
        println( "CONDUCTION COUPLING TERMS" )
        println( "Codiagonals (hoppings)")
        @printf "  %-8s  %-10s  %-s\n" "orbital" "multiplet" "amplitude"
        for (J,hop_matrix) in hop_symparams,
            cartesian_index in CartesianIndices(hop_matrix)

            i,j = Tuple(cartesian_index)
            isapprox(hop_matrix[i,j],0.0) && continue

            @printf "  %-8s  %-3i => %-3i  %-.3f\n" J i j hop_matrix[i,j]

        end
        println( "Diagonals (occupation energies)")
        @printf "  %-8s  %-10s  %-s\n" "irrep" "multiplet" "amplitude"
        for (J_o,J_o_multiplets_diagonals) in channels_diagonals[n],
            (r_o,r_o_multiplet_diagonal) in enumerate(J_o_multiplets_diagonals)

            @printf "  %-8s  %-3i => %-.3f\n" J_o r_o r_o_multiplet_diagonal

        end

        # compute new iaj lehmann amplitudes
        #
        #   <i||f†_a||j>_α
        #
        lehmann_iaj = compute_lehmann_iaj( 
            collect(multiplets_a), 
            ksum ,
            lehmann_muanu ,
            combinations_Gu_muiualpha,
            irrEU ,
            cg_o_comb2A
        )

        # construct and diagonalize ⟨u||H||v⟩
        diagonalization_performance = @timed (irrEU,combinations_Gu_muiualpha) = construct_and_diagonalize_uHv_totalangularmomentum( 
                    multiplets_block, 
                    multiplets_shell,
                    irrEU, 
                    hop_symparams, 
                    dsum,
                    lehmann_iaj,
                    lehmann_muanu,
                    multiplets_a,
                    multiplets_a;
                    conduction_diagonals=channels_diagonals[1]
        )

        # information
        maximum_spin2 = maximum(collect( G[3] for (G,(E,U)) in irrEU ))
        if maximum_spin2>max_spin2
            error( """
                ERROR: max_spin2 too low!
                max_spin2 = $max_spin2
                maximum spin2 found in calculation = $maximum_spin2
            """)
        end
        maximum_spin2_all_iterations = maximum([ maximum_spin2 , maximum_spin2_all_iterations ])
        maximum_eigenenergy_after_diagonalization = maximum([ maximum(E) for (_,(E,_)) in irrEU ])
        println( "DIAGONALIZATION" )
        println( "  time: $(diagonalization_performance.time)s" )
        println( "  memory: $(diagonalization_performance.bytes*10^-6)Mb" )
        println( "  garbage collection time: $(diagonalization_performance.gctime)s" )
        println( "  maximum 2S needed: $maximum_spin2")
        println( "  maximum eigenenergy obtained: $maximum_eigenenergy_after_diagonalization" )

        # store performance
        push!( diagonalization_performances , diagonalization_performance )

        # dmnrg information
        if dmnrg && dmnrg_run==1

            # info gathered from cutoff to previous state
            push!( multiplets_kept_train , multiplets_block )
            push!( multiplets_disc_train , Set(discarded) )

            # info gathered from diagonalization in this step
            push!( combinations_uprima_train , combinations_uprima )
            push!( diagonalizers , Dict( G=>U for (G,(_,U)) in irrEU ) )
            push!( spectrum_train , Dict( G=>E for (G,(E,_)) in irrEU ) )
        end

        # spectrum 
        if write_spectrum
            (n%2==0)  && save_spectrum!( spectrum_even , irrEU )
            (n%2!==0) && save_spectrum!( spectrum_odd  , irrEU )
        end

        # impurity info
        if compute_impmults
            mm_i = imp_mults( irrEU ,
                              oindex2dimensions ,
                              combinations_uprima ,
                              mm_i )
            m_imp::Vector{Float64} = mult_thermo( irrEU ,
                                 betabar ,
                                 oindex2dimensions ,
                                 mm_i )
            push!( impmults , m_imp )
        end

        # thermodynamics 
        #
        # iteration temperature
        temperature = compute_temperature_newdiscretization(n,L,betabar,scale)
        # compute thermodynamic variables
        magnetic_susceptibility = compute_magnetic_susceptibility( irrEU , betabar , oindex2dimensions )
        entropy = compute_entropy( irrEU , betabar , oindex2dimensions )
        heat_capacity = compute_heat_capacity( irrEU , betabar , oindex2dimensions )
        free_energy = compute_free_energy( irrEU , betabar , oindex2dimensions )
        number_particles = compute_average_number_of_particles( irrEU , betabar , oindex2dimensions )
        energy = compute_energy( irrEU , betabar , oindex2dimensions )
        partition_function = compute_partition_function( irrEU , betabar , oindex2dimensions )
        thermodynamic_matrix = [temperature magnetic_susceptibility entropy heat_capacity free_energy number_particles energy partition_function]
        # store even/odd therodynamic results
        if n%2==0
            thermo_even = vcat( thermo_even , thermodynamic_matrix )
        else
            thermo_odd  = vcat( thermo_odd , thermodynamic_matrix )
        end
        thermo_evenodd = vcat( thermo_evenodd , thermodynamic_matrix )
        # information
        println( "THERMODYNAMICS" )
        @printf "  %s = %.3e\n" "temperature" temperature
        @printf "  %s = %.3f\n" "magnetic susceptibility" magnetic_susceptibility
        @printf "  %s = %.3f\n" "entropy" entropy
        @printf "  %s = %.3f\n" "heat capacity" heat_capacity
        @printf "  %s = %.3f\n" "free energy" free_energy
        @printf "  %s = %i\n"   "average number of particles" number_particles
        @printf "  %s = %.3f\n" "energy" energy
        @printf "  %s = %.3e\n" "Z" partition_function

        # spectral function calculation
        if spectral 

            spectral_performance = @timed begin
                impurity_operators["particle"] = update_operator( impurity_operators["particle"], 
                                                                  collect(multiplets_atomhop) ,
                                                                  Karray_orbital ,
                                                                  Karray_spin ,
                                                                  combinations_uprima ,
                                                                  irrEU )
                multiplets_kept = Set{IntMultiplet}()
                multiplets_discarded = Set{IntMultiplet}()
                if dmnrg_run==2
                    irrEU_copy = copy(irrEU)
                    cut = n==iterations ? sum( (iszero(e) ? 1 : 0) for (_,(E,_)) in irrEU for e in E ) : cutoff_magnitude
                    (multiplets_kept,multiplets_discarded) = cut_off!( irrEU ; 
                                                                       type=cutoff_type , 
                                                                       cutoff=cut , 
                                                                       safeguard=(n!==iterations) ,
                                                                       minmult=minmult ,
                                                                       mine=mine )
                    irrEU = irrEU_copy
                end
                add_correlation_contribution!(
                    spectral_functions["spectral"],
                    impurity_operators["particle"],
                    impurity_operators["particle"],
                    oindex2dimensions,
                    irrEU,
                    n ,
                    broadening_distribution ,
                    spectral_broadening ,
                    iterscale(scale,L,n) ,
                    K_factor ; 
                    correlation_type="spectral",
                    T=spectral_temperature ,
                    limit_shell = n==iterations ,
                    extra_iterations=extra_iterations ,
                    density_matrix = dmnrg_run==2 ? density_matrices[n] : Dict{IntIrrep,Matrix{Float64}}() ,
                    multiplets_kept=collect(multiplets_kept) ,
                    multiplets_discarded=collect(multiplets_discarded) ,
                    L=L,
                    scale=scale,
                    half_weight_idx=half_weight_idx,
                    half_weight_energy=half_weight_energy
                )
            end

            # information
            maximum_irrep_spin2 = irreps -> maximum((irreps[1][3],irreps[2][3],irreps[3][3]))
            maximum_spin2 = maximum([ maximum_irrep_spin2(irreps) for irreps in keys(impurity_operators["particle"]) ])
            maximum_spin2_all_iterations = maximum([ maximum_spin2_all_iterations , maximum_spin2 ])
            println( "EXCITATION MATRIX CALCULATION" )
            println( "  time: $(spectral_performance.time)s" )
            println( "  memory: $(spectral_performance.bytes*10^-6)Mb" )
            println( "  garbage collection time: $(spectral_performance.gctime)s" )
            println( "  maximum 2S needed: $maximum_spin2")
            push!( spectral_performances , spectral_performance )

        end

        println()
    end # NRG iterations finished
    println()


    # average even and odd thermodynamic results
    #
    # reverse thermodynamic matrices for interpolation
    reverse!( thermo_even , dims=1 )
    reverse!( thermo_odd  , dims=1 )
    reverse!( thermo_evenodd , dims=1 )
    # collect temperatures from even and odd
    temperatures_even = thermo_even[:,1]
    temperatures_odd  = thermo_odd[:,1]
    temperatures_evenodd = sort(vcat(temperatures_even,temperatures_odd))
    # interpolate even and odd data to new temperatures
    thermo_even_interpolated = interpolate_thermo_matrix( thermo_even , temperatures_evenodd )
    thermo_odd_interpolated  = interpolate_thermo_matrix( thermo_odd  , temperatures_evenodd )
    # average even and odd
    thermo_average = copy(thermo_even_interpolated)
    thermo_average[:,2:end] = 0.5*( thermo_even_interpolated[:,2:end] + thermo_odd_interpolated[:,2:end] )

    spectral && save_correlation_spectral_decomposition(spectral_functions,label,z)

    # reduced density matrices for all steps
    if dmnrg && dmnrg_run==1

        if iszero(spectral_temperature)

            Z = sum( (iszero(e) ? oindex2dimensions[I]*(S+1) : 0.0)
                     for ((_,I,S),(E,_)) in irrEU
                     for e in E )
            push!( 
                density_matrices ,
                Dict(
                    (N,I,S)=>diagm([ ( iszero(e) ? 1.0/Z : 0.0 ) for e in E ])
                    for ((N,I,S),(E,_)) in irrEU
                    if iszero(E[1])
                )
            )
            for i in reverse(eachindex(diagonalizers))
                diagonalizer = diagonalizers[i]
                combinations_uprima = combinations_uprima_train[i]
                multiplets_kept = multiplets_kept_train[i]
                insert!( 
                    density_matrices ,
                    1 ,
                    backwards_dm_reduced_T0(density_matrices[1],
                                            diagonalizer,
                                            combinations_uprima,
                                            oindex2dimensions) 
                )
            end
            #for (i,dm) in enumerate(density_matrices)
            #    println( "DM for iteration $(i) with (i+1)=$(i+1)" )
            #    @show keys(dm)
            #    trace = sum( tr(d)*oindex2dimensions[I]*(S+1) for ((_,I,S),d) in dm )
            #    tr_plus = sum(  (N>(i+1) ? tr(d) : 0.0) for ((N,_,_),d) in dm )
            #    tr_minus = sum( (N<(i+1) ? tr(d) : 0.0) for ((N,_,_),d) in dm )
            #    @show tr_plus,tr_minus
            #    @show trace
            #    println()
            #end
        else # T!==0

            # add last contribution to partition function
            partition_dmnrg += sum( 
                oindex2dimensions[I]*(S+1)*exp(-e*iterscale(scale,L,iterations)/spectral_temperature)
                for ((_,I,S),(E,_)) in irrEU
                for e in E
            )
            Z_last = sum( 
                oindex2dimensions[I]*(S+1)*exp(-e*iterscale(scale,L,iterations)/spectral_temperature)
                for ((_,I,S),(E,_)) in irrEU
                for e in E
            )
            Z = partition_dmnrg
            @show Z
            @show Z_last
            @show 

            # set first density matrix for n=iterations
            push!( 
                density_matrices ,
                Dict(
                    (N,I,S)=>diagm([ exp(-e*iterscale(scale,L,iterations)/spectral_temperature)/Z for e in E ])
                    for ((N,I,S),(E,_)) in irrEU
                )
            )

            # weights of density matrices for each iteration
            partial_dm_weight_train = Float64[]
            for i in eachindex(spectrum_train)
                n = i+1
                if n!==iterations
                    multiplets_disc = multiplets_disc_train[i+1]
                    if length(multiplets_disc)==0
                        push!(partial_dm_weight_train,0.0)
                        continue
                    end
                    spectrum = spectrum_train[i]
                    #@show multiplets_disc
                    #print_dict(spectrum)
                    push!( 
                        partial_dm_weight_train ,
                        shell_dimension^(iterations-n)/
                        Z*
                        sum(oindex2dimensions[I]*(S+1)*exp(-spectrum[N,I,S][r]*iterscale(scale,L,n)/spectral_temperature) for (N,I,S,r) in multiplets_disc)
                    )
                else
                    spectrum = spectrum_train[end]
                    push!( 
                        partial_dm_weight_train ,
                        sum(oindex2dimensions[I]*(S+1)*exp(-e*iterscale(scale,L,n)/spectral_temperature) for ((_,I,S),E) in spectrum for e in E)/Z
                    )
                end
            end
            @show partial_dm_weight_train
            max_weight, max_weight_idx = findmax(partial_dm_weight_train)
            _, half_weight_idx = findmin( abs2.(max_weight/2 .- partial_dm_weight_train[max_weight_idx:end]) )
            half_weight_idx += max_weight_idx
            half_weight = partial_dm_weight_train[half_weight_idx]
            @show max_weight, max_weight_idx
            @show half_weight, half_weight_idx
            @show iterscale(scale,L,max_weight_idx+1)
            @show iterscale(scale,L,half_weight_idx+1)
            @show spectral_temperature
            half_weight_idx = max_weight_idx
            half_weight_energy = iterscale(scale,L,half_weight_idx+1)

            for i in reverse(eachindex(diagonalizers))
                iteration = i+1
                diagonalizer = diagonalizers[i]
                combinations_uprima = combinations_uprima_train[i]
                multiplets_kept = multiplets_kept_train[i]
                multiplets_disc = multiplets_disc_train[i]
                spectrum = i-1<=0 ? nothing : spectrum_train[i-1]
                @show i-1,length(spectrum_train)
                insert!( 
                    density_matrices ,
                    1 ,
                    i-1<=0 ? Dict{IntIrrep,Matrix{ComplexF64}}() : backwards_dm_reduced_Tnonzero(
                        density_matrices[1],
                        diagonalizer,
                        combinations_uprima,
                        oindex2dimensions,
                        multiplets_kept,
                        multiplets_disc,
                        spectral_temperature,
                        Z,
                        shell_dimension,
                        iterations,
                        iteration,
                        iterscale(scale,L,iteration),
                        spectrum
                    ) 
                )
                #for (G,m) in density_matrices[1]
                #    @show G
                #    @show m[1,1]
                #    println()
                #end
            end
        end
    end

    # print summary information
    println( "=========================" )
    println( "SUMMARY OF NRG ITERATIONS" )
    println( "=========================" )
    # energy, kept states and spin2
    cutoff_eigenenergy_average = sum(cutoff_eigenenergies)/length(cutoff_eigenenergies)
    kept_discarded_ratio_average = sum(kept_discarded_ratios)/length(kept_discarded_ratios)
    println( "CHECK" )
    println( "  maximum 2S needed: $maximum_spin2_all_iterations" )
    println( "  average cutoff eigenenergy: $cutoff_eigenenergy_average" )
    println( "  average kept/discarded ratio: $kept_discarded_ratio_average" )
    println( "PERFORMANCE SUMMARY" )
    println( "  Diagonalization" )
    print_performance_summary( diagonalization_performances )
    if spectral
        println( "  Excitation matrix calculation" )
        print_performance_summary( spectral_performances )
    end
    println( "  Global" )
    println( "    time: $(nrg_performance.time)s" )
    println( "    memory: $(nrg_performance.bytes*10^-6)Mb" )
    println( "    garbage collection time: $(nrg_performance.gctime)s" )

    # save data to file (spectral is saved from within the function)
    #
    # thermodynamics
    if !spectral

        # directory
        isdir("thermodata") || mkdir("thermodata")

        # thermodynamic data for this given value of z
        write_thermodata_onez( thermo_average , calculation , label , z )
        write_thermodata_onez( thermo_evenodd , calculation , label*"_evenodd" , z )

        # impurity contribution (diff)
        if calculation=="IMP"
            thermo_clean_filename = thermo_filename_one_z( label , "clean" , z )
            println()
            if length(glob(thermo_clean_filename))!==0
                println( "Saving thermodynamic impurity contribution to $(thermo_filename_one_z(label,"diff",z))..." )
                write_thermodiff( label , z )
                write_thermodiff( label*"_evenodd" , z )
            end
        end

    end
    # impurity properties 
    if compute_impmults

        println( "Saving impurity projections..." )
        # directory
        isdir("impurityprojections") || mkdir("impurityprojections")

        if calculation=="IMP" 
            write_impurity_info( impmults , orbital_multiplets , mult2index , label , z )
        end

    end
    # spectra at each step
    if write_spectrum
        write_nrg_spectra( spectrum_even , spectrum_odd )
    end

    return ( 
        thermo = thermo_average ,
        diagonalization_performances = diagonalization_performances ,
        impmults = compute_impmults ? impmults : nothing ,
        spectral_performances = spectral ? spectral_performances : nothing ,
        density_matrices = density_matrices ,
        half_weight_idx = half_weight_idx ,
        half_weight_energy = half_weight_energy
    )
end

# full nrg run: initialize system and run nrg iterations
function nrg_full_totalangularmomentum( 
            label::String ,
            calculation::String ,
            L::Float64 ,
            iterations::Int64 ,
            cutoff_type::String ,
            cutoff_magnitude ,
            multiplets_dir::String ,
            impurity_config::Dict{Float64,Int64} ,
            shell_config::Dict{Float64,Int64} ,
            epsilon_symparams::Dict{ Float64 , Vector{ComplexF64} } ,
            u_symparams::Dict{ Tuple{String,Float64} , Matrix{ComplexF64} } ,
            hop_symparams::Dict{ Float64 , Matrix{ComplexF64} } ;
            distributed::Bool=false ,
            z::Float64=0.0 ,
            max_J2::Int64=10 , # for J=1/2, otherwise pick larger
            channels_dos::Dict{ Float64 , Vector{Function} }=Dict{ Float64 , Vector{Function} }() ,
            discretization::String=discretization_default ,
            tridiagonalization::String=tridiagonalization_default ,
            enforce_particle_hole_symmetry::Bool=true,
            mine::Float64=0.0 ,
            betabar::Float64=1.0 ,
            spectral::Bool=false ,
            spectral_broadening::Float64=0.5 ,
            broadening_distribution::String="loggaussian" ,
            K_factor::Float64=2.0 ,
            orbitalresolved::Bool=false ,
            spectral_temperature::Float64=0.0 ,
            extra_iterations::Int64=0 ,
            dmnrg::Bool=false ,
            compute_impmults::Bool=false ,
            scale_asymptotic::Bool=true ,
            band_width::Float64=1.0 ) 

    # impmults only with imp
    compute_impmults = compute_impmults && (calculation=="IMP")

    println( "********************************" )
    println( "Full NRG calculation with z=$(z)" )
    println( "********************************" )
    println()

    if (spectral && calculation=="CLEAN") 
        error( "Calculation must be IMP for computing the spectral function" )
        return nothing 
    end

    # orbital irreps present in the atom
    atom_orbital_irreps::Vector{Float64} = collect(keys(impurity_config))

    println( "====================" )
    println( "SETUP AND PARAMETERS" )
    println( "====================" )
    @show calculation
    @show distributed 
    @show discretization
    @show L
    @show iterations
    @show z
    @show betabar
    @show cutoff_type
    @show cutoff_magnitude
    @show max_J2
    @show spectral
    println()
    println( "OCCUPATION ENERGIES" )
    print_dict( epsilon_symparams ) 
    println( "COULOMB PARAMETERS" )
    print_dict( u_symparams ) 
    println( "HYBRIDIZATION PARAMETERS" )
    print_dict( hop_symparams )
    println()

    #   ==========================   #
    #%% SYMMETRY-RELATED VARIABLES %%#
    #   ==========================   #

    # orbital symmetry
    #(oirreps::Vector{String},
    # oirreps2indices::Dict{String,Int64},
    # oirreps2dimensions::Dict{String,Int64},
    # oindex2dimensions::Vector{Int64},
    # cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}}) = get_cg_o_info_nonsimple( cg_o_dir , atom_orbital_irreps )
    identityrep::String = "A"
    oirreps::Vector{String} = [identityrep]
    oirreps2indices::Dict{String,Int64} = Dict(identityrep=>1)
    oirreps2dimensions::Dict{String,Int64} = Dict(identityrep=>1)
    oindex2dimensions::Dict{Int64,Int64} = Dict(1=>1)
    cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}} = Dict(
        (1,1,1) => ones(ComplexF64,1,1,1,1)
    )

    # total angular momentum symmetry (handled as spin symmetry)
    #cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} = Dict( (0,0,0)=>[1.0;;;] )
    # SU(2) CLEBSCH-GORDAN COEFFICIENTS NOT STORED

    # for dmnrg
    shell_dimension = reduce( * , [jdim(J)*R for (J,R) in shell_config] )

    # G1 x G2 = n1×G1' + n2×G2' + ...
    # (G1,G2) => [G1',G2',...]
    cg_o_fullmatint_keys = keys(cg_o_fullmatint)
    keys_as_dict_o::Dict{NTuple{2,Int64},Vector{NTuple{2,Int64}}} = Dict(
        (I1,I2)=>collect(Set((I_3,size(arr,1)) for ((I_1,I_2,I_3),arr) in cg_o_fullmatint if (I_1,I_2)==(I1,I2)))
        for (I1,I2,_) in cg_o_fullmatint_keys 
    )
    cg_o_comb2A::Dict{ NTuple{3,Int64} , Int64 } = Dict(
        (I1,I2,I3)=>size(arr,1)
        for ((I1,I2,I3),arr) in cg_o_fullmatint
    )
    # REMOVED SPIN PART
    #cg_s_fullmatint_keys = keys(cg_s_fullmatint)
    #keys_as_dict_s::Dict{NTuple{2,Int64},Vector{Int64}} = Dict(
    #    (S1,S2)=>collect(Set(x[3] for x in cg_s_fullmatint_keys if x[1:2]==(S1,S2)))
    #    for (S1,S2,_) in cg_s_fullmatint_keys 
    #)

    #   ===========   #
    #%% ATOMIC PART %%#
    #   ===========   #
    println()
    println( ":::::::::::::::::::" )
    println( "--- ATOMIC PART ---" )
    println( ":::::::::::::::::::" )
    println()

    #   --------------   #
    #%% discretization %%#
    #   --------------   #

    # default behavior: eta=x->1/2
    if length(channels_dos)==0 
        channels_dos = Dict{Float64,Vector{Function}}( 
            J=>Function[x->0.5 for i in 1:size(hop_matrix,1)]
            for (J,hop_matrix) in hop_symparams 
        )
    end

    # channel coupling parameters
    channels_tridiagonal = discretize_bands( channels_dos ,
                                             L ,
                                             z , 
                                             iterations ;
                                             discretization=discretization ,
                                             tridiagonalization=tridiagonalization ,
                                             enforce_particle_hole_symmetry=enforce_particle_hole_symmetry )
    channels_tridiagonal_int::Dict{Int64,Vector{Tuple{Vector{Float64},Vector{Float64}}}} = Dict(
        doublej(J)=>v 
        for (J,v) in channels_tridiagonal
    )
    # [ n -> { I_o => [ r_o -> diagonal_element ] } ]
    channels_codiagonals::Vector{Dict{Int64,Vector{Float64}}} = [
        # n (iterations) loop
        Dict(# J2 loop
            J2 => [# r_o loop
                r_o_multiplet_codiagonals[n]
                for (r_o,(_,r_o_multiplet_codiagonals)) in enumerate(J2_multiplets_couplings)
            ]
            for (J2,J2_multiplets_couplings) in channels_tridiagonal_int
        )
        for n in 1:(iterations-1)
    ]

    # scaling:
    #
    #   H0 = L^(-0.5) * ( Himp + Hhyb + Hshell0 )
    #
    #   - L^(-0.5) factor already included in rescale function (below)
    #     for Himp and Hhyb
    #   - for Hshell0, the factor is already included in the rescaled constants
    # 
    # in the present method, the scale, which otherwise corresponds to the 
    # largest among the first asymptotic codiagonal coupling terms, is just 
    # 1.0 because for the sake of generality we do not assume the asymptotic
    # form (yet?).
    scale::Float64 = band_width
    if scale_asymptotic
        codiagonals_first = collect(values(channels_codiagonals[10]))[1][1]
        scale *= codiagonals_first
        channels_codiagonals = [
            # n (iterations) loop
            Dict(# J2 loop
                J2 => [# r_o loop
                    r_o_multiplet_codiagonals[n]/codiagonals_first
                    for (r_o,(_,r_o_multiplet_codiagonals)) in enumerate(J2_multiplets_couplings)
                ]
                for (J2,J2_multiplets_couplings) in channels_tridiagonal_int
            )
            for n in 1:(iterations-1)
        ]
    end
    println( "Scale: D(× asymptotic hopping) = $scale")
    println()

    # adapt iterations to T-dependent spectral function calculation
    if spectral && !iszero(spectral_temperature) && !dmnrg
        println()
        #_,n_limit = findmin(abs.([iterscale(scale,L,n) for n in 0:iterations if iterscale(scale,L,n)>betabar*spectral_temperature].-betabar*spectral_temperature))
        ispositive(x) = x>0
        _,n_limit = findmin(filter( ispositive , [iterscale(scale,L,n) for n in 0:iterations] .- betabar*spectral_temperature ))
        n_limit -= 1
        if n_limit>iterations
            println( "WARNING: Not enough iterations to reach temperature limit." )
        end
        println( "Defined T: $spectral_temperature")
        println( "Iterations cut from $iterations to $n_limit for the effective chain to reach the energy scale $(iterscale(scale,L,n_limit))." )
        println()
        iterations = n_limit>iterations ? iterations : n_limit
    end

    #   ------------------- #
    #%% rescaled parameters #
    #   ------------------- #
    oindices2irreps = Dict( v=>k for (k,v) in oirreps2indices )
    epsilon_symparams = Dict( k=>@.rescale(v,L,z,scale) for (k,v) in epsilon_symparams )
    u_symparams       = Dict( k=>@.rescale(v,L,z,scale) for (k,v) in u_symparams )
    A_L = 0.5*(L+1)/(L-1)*log(L)  # correction factor
    hopscale = discretization=="yoshida1990" ? scale/sqrt(A_L) : scale
    hop_symparams_int = Dict{Int64,Matrix{ComplexF64}}( doublej(J)=>(@.rescale(v,L,z,hopscale)) for (J,v) in hop_symparams )
    println( "RESCALED PARAMETERS FOR H0" )
    @show epsilon_symparams 
    @show u_symparams 
    @show hop_symparams
    println()

    
    #   ------------------------------- #
    #%% symstates, basis and multiplets #
    #   ------------------------------- #
    symstates_atom_noint::Dict{Tuple{Int64,String,Float64,Int64,Float64,Int64},State} = Dict()
    multiplets_atom_noint::Set{Tuple{Int64,String,Float64,Int64}} = Set()
    multiplets_a_atom_noint::Set{Tuple{Int64,String,Float64,Int64}} = Set()
    mult2index::Dict{ClearMultiplet,Int64} = Dict()
    orbital_multiplets::Vector{Tuple{Int64,String,Float64,Int64}} = []
    if calculation=="IMP"
        symstates_atom_noint,
        basis_atom,
        multiplets_atom_noint,
        multiplets_a_atom_noint = 
            get_symstates_basis_multiplets_totalangularmomentum( 
                    impurity_config,
                    multiplets_dir;
                    verbose=true )
        orbital_multiplets = ordered_multiplets(multiplets_atom_noint)
        mult2index = Dict( m=>i for (i,m) in enumerate(orbital_multiplets))
        multiplets_atom::Set{NTuple{4,Int64}} = multiplets2int( multiplets_atom_noint , 
                                                                oirreps2indices )
        multiplets_a_atom::Set{NTuple{4,Int64}} = multiplets2int( multiplets_a_atom_noint , 
                                                                  oirreps2indices )
    else 
        multiplets_atom_noint = Set([(0,joirrep,0.0,1)]) 
        multiplets_atom = multiplets2int(multiplets_atom_noint,
                                         oirreps2indices)
        multiplets_a_atom = multiplets_atom
    end

    #   ------------------------ #
    #%% reduced pcg coefficients #
    #   ------------------------ #
    lehmann_iaj::IntIrrepPCGNS = 
        calculation=="CLEAN" ? 
        IntIrrepPCG() :
        get_lehmann_reduced_totalangularmomentum( 
            basis_atom ,
            symstates_atom_noint::ClearSymstateDict ,
            multiplets_atom::IntMultipletSet ;
            verbose=false )::IntIrrepPCGNS
    print_dict(lehmann_iaj)

    #   ------------- #
    #%% impurity atom #
    #   ------------- #
    if calculation=="IMP"

        # operators
        epsilon::Operator{typeof(basis_atom)} = epsilon_sym_totalangularmomentum( symstates_atom_noint , epsilon_symparams ; verbose=false )
        coulomb::Operator{typeof(basis_atom)} = u_sym_totalangularmomentum( symstates_atom_noint , u_symparams ; verbose=false )

        # hamiltonian 
        H::Operator{typeof(basis_atom)} = epsilon + coulomb 

    end

    #   -----   #
    #%% irreu %%#
    #   -----   #
    irrEU_clear::Dict{ Tuple{Int64,String,Float64} , Tuple{Vector{Float64},Matrix{ComplexF64}} } = Dict()
    if calculation=="IMP"
        irrEU_clear = 
            get_irrEU_initial(symstates_atom_noint,H;verbose=true)::Dict{ Tuple{Int64,String,Float64} , Tuple{Vector{Float64},Matrix{ComplexF64}} }
    elseif calculation=="CLEAN" 
        irrEU_clear = 
            get_irrEU_initial(identityrep,oirreps2indices)::Dict{ Tuple{Int64,String,Float64} , Tuple{Vector{Float64},Matrix{ComplexF64}} }
    end
    println( "------------------------------------" )
    println( "IMPURITY SPECTRUM" )
    println()
    println( "Before rescale:" )
    println()
    print_spectrum( Dict((G,(E.*(scale*sqrt(L)),U)) for (G,(E,U)) in irrEU_clear) )
    println()
    println( "After rescale:" )
    println()
    print_spectrum( irrEU_clear )
    println( "------------------------------------" )
    println()
    irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} } = 
        irrEU2int( irrEU_clear , oirreps2indices ) 


    #   ------------------------   #
    #%% impurity quantum numbers %%#
    #   ------------------------   #
    mm_i::Dict{IntMultiplet,Vector{Float64}} = Dict()
    if compute_impmults
        mm_i,m_imp::Vector{Float64} = 
            setup_impmultinfo( multiplets_atom ,
                               irrEU ,
                               betabar ,
                               oindex2dimensions )
        println( "IMPURITY COMPOSITION" )
        println()
        @show mm_i
        println()
        @show m_imp
        println()
    end


    #   ==================   #
    #%% SHELL CONSTRUCTION %%#
    #   ==================   #

    println()
    println( ":::::::::::::::::::::::" )
    println( "--- SHELL STRUCTURE ---" )
    println( ":::::::::::::::::::::::" )
    println()

    #   ------------------------------- #
    #%% symstates, basis and multiplets #
    #   ------------------------------- #
    symstates_shell_noint::Dict{Tuple{Int64,String,Float64,Int64,Float64,Int64},State},
    basis_shell,
    multiplets_shell_noint::Set{Tuple{Int64,String,Float64,Int64}},
    multiplets_a_shell_noint::Set{Tuple{Int64,String,Float64,Int64}} = 
        get_symstates_basis_multiplets_totalangularmomentum(
                shell_config,
                multiplets_dir;
                verbose=true )
    multiplets_shell::Set{NTuple{4,Int64}} = multiplets2int( multiplets_shell_noint , 
                                       oirreps2indices )
    multiplets_a_shell::Set{NTuple{4,Int64}} = multiplets2int( multiplets_a_shell_noint , 
                                         oirreps2indices )

    #   ------------------------------------------------------   #
    #%% shell hamiltonian for particle-hole asymmetric systems %%#
    #   ------------------------------------------------------   #
    #
    # occupation operator for the various shell orbitals
    #   { orbital_irrep_int => [{ G_mu=>[N_I_mu_r_mu] }] }
    orbitals2diagonalcounter = Dict{ Int64 , Vector{Dict{ IntIrrep , Vector{Float64} }} }()
    for (J,orbital_multiplets_diagonals) in channels_tridiagonal

        # double J
        J2 = doublej(J)

        # occupations (irrEU format) for each orbital multiplet
        # belonging to the orbital_irrep
        orbital_multiplets_occupations = []

        # iterate through orbital multiplets
        for (r_o,_) in enumerate(orbital_multiplets_diagonals) 

            # one-electron orbital multiplet for which to compute the occupations
            orbital = (J,r_o)

            # operator for the chosen one-electron orbital
            counter_operator = electron_counter_sym_totalangularmomentum( 
                symstates_shell_noint ,
                orbital 
            )

            # diagonalize operator
            orbital_count_diagonalization = irrEU2int(get_irrEU_initial( symstates_shell_noint , counter_operator ),oirreps2indices)

            # introduce eigenvalues (number of particles) into dictionary
            push!( orbital_multiplets_occupations , Dict( G_mu=>N_mu for (G_mu,(N_mu,_)) in orbital_count_diagonalization) )

        end

        # store multiplet occupations (irrEU format) in main dictionary
        orbitals2diagonalcounter[J2] = copy(orbital_multiplets_occupations)

    end
    # occupations for the shell symstates
    #   { G_mu => [ r_mu -> [ J_o => [ J_o -> number_of_particles] ] ] }
    G2R_mu = Dict( 
        G_mu => length([m for m in multiplets_shell if m[1:3]==G_mu]) 
        for G_mu in Set( m[1:3] for m in multiplets_shell )
    )
    shell_sym_occupations = Dict{ IntIrrep , Vector{Dict{Int64,Vector{Float64}}} }(
        # (G_mu,R_mu) iteration
        G_mu => [# r_mu in 1:R_mu iteration
                    Dict(# (I_o,I_o_multiplets_occupations) iteration
                        J2_o => [# r_o_occupations iteration
                            r_o_occupations[G_mu][r_mu]
                            for r_o_occupations in J2_o_multiplets_occupations
                        ] 
                        for (J2_o,J2_o_multiplets_occupations) in orbitals2diagonalcounter
                    )
                    for r_mu in 1:R_mu
                ]
        for (G_mu,R_mu) in G2R_mu
    )
    # diagonal shell parameters
    #   [ iteration -> { G_mu => [ r_mu -> [ I_o -> [ r_o -> diagonal_element ]]] } ]
    channels_diagonals::Vector{Dict{ IntIrrep , Vector{Float64} }} = [
        # n (iteration) loop 
        Dict(# G_mu loop
            G_mu => [# r_mu loop
                        sum(# J_o loop
                            sum(# r_o loop
                                r_o_couplings[n]*r_mu_multiplet_occupations[J_o][r_o] 
                                for (r_o,(r_o_couplings,_)) in enumerate(J_o_multiplets_couplings)
                            )
                            for (J_o,J_o_multiplets_couplings) in channels_tridiagonal_int
                        )
                        for r_mu_multiplet_occupations in G_mu_multiplets_occupations
                    ]
            for (G_mu,G_mu_multiplets_occupations) in shell_sym_occupations
        )
        for n in 1:iterations
    ]

    #   ------------------------ #
    #%% reduced pcg coefficients #
    #   ------------------------ #
    lehmann_muanu::Dict{NTuple{3,NTuple{3,Int64}},Array{ComplexF64,4}} = get_lehmann_reduced_totalangularmomentum( 
        basis_shell,
        symstates_shell_noint,
        multiplets_shell;
        verbose=true
    )


    #   ================================   #
    #%% COUPLING ATOM TO INNERMOST SHELL %%#
    #   ================================   #

    println()
    println( "::::::::::::::::::::::::::::::::::::::::" )
    println( "--- COUPLING ATOM TO INNERMOST SHELL ---" )
    println( "::::::::::::::::::::::::::::::::::::::::" )
    println()

    #   ------------------------------ #
    #%% precompute clebsch-gordan sums #
    #   ------------------------------ #
    #
    dsum,ksum,fsum = ClebschGordanSums(
        maximum([m[3] for m in vcat(collect(multiplets_a_shell),collect(multiplets_a_atom))]), # maximum_J2_onelectron
        max_J2,
        maximum([m[3] for m in vcat(collect(multiplets_shell),collect(multiplets_atom))]); # maximum_J2_onelectron
        verbose=false 
    )

    #   -------- #
    #%% spectral #
    #   -------- #
    impurity_operators = Dict{String,Dict{IntTripleG,Array{ComplexF64,3}}}()
    spectral_functions = Dict{String,Dict{IntMultiplet,Matrix{Float64}}}()
    if spectral

        impurity_operators["particle"] = lehmann_iaj
        GG_a  = Set(G_a for (_,G_a,_) in keys(lehmann_iaj))
        G2R_a = Dict( 
            G_a=>size(mat,2)
            for G_a in GG_a
            for ((_,Ga,_),mat) in impurity_operators["particle"]
            if G_a==Ga 
        )

        extra_iterations = (dmnrg || iszero(spectral_temperature)) ? 0 : extra_iterations
        spectral_functions = Dict{String,Dict{IntMultiplet,Matrix{Float64}}}(
            "spectral"=>Dict(
            (G_a...,r_a)=>reduce(vcat,sort([[K_factor*sign*iterscale(scale,L,n) 0.0] for n in 0:(iterations+extra_iterations) for sign in [-1,1]],by=x->x[1]))
                for (G_a,R_a) in G2R_a for r_a in 1:R_a
            )
        )

        add_correlation_contribution!(
            spectral_functions["spectral"],
            impurity_operators["particle"],
            impurity_operators["particle"],
            oindex2dimensions,
            irrEU,
            0 ,
            broadening_distribution ,
            spectral_broadening ,
            iterscale(scale,L,0) ,
            K_factor ; 
            correlation_type="spectral",
            T=spectral_temperature ,
            limit_shell = iterations==0 ,
            extra_iterations=extra_iterations
        )
        #compute_correlation_peaks(
        #    impurity_operators["particle"],
        #    impurity_operators["particle"],
        #    oindex2dimensions,
        #    irrEU,
        #    z,
        #    0 ;
        #    correlation_type="spectral",
        #    T=spectral_temperature,
        #    iteration_scale=iterscale(scale,L,0)
        #)

        #M = pcgred_atom 
        #part0 = get_partition0(irrEU,oindex2dimensions)
        #if orbitalresolved 
        #    A = redM2A_orbitalresolved( 
        #            M,
        #            collect(multiplets_a_atom),
        #            cg_o_fullmatint,
        #            cg_s_fullmatint,
        #            irrEU,
        #            part0
        #    )
        #else
        #    A = redM2A( 
        #            M,
        #            collect(multiplets_a_atom),
        #            cg_o_fullmatint,
        #            cg_s_fullmatint,
        #            irrEU,
        #            part0
        #    )
        #end
        #AA = [A]

        Mo_tot = length(oirreps2indices) 
        II_a = collect(Set([G[2] for G in get_irreps( multiplets_a_atom )]))
        Ms_atomspin = maximum([m[3] for m in multiplets_atom])
        Ms_shellspin = maximum([m[3] for m in multiplets_shell]) 
        Ms_tot = maximum((max_spin2,Ms_atomspin,Ms_shellspin))
        Ms_shell = maximum((Ms_atomspin,Ms_shellspin))
        Karray_orbital,Karray_spin = 
                    compute_Ksum_arrays(
                        oindex2dimensions,
                        cg_o_fullmatint,
                        cg_s_fullmatint,
                        Mo_tot ,
                        II_a ,
                        Ms_tot ,
                        Ms_shell)
    end
        
    #   ------------------   #
    #%% combinations μ,i→u %%#
    #   ------------------   #
    combinations_Gu_muiualpha = get_combinations_Gu_muiualpha_initial(
                identityrep ,
                calculation ,
                multiplets_atom ,
                oirreps2indices ;
                verbose=true )

    #   ---------------------------------------   #
    #%% matrix construction and diagonalization %%#
    #   ---------------------------------------   #
    (irrEU,combinations_Gu_muiualpha) = construct_and_diagonalize_uHv_totalangularmomentum( 
                    multiplets_atom , 
                    multiplets_shell ,
                    irrEU , 
                    hop_symparams_int , 
                    dsum ,
                    lehmann_iaj ,
                    lehmann_muanu ,
                    collect(multiplets_a_atom) , 
                    collect(multiplets_a_shell) ;
                    conduction_diagonals=channels_diagonals[1],
                    verbose=false );
    println( "-----------------------------------------------" )
    println( "SPECTRUM OF ATOM + INNERMOST SHELL (NORMALIZED)" )
    println()
    print_spectrum( irrEU )
    println()
    println( "-----------------------------------------------" )
    println()

    #   --------------------------- #
    #%% update impurity information # 
    #   --------------------------- #
    if compute_impmults
        mm_i,m_imp = update_impmultinfo( 
                        mm_i ,
                        irrEU ,
                        betabar ,
                        oindex2dimensions ,
                        combinations_Gu_muiualpha )
    end

    #   --------------------------- #
    #%% update spectral information # 
    #   --------------------------- #
    if spectral 

        impurity_operators["particle"] = update_operator( impurity_operators["particle"], 
                                                          collect(multiplets_a_atom) ,
                                                          Karray_orbital ,
                                                          Karray_spin ,
                                                          combinations_uprima ,
                                                          irrEU )
        add_correlation_contribution!(
            spectral_functions["spectral"],
            impurity_operators["particle"],
            impurity_operators["particle"],
            oindex2dimensions,
            irrEU,
            1 ,
            broadening_distribution ,
            spectral_broadening ,
            iterscale(scale,L,1) ,
            K_factor ; 
            correlation_type="spectral",
            T=spectral_temperature ,
            limit_shell = iterations==1 ,
            extra_iterations=extra_iterations
        )
        #compute_correlation_peaks(
        #    impurity_operators["particle"],
        #    impurity_operators["particle"],
        #    oindex2dimensions,
        #    irrEU,
        #    z,
        #    1 ;
        #    correlation_type="spectral",
        #    T=spectral_temperature,
        #    iteration_scale=iterscale(scale,L,1)
        #)

        #if orbitalresolved 
        #    M, AA = update_redmat_AA_CGsummethod_orbitalresolved(
        #            M,
        #            irrEU ,
        #            combinations_uprima ,
        #            collect(multiplets_a_atom) ,
        #            cg_o_fullmatint ,
        #            cg_s_fullmatint ,
        #            Karray_orbital ,
        #            Karray_spin ,
        #            AA ,
        #            oindex2dimensions ;
        #            verbose=false )
        #else
        #    M, AA = update_redmat_AA_CGsummethod(
        #            M,
        #            irrEU ,
        #            combinations_uprima ,
        #            collect(multiplets_a_atom) ,
        #            cg_o_fullmatint ,
        #            cg_s_fullmatint ,
        #            Karray_orbital ,
        #            Karray_spin ,
        #            AA ,
        #            oindex2dimensions ;
        #            verbose=false )
        #end
    end


    #   =============   #
    #%% NRG PROCEDURE %%# 
    #   =============   #
    println()
    println( ":::::::::::::::::::::" )
    println( "--- NRG PROCEDURE ---" )
    println( ":::::::::::::::::::::" )
    println()
    if !spectral

        nrg = NRG_totalangularmomentum( 
            label,
            calculation,
            iterations, 
            cutoff_type, 
            cutoff_magnitude,
            L,
            irrEU,
            multiplets_shell, 
            cg_o_comb2A,
            dsum,
            ksum,
            lehmann_muanu,
            collect(multiplets_a_shell) ,
            combinations_Gu_muiualpha,
            betabar,
            [1],
            channels_codiagonals,
            max_J2;
            mine=mine ,
            z=z ,
            compute_impmults=compute_impmults ,
            mult2index=mult2index ,
            orbital_multiplets=orbital_multiplets ,
            mm_i=mm_i ,
            channels_diagonals=channels_diagonals 
        )
        #nrg = NRG( label ,
        #           calculation ,
        #           iterations,
        #           cutoff_type,
        #           cutoff_magnitude,
        #           L,
        #           hop_symparams_int,
        #           irrEU,
        #           multiplets_shell,
        #           cg_o_fullmatint,
        #           cg_s_fullmatint,
        #           keys_as_dict_o ,
        #           keys_as_dict_s ,
        #           Csum_o_array ,
        #           Csum_s_array ,
        #           Bsum_o_array ,
        #           Bsum_s_array ,
        #           pcgred_shell,
        #           collect(multiplets_a_shell), 
        #           combinations_uprima,
        #           betabar,
        #           oindex2dimensions,
        #           channels_codiagonals ,
        #           max_spin2 ;
        #           mine=mine ,
        #           distributed=distributed ,
        #           z=z ,
        #           verbose=false ,
        #           precompute_iaj=precompute_iaj ,
        #           compute_impmults=compute_impmults ,
        #           mult2index=mult2index ,
        #           orbital_multiplets=orbital_multiplets ,
        #           mm_i=mm_i ,
        #           channels_diagonals=channels_diagonals )
    elseif dmnrg

        nrg = NRG( label ,
                   calculation ,
                   iterations,
                   cutoff_type,
                   cutoff_magnitude,
                   L,
                   hop_symparams_int,
                   copy(irrEU),
                   multiplets_shell,
                   cg_o_fullmatint,
                   cg_s_fullmatint,
                   keys_as_dict_o ,
                   keys_as_dict_s ,
                   Csum_o_array ,
                   Csum_s_array ,
                   Bsum_o_array ,
                   Bsum_s_array ,
                   pcgred_shell ,
                   collect(multiplets_a_shell), 
                   copy(combinations_uprima),
                   betabar,
                   oindex2dimensions,
                   channels_codiagonals ,
                   max_spin2 ;
                   mine=mine ,
                   distributed=distributed ,
                   z=z ,
                   dmnrg=dmnrg ,
                   dmnrg_run=1 ,
                   spectral_temperature=spectral_temperature ,
                   shell_dimension=shell_dimension ,
                   scale=Float64(scale) ,
                   precompute_iaj=precompute_iaj ,
                   compute_impmults=compute_impmults ,
                   mult2index=mult2index ,
                   orbital_multiplets=orbital_multiplets ,
                   mm_i=mm_i ,
                   channels_diagonals=channels_diagonals )

        NRG( label ,
             calculation ,
             iterations,
             cutoff_type,
             cutoff_magnitude,
             L,
             hop_symparams_int,
             irrEU,
             multiplets_shell,
             cg_o_fullmatint,
             cg_s_fullmatint,
             keys_as_dict_o ,
             keys_as_dict_s ,
             Csum_o_array ,
             Csum_s_array ,
             Bsum_o_array ,
             Bsum_s_array ,
             pcgred_shell ,
             collect(multiplets_a_shell), 
             combinations_uprima,
             betabar,
             oindex2dimensions,
             channels_codiagonals ,
             max_spin2 ;
             mine=mine ,
             distributed=distributed ,
             z=z ,
             dmnrg=dmnrg ,
             dmnrg_run=2 ,
             density_matrices=nrg.density_matrices ,
             spectral=true ,
             spectral_functions=spectral_functions ,
             spectral_broadening=spectral_broadening ,
             broadening_distribution=broadening_distribution ,
             K_factor=K_factor ,
             orbitalresolved=orbitalresolved ,
             impurity_operators=impurity_operators ,
             spectral_temperature=spectral_temperature ,
             #M=M,
             #AA=AA , 
             Karray_orbital=Karray_orbital ,
             Karray_spin=Karray_spin ,
             multiplets_atomhop=collect(multiplets_a_atom) ,
             scale=Float64(scale) ,
             precompute_iaj=precompute_iaj ,
             compute_impmults=compute_impmults ,
             mult2index=mult2index ,
             orbital_multiplets=orbital_multiplets ,
             mm_i=mm_i ,
             channels_diagonals=channels_diagonals ,
             half_weight_idx=nrg.half_weight_idx ,
             half_weight_energy=nrg.half_weight_energy )

    elseif spectral

        nrg = NRG( label ,
                   calculation ,
                   iterations,
                   cutoff_type,
                   cutoff_magnitude,
                   L,
                   hop_symparams_int,
                   irrEU,
                   multiplets_shell,
                   cg_o_fullmatint,
                   cg_s_fullmatint,
                   keys_as_dict_o ,
                   keys_as_dict_s ,
                   Csum_o_array ,
                   Csum_s_array ,
                   Bsum_o_array ,
                   Bsum_s_array ,
                   pcgred_shell ,
                   collect(multiplets_a_shell), 
                   combinations_uprima,
                   betabar,
                   oindex2dimensions,
                   channels_codiagonals ,
                   max_spin2 ;
                   mine=mine ,
                   distributed=distributed ,
                   z=z ,
                   verbose=false ,
                   spectral=true ,
                   spectral_functions=spectral_functions ,
                   spectral_broadening=spectral_broadening ,
                   broadening_distribution=broadening_distribution ,
                   K_factor=K_factor ,
                   orbitalresolved=orbitalresolved ,
                   impurity_operators=impurity_operators ,
                   spectral_temperature=spectral_temperature ,
                   extra_iterations=extra_iterations ,
                   #M=M,
                   #AA=AA , 
                   Karray_orbital=Karray_orbital ,
                   Karray_spin=Karray_spin ,
                   multiplets_atomhop=collect(multiplets_a_atom) ,
                   scale=Float64(scale) ,
                   precompute_iaj=precompute_iaj ,
                   compute_impmults=compute_impmults ,
                   mult2index=mult2index ,
                   orbital_multiplets=orbital_multiplets ,
                   mm_i=mm_i ,
                   channels_diagonals=channels_diagonals )
    end

    println()
    println( "END OF FULL NRG CALCULATION WITH z=$(z)" )
end
