function nrg_butterfly( 
            label::String ,
            calculation::String ,
            L::Float64 ,
            iterations::Int64 ,
            cutoff_type::String ,
            cutoff_magnitude::R ,
            cg_o_dir::String ,
            multiplet_dir::String ,
            input_file::String ,
            chemical_potential::Float64 ,
            hop_symparams::Dict{ String , Matrix{ComplexF64} } ;
            z::Float64=0.0 ,
            max_spin2::Int64=10 ,
            channel_dos::Dict{ String , Vector{Function} }=Dict{ String , Vector{Function} }() ,
            discretization="standard" ,
            minmult::Int64=0 ,
            mine::Float64=0.0 ,
            betabar::Float64=1.0 ,
            spectral::Bool=false ,
            K_factor::Float64=2.0 ,
            etafac::Float64=1.0 ,
            orbitalresolved::Bool=false,
            compute_impmults::Bool=false ) where {R<:Real}

    # strict defaults
    precompute_iaj::Bool=true
    
    # impmults only with imp
    compute_impmults = compute_impmults && (calculation=="IMP")

    println( "********************************" )
    println( "Full NRG calculation with z=$(z)" )
    println( "********************************" )
    println()

    # --------------------------- #
    # Clebsch-Gordan coefficients #
    # --------------------------- #
    
    print( "Obtaining Clebsch-Gordan coefficients... " )
    @time begin 
    (oirreps,
     oirreps2indices,
     oirreps2dimensions,
     oindex2dimensions,
     cg_o_fullmatint) = get_cg_o_info( cg_o_dir , ["A1g"]  )
    cg_s_fullmatint = get_cg_s_fullmatint( max_spin2 );
    end; println()
    
    # G1 x G2 = G1' + G2' + ...
    # (G1,G2) => [G1',G2',...]
    cg_o_fullmatint_keys = keys(cg_o_fullmatint)
    cg_s_fullmatint_keys = keys(cg_s_fullmatint)
    keys_as_dict_o::Dict{NTuple{2,Int64},Vector{Int64}} = Dict(
        (I1,I2)=>collect(Set(x[3] for x in cg_o_fullmatint_keys if x[1:2]==(I1,I2)))
        for (I1,I2,_) in cg_o_fullmatint_keys 
    )
    keys_as_dict_s::Dict{NTuple{2,Int64},Vector{Int64}} = Dict(
        (S1,S2)=>collect(Set(x[3] for x in cg_s_fullmatint_keys if x[1:2]==(S1,S2)))
        for (S1,S2,_) in cg_s_fullmatint_keys 
    )

    #   ==================================   #
    #&& ONE-SHELL SYMSTATES AND MULTIPLETS &&#
    #   ==================================   #

    # read impurity input
    number_of_impurity_orbitals,
    impurity_spectrum,
    impurity_symstates = read_impurity_input(input_file)

    # read shell input
    number_of_shell_orbitals = read_shell_input(input_file)

    # read orbital rotation input
    orbital_rotation = read_rotation_input( input_file,
                                            number_of_impurity_orbitals )

    #   --------------------- #
    #%% atomic symstate basis #
    #   --------------------- #

    # define dictionary
    identityrep = "A1g"
    hiztegia = merge(
        Dict{String,String}( 
            "o$i" => "A1g" for i in 1:number_of_impurity_orbitals
        ),
        Dict{String,Float64}(
            "u" =>  0.5,
            "d" => -0.5
        )
    )

    # inpurity states in simple format
    impurity_orbital_states = reduce(vcat , 
        [
        [(0,"o$i",1,"u"),(0,"o$i",1,"d")]
        for i in 1:number_of_impurity_orbitals
        ] 
    )

    # construct hilbert space and basis
    hilbert = HilbertSpace(impurity_orbital_states)
    basis_0 = CanonicalBasis( hilbert )

    ##   ------------------ #
    ##%% selected symstates #
    ##   ------------------ #

    # atomic symstates
    symstates_imp = Dict{Tuple{Int64,String,Float64,Int64,Float64,Int64},State}()

    # vacuum (seed) states
    vac = State( basis_0.states[1] , basis_0 )

    # creation operators 
    #
    # define eigenoperators c^\dagger_a wrt J
    impurity_creation_operators_up = [
        Operator( SymbolCreationOperator(o_s) , basis_0 )
        for o_s in filter( x->x[4]=="u" , impurity_orbital_states )
    ]
    impurity_creation_operators_do = [
        Operator( SymbolCreationOperator(o_s) , basis_0 )
        for o_s in filter( x->x[4]=="d" , impurity_orbital_states )
    ]
    impurity_creation_operators = reduce(vcat,[
        [u,d] for (u,d) in zip(impurity_creation_operators_up,impurity_creation_operators_do)
    ])
    # transformed operators c^\dagger_\alpha, old basis
    impurity_creation_operators_up_rotated = map( 
        operator->transform_creation_operator(operator,impurity_creation_operators_up,orbital_rotation) ,
        impurity_creation_operators_up
    )
    impurity_creation_operators_do_rotated = map( 
        operator->transform_creation_operator(operator,impurity_creation_operators_do,orbital_rotation) ,
        impurity_creation_operators_do
    )
    impurity_creation_operators_rotated = reduce(vcat,[
        [u,d] for (u,d) in zip(impurity_creation_operators_up_rotated,impurity_creation_operators_do_rotated)
    ])
    
    vac = State( basis_0.states[1] , basis_0 )
    symstates_imp = Dict{ClearQNums,State}()
    for (symQNums,slater_composition) in impurity_symstates

        # zero state
        symstates_imp[symQNums] = State(basis_0)

        for (coeff,slater_determinant) in slater_composition
            symstates_imp[symQNums] += coeff*slater_state(slater_determinant,impurity_creation_operators_rotated,vac)
        end
    end
    print_dict(symstates_imp)
    multiplets_imp = get_multiplets( symstates_imp )
    @show multiplets_imp


    println( "***********************************************************************")
    println( "IMPURITY MULTIPLETS" )
    print_multiplets_Nordered( multiplets_imp )
    #println( "ATOMIC SELECTED SYMSTATES" )
    #print_symstates_Nordered( symstates_atom )
    #println( "***********************************************************************")
    println()

    #   --------------   #
    #%% discretization %%#
    #   --------------   #

    # default behavior: eta=x->1/2
    if length(channels_dos)==0 
        channels_dos = Dict{String,Vector{Function}}( 
            orbital_irrep=>Function[x->0.5 for i in 1:size(hop_matrix,1)]
            for (orbital_irrep,hop_matrix) in hop_symparams 
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
        oirreps2indices[oirrep]=>v 
        for (oirrep,v) in channels_tridiagonal
    )
    # [ n -> { I_o => [ r_o -> diagonal_element ] } ]
    channels_codiagonals::Vector{Dict{Int64,Vector{Float64}}} = [
        # n (iterations) loop
        Dict(# I_o loop
            I_o => [# r_o loop
                r_o_multiplet_codiagonals[n]
                for (r_o,(_,r_o_multiplet_codiagonals)) in enumerate(I_o_multiplets_couplings)
            ]
            for (I_o,I_o_multiplets_couplings) in channels_tridiagonal_int
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
    scale::Float64 = 1.0 

    ##   ------------   #
    ##%% channel etas %%#
    ##   ------------   #
    #channel_etas = Dict{String,Vector{Function}}( k=>Function[x->0.5 for i in 1:size(v)[1]]
    #    for (k,v) in hop_symparams )

    ## symmetry-wise channel structure 
    ##       irrep => [ ϵ , ̄ϵ , ξ ]
    #channel_symstructure = Dict{String,Vector{Tuple{Vector{Float64},Vector{Float64},Vector{Float64},Float64}}}( 
    #        k=>[get_hoppings(iterations,L,z,coupling)
    #            for coupling in couplings]
    #        for (k,couplings) in channel_etas
    #)

    ## symmetry-wise etabar 
    #etabar_sym = Dict{String,Vector{Float64}}( 
    #        k=>[s[4] for s in v]
    #        for (k,v) in channel_symstructure 
    #)
    #@show etabar_sym
    #println()

    ## scale parameter: first asymptotic hopping element 
    ##       irrep => [ ̄ϵ[1] ]
    #scale_symparams = Dict{String,Vector{Float64}}( 
    #        k=>[h[2][1] for h in s] 
    #        for (k,s) in channel_symstructure 
    #)
    #@show scale_symparams
    #println()
    #scale::Float64 = maximum([v for (k,V) in scale_symparams for v in V])
    #factor_symparams = Dict{String,Vector{Float64}}( k=>v./scale for (k,v) in scale_symparams )
    #@show factor_symparams
    #println()

    ## hopping parameters: ξ = ϵ / ̄ϵ
    ##       irrep => [ ξ ]
    #xi_symparams::Dict{Int64,Vector{Vector{ComplexF64}}} = Dict( 
    #        oirreps2indices[k]=>[ComplexF64.(h[3].*factor_symparams[k][i]) for (i,h) in enumerate(s)] 
    #        for (k,s) in channel_symstructure
    #)
    #@show xi_symparams;

    ##   =================== #
    ##&& INITIAL HAMILTONIAN #
    ##   =================== #

    println( "Preparing initial Hamiltonian... " )
    @time begin

    #   ----- #
    #%% irreu #
    #   ----- #
    # clean spectrum
    irrEU_clean = get_irrEU_clean( "A1g" )
    # impurity spectrum formatting
    impurity_irreps_and_multiplicities = get_irreps( multiplets_imp , multiplicity=true )
    N_0 = minimum([ N for ((N,_,_),_) in impurity_irreps_and_multiplicities ])+1
    GG_neutral = filter( x->x[1][1]==N_0 , collect(impurity_irreps_and_multiplicities) )[1]
    E_0 = minimum(filter( x->(x[1][1:3] in GG_neutral) , collect(impurity_spectrum) ))[2]
    irrEU_imp = Dict{ClearIrrep,Tuple{Vector{Float64},Matrix{ComplexF64}}}( )
    for (G,rr) in impurity_irreps_and_multiplicities
        N,I,S = G
        U = ComplexF64.(diagm(map(x->1,1:rr)))
        if N==N_0
            irrEU_imp[G] = ([impurity_spectrum[G...,r]-E_0 for r in 1:rr],U)
        elseif N==N_0+1
            irrEU_imp[G] = ([impurity_spectrum[G...,r]-E_0+chemical_potential for r in 1:rr],U)
        elseif N==N_0-1
            irrEU_imp[G] = ([impurity_spectrum[G...,r]-E_0-chemical_potential for r in 1:rr],U)
        end
    end
    print_dict(irrEU_imp)
    irrEU = calculation=="IMP" ? irrEU_imp : irrEU_clean
    multiplets_block = calculation=="IMP" ? multiplets_imp :
                                            Set([ (0,"A1g",0.0,1) ])
    omults::Vector{ClearMultiplet} = ordered_multiplets(multiplets_block)
    mult2index::Dict{ClearMultiplet,Int64} = Dict( m=>i for (i,m) in enumerate(omults))

    end #timing
    println()

    print_spectrum( irrEU )

    #   ------------------- #
    #%% rescaled parameters #
    #   ------------------- #

    # orbital indices => irreps irreps
    oindices2irreps = Dict( v=>k for (k,v) in oirreps2indices )

    # correct hop scaling
    A_L = discretization=="yoshida1990" ? 0.5*(L+1)/(L-1)*log(L) : 1.0

    # rescale hoppings
    hop_symparams_int = Dict{Int64,Matrix{ComplexF64}}( 
        oirreps2indices[k]=>(@.rescale(v,L,z,scale/sqrt(A_L))) 
        for (k,v) in hop_symparams 
    )

    # rescale spectrum 
    irrEU = Dict{ClearIrrep,Tuple{Vector{Float64},Matrix{ComplexF64}}}(
        G=>( @.rescale(E,L,z,scale) , U ) for (G,(E,U)) in irrEU )
    )

    println( "RESCALED IMPURITY HAMILTONIAN" )
    print_spectrum( irrEU )
    print_dict(hop_symparams)
    return

    # ---------- #
    # atomic pcg #
    # ---------- #
    print( "Constructing atomic pcg... " )
    @time begin

    # symops
    symcreops_imp = Dict(
        (1,hiztegia[orbital],0.5,i,hiztegia[spinproj],parse(Int64,orbital[2])) => impurity_creation_operators[n]
        for (n,(_,orbital,i,spinproj)) in enumerate(impurity_orbital_states)
    )

    # hoppers
    pcg_block = Dict{NTuple{3,Tuple{Int64,String,Float64,Int64,Float64,Int64}},ComplexF64}()
    qq_a = sort( collect(keys(symcreops_imp)) , by=x->(x[end]*1000-x[3]) )[1:2*number_of_shell_orbitals]
    if calculation=="IMP"

        # pcg
        for (q1,s1) in symstates_imp,
            (q2,s2) in symstates_imp,
            qo  in qq_a

            lehmann = s1 * symcreops_imp[qo] * s2
            isapprox(lehmann,0.0) || (pcg_block[(q1,qo,q2)]=lehmann)

        end
    end

    end #timing
    println()


    #   ----------------------------------- #
    #%% shell symstates, multiplets and pcg #
    #   ----------------------------------- #

    print( "Obtaining shell info... " )

    # symstates
    symstates_shell, basis_shell, multiplets_shell_noint, multiplets_a_shell= get_symstates_basis_multiplets( 
        Dict{String,Int64}("A1g"=>number_of_shell_orbitals),
        oirreps2dimensions,
        identityrep,
        multiplet_dir,
        cg_o_dir ;
        verbose=true 
    )
    qq_a = collect(filter(x->x[1]==1,keys(symstates_shell)))
    multiplets_shell = multiplets2int(multiplets_shell_noint, oirreps2indices)
    hiztegia_shell = Dict{String,Any}(
        "A1g"=>"A1g", "u"=>0.5, "d"=>-0.5 
    )
    pcgred_shell = get_pcgred( 
        basis_shell ,
        symstates_shell::ClearSymstateDict ,
        multiplets_shell::IntMultipletSet ,
        hiztegia_shell ,
        oirreps2indices::Dict{String,Int64} ,
        cg_o_fullmatint::IntCG ,
        cg_s_fullmatint::IntCG 
    )
    print_dict(pcgred_shell)

    println()

    #   ------------------------ #
    ##% conversion to int format #
    #   ------------------------ #
    multiplets_block = Set([ convert_to_int(m,oirreps2indices) 
                             for m in multiplets_block ])
    multiplets_imp = Set([ convert_to_int(m,oirreps2indices) 
                             for m in multiplets_imp ])
    irrEU = Dict( convert_to_int(G,oirreps2indices)=>x 
                             for (G,x) in irrEU )
    qq_a = [ convert_to_int(q,oirreps2indices) 
                             for q in qq_a ]
    pcg_block = Dict( (convert_to_int(k[1],oirreps2indices),
                       convert_to_int(k[2],oirreps2indices),
                       convert_to_int(k[3],oirreps2indices))=>v 
                     for (k,v) in pcg_block )

    #   ---------------   #
    #%% combinations u' %%#
    #   ---------------   #
    print( "Computing combinations_uprima... " )
    @time combinations_uprima = get_combinations_uprima_initial(
                identityrep ,
                calculation ,
                multiplets_block ,
                oirreps2indices ;
                verbose=true )
    print_dict(combinations_uprima)

    #   ------------------- #
    #%% impurity multiplets #
    #   ------------------- #
    print( "Impurity multiplet setup... " )
    @time begin 
    mm_i,m_imp = setup_impmultinfo( 
                    multiplets_block ,
                    irrEU ,
                    betabar ,
                    oindex2dimensions )
    end #timing
    println()

    #   --------------------   #
    #%% reduced pcg matrices %%#
    #   --------------------   #
    print( "Computing reduced pcg... " )
    @time begin

    multiplets_a = Set( (q[1:3]...,q[end]) for q in qq_a )
    print_dict( pcg_block )
    @show multiplets_block
    @show multiplets_a
    pcgred_block = get_redmat3(
                    pcg_block ,
                    multiplets_block ,
                    multiplets_a ,
                    cg_o_fullmatint ,
                    cg_s_fullmatint ;
                    verbose=false )
    end
    println()
    println( "PCGRED BLOCK" )
    print_dict( pcgred_block )
    println( "PCGRED SHELL" )
    print_dict( pcgred_shell )
    println()




    #   -----------------   #
    #%% excitation matrix %%#
    #   -----------------   #

    if spectral
        print( "Computing excitation matrix... ") 

        atom_creops = symcreops_imp
        pcg_atomic_excitations = Dict()
        for (qbra,sbra) in symstates_imp,
            (qket,sket) in symstates_imp,
            (qa,oa)     in atom_creops 

            c = sbra*oa*sket 
            if !isapprox( abs2(c) , 0.0 )
                pcg_atomic_excitations[(qbra,qa,qket)] = c 
            end

        end
        pcg_atomic_excitations = Dict( 
            (convert_to_int(k[1],oirreps2indices),
             convert_to_int(k[2],oirreps2indices),
             convert_to_int(k[3],oirreps2indices))=>v
             for (k,v) in pcg_atomic_excitations )
        #multiplets_a = Set((q[1:3]...,q[6]) for (s1,q,s2) in keys(pcg_atomic_excitations))
        multiplets_atomhop = Set((convert_to_int(q,oirreps2indices)[1:3]...,q[6]) 
                           for q in keys(atom_creops))
        if orbitalresolved 
            Mred, AA = setup_redmat_AA_orbitalresolved(
                        pcg_atomic_excitations ,
                        multiplets_block ,
                        multiplets_atomhop ,
                        cg_o_fullmatint ,
                        cg_s_fullmatint ,
                        irrEU ,
                        oindex2dimensions ;
                        verbose=false )
        else 
            Mred, AA = setup_redmat_AA(
                        pcg_atomic_excitations ,
                        multiplets_block ,
                        multiplets_atomhop ,
                        cg_o_fullmatint ,
                        cg_s_fullmatint ,
                        irrEU ,
                        oindex2dimensions ;
                        verbose=false )
        end
        #G0 = (2,1,2)
        #GM = (1,1,1)
        #GM_map = Dict(
        #    ((1,1,1),1) => "M_0-",
        #    ((1,1,1),2) => "M_1-",
        #    ((4,1,2),1) => "M_0+",
        #    ((4,1,0),1) => "M_1+",

        #)
        #for ((G1,_,G2),matrix) in Mred
        #    G1==G0 && (GM = G2; idx=3)
        #    G2==G0 && (GM = G1; idx=1)
        #    for c in CartesianIndices(matrix)
        #        r = Tuple(c)[idx]
        #        a = Tuple(c)[2] 
        #        println( "$(GM_map[(GM),r]) , a = $a , lehmann = $(matrix[c])")
        #        println()
        #    end

        #end
        #println( "====" )
        #println( "Mred" )
        #println( "====" )
        #print_dict(Mred)
        #println()

        J = compute_J_matrix( Mred , irrEU , (2,1,2,1) , number_of_impurity_orbitals )
        jj,U = diagonalize_J( J )
        @show J
        @show jj
        @show U
        println("U")
        for i in 1:size(U,1)

            for j in 1:size(U,2)
                print( " " )
                print( real(U[i,j]) )
            end
            println()

        end
        println()
        println("orbital rotation")
        for i in 1:size(U,1)

            for j in 1:size(U,2)
                print( " " )
                print( real(orbital_rotation[i,j]) )
            end
            println()

        end
        return

        Mo_tot = length(oirreps2indices) 
        II_a = collect(Set([G[2] for G in get_irreps( multiplets_atomhop )]))
        Ms_atomspin = maximum([m[3] for m in multiplets_block])
        Ms_shellspin = maximum([m[3] for m in multiplets_shell]) 
        Ms_tot = Int64(maximum((max_spin2,Ms_atomspin,Ms_shellspin)))
        Ms_shell = Int64(maximum((Ms_atomspin,Ms_shellspin)))
        Karray_orbital,Karray_spin = 
                    compute_Ksum_arrays(
                        oindex2dimensions,
                        cg_o_fullmatint,
                        cg_s_fullmatint,
                        Mo_tot ,
                        II_a ,
                        Ms_tot ,
                        Ms_shell )

        println()
    end
    alpha = compute_ebar0_z( z , L ; discretization=discretization )


    #   -------------------   #
    #%% clebsch-gordan sums %%#
    #   -------------------   #
    print( "Precomputing Clebsch-Gordan sums..." )
    @time begin
    Bsum_o_dict,Bsum_s_dict,Csum_o_dict,Csum_s_dict =
        precompute_CGsums(
                oirreps ,
                multiplets_a,
                calculation=="IMP" ? multiplets_imp : multiplets_shell ,
                max_spin2 ,
                oindex2dimensions ,
                cg_o_fullmatint ,
                cg_s_fullmatint ;
                verbose=false )
    Bsum_o_array,Bsum_s_array,Csum_o_array,Csum_s_array = 
        CGsums_dict2array( Bsum_o_dict,
                           Bsum_s_dict,
                           Csum_o_dict,
                           Csum_s_dict ) 
    end #timing 
    println()


    #   --------------------------------------- #
    #%% diagonalization: atom + innermost shell #
    #   --------------------------------------- #
    print( "Diagonalizing atom + innermost shell... " )
    print_dict( combinations_uprima )
    @time (irrEU,combinations_uprima) = matdiag_redmat_old_discretization(
                    multiplets_block , 
                    multiplets_shell ,
                    irrEU , 
                    hop_symparams_int , 
                    cg_o_fullmatint , 
                    cg_s_fullmatint ,
                    keys_as_dict_o ,
                    keys_as_dict_s ,
                    Csum_o_array ,
                    Csum_s_array ,
                    Bsum_o_array ,
                    Bsum_s_array ,
                    pcgred_block ,
                    pcgred_shell ,
                    collect(multiplets_a) , 
                    collect(multiplets_a) , 
                    combinations_uprima , 
                    oindex2dimensions ;
                    verbose=false )
    println()
    println( "SPECTRUM AFTER ADDING FIRST SHELL" )
    print_spectrum(irrEU)

    # impurity thermodynamics 
    if compute_impmults
        mm_i,m_imp = update_impmultinfo( 
                        mm_i ,
                        irrEU ,
                        betabar ,
                        oindex2dimensions ,
                        combinations_uprima )
    end

    if spectral
        if orbitalresolved
            Mred, AA = update_redmat_AA_CGsummethod_orbitalresolved(
                        Mred ,
                        irrEU ,
                        combinations_uprima ,
                        collect(multiplets_atomhop) ,
                        cg_o_fullmatint ,
                        cg_s_fullmatint ,
                        Karray_orbital ,
                        Karray_spin ,
                        AA ,
                        oindex2dimensions ;
                        verbose=false )
            print_A( AA[end] )

        else
            Mred, AA = update_redmat_AA_CGsummethod(
                        Mred ,
                        irrEU ,
                        combinations_uprima ,
                        collect(multiplets_atomhop) ,
                        cg_o_fullmatint ,
                        cg_s_fullmatint ,
                        Karray_orbital ,
                        Karray_spin ,
                        AA ,
                        oindex2dimensions ;
                        verbose=false )
        end
    end

    #   =============== #
    #%% NRG CALCULATION #
    #   =============== #

    if spectral
        nrg = NRG_old_discretization( label ,
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
                   collect(multiplets_a) , 
                   combinations_uprima ,
                   betabar ,
                   oindex2dimensions ,
                   xi_symparams ;
                   z=z ,
                   discretization=discretization ,
                   spectral=true ,
                   K_factor=K_factor ,
                   orbitalresolved=orbitalresolved,
                   M=Mred ,
                   AA=AA , 
                   alpha=alpha,
                   spectral_broadening=etafac ,
                   Karray_orbital=Karray_orbital ,
                   Karray_spin=Karray_spin ,
                   multiplets_atomhop=collect(multiplets_atomhop) ,
                   compute_impmults=compute_impmults,
                   mult2index=mult2index,
                   orbital_multiplets=omults,
                   mm_i=mm_i)
    else
        nrg = NRG_old_discretization( label ,
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
                   pcgred_shell,
                   collect(multiplets_a), 
                   combinations_uprima,
                   betabar,
                   oindex2dimensions,
                   xi_symparams;
                   z=z ,
                   alpha=alpha,
                   discretization=discretization ,
                   verbose=false ,
                   spectral=false ,
                   compute_impmults=compute_impmults,
                   mult2index=mult2index,
                   orbital_multiplets=omults)
    end

end

# -------------------------- #
# File reading functionality #
# -------------------------- #

# read lines corresponding to dictionary entries
function read_dictionary_lines( split_lines::Vector{Vector{String}} , 
                                line_count::Int64 )

    parameter_dictionary = Dict()

    internal_line_count = line_count
    while true

        # increase line count
        internal_line_count += 1
        current_split_line = split_lines[internal_line_count]

        # stop if line is empty
        length(current_split_line)==0 && break

        # assign dictionary entry
        parameter_dictionary[current_split_line[1]] = current_split_line[3]

    end

    read_line_count = internal_line_count-line_count

    return parameter_dictionary,read_line_count

end

# read lines corresponding to impurity info
function read_impurity_input( file_name::String )

    multiplet_spectrum = Dict{ClearMultiplet,Float64}()
    symstates = Dict{ClearQNums,Vector{Tuple{ComplexF64,String}}}()
    number_of_orbitals = 0
    impurity_config = Dict{String,Int64}()

    reading_impurity_config = false
    impurity_found = false

    number_of_orbitals = 0
    N = 0
    I = ""
    S = 0.
    r = 0
    i = 0
    s = 0.
    E = 0.
    multipletQNums = (N,I,S,r)
    stateQNums = (N,I,S,i,s,r)

    # read file
    open( file_name , "r" ) do f

        split_lines = split.(strip.(readlines(f)))

        for (line_count,current_split_line) in enumerate(split_lines)

            # if we are not yet in the impurity section, continue
            if !impurity_found

                #empty string before impurity info -> continue
                length(current_split_line)==0 && continue
                current_split_line[1]=="IMPURITY" && (impurity_found=true)
                continue

            else 

                # empty string after impurity info -> break
                length(current_split_line)==0 && break
            end

            # skip comment
            current_split_line[1][1]=='#' && continue

            # number of orbitals
            if length(current_split_line)==1
                number_of_orbitals = parse(Int64,current_split_line[1])
            end

            ## impurity configuration
            #(reading_impurity_config && current_split_line[1]=="}") && (reading_impurity_config=false)
            #if reading_impurity_config
            #    impurity_config[current_split_line[1]] = parse(Int64,current_split_line[3]) 
            #end
            #current_split_line[1]=="impurity_config" && (reading_impurity_config = true; continue)

            # multiplet information line
            if length(current_split_line)==5

                N = parse(Int64,current_split_line[1])
                I = current_split_line[2]
                S = parse(Float64,current_split_line[3])
                r = parse(Int64,current_split_line[4])
                E = parse(Float64,current_split_line[5])

                multipletQNums = ( N , I , S , r )

                multiplet_spectrum[multipletQNums] = E

            end

            # state declaration
            if length(current_split_line)==2

                i = parse(Int64,current_split_line[1])
                s = parse(Float64,current_split_line[2])

                stateQNums = ( N , I , S , i , s , r )
                
                symstates[stateQNums] = []

            end

            # state coefficients in Slater basis
            if length(current_split_line)==3

                coefficient = parse(ComplexF64,current_split_line[1])
                slater_state = current_split_line[3]
                push!( symstates[stateQNums] , (coefficient,slater_state) )

            end

        end

    end

    return number_of_orbitals,multiplet_spectrum,symstates

end

# read lines corresponding to shell info
function read_shell_input( file_name::String )

    shell_found = false

    number_of_orbitals::Int64 = 0

    # read file
    open( file_name , "r" ) do f

        split_lines = split.(strip.(readlines(f)))

        for (line_count,current_split_line) in enumerate(split_lines)

            # if we are not yet in the impurity section, continue
            if !shell_found

                #empty string before impurity info -> continue
                length(current_split_line)==0 && continue
                current_split_line[1]=="SHELL" && (shell_found=true)
                continue

            else 

                # empty string after impurity info -> break
                length(current_split_line)==0 && break

            end

            # skip comment
            current_split_line[1][1]=='#' && continue

            # number of orbitals
            length(current_split_line)==1 && (number_of_orbitals=parse(Int64,current_split_line[1]))

        end

    end

    return number_of_orbitals

end

# read lines corresponding to rotation as column matrix
function read_rotation_input( file_name::String ,
                              number_of_impurity_orbitals::Int64 )


    rotation_found = false

    rotation_matrix = zeros(ComplexF64,number_of_impurity_orbitals,number_of_impurity_orbitals)

    # read file
    open( file_name , "r" ) do f

        split_lines = split.(strip.(readlines(f)))

        for (line_count,current_split_line) in enumerate(split_lines)

            # if we are not yet in the impurity section, continue
            if !rotation_found

                #empty string before impurity info -> continue
                length(current_split_line)==0 && continue
                current_split_line[1]=="ROTATION" && (rotation_found=true)
                continue

            else 

                # empty string after impurity info -> break
                length(current_split_line)==0 && break

            end

            # read rotation matrix
            for i in 1:number_of_impurity_orbitals,
                j in 1:number_of_impurity_orbitals

                rotation_matrix[i,j] = parse(ComplexF64,split_lines[line_count+i-1][j])

            end
            break
        end

    end

    return rotation_matrix

end

function slater_state( slater_string::String ,
                       impurity_creation_operators::Vector{<:Operator} ,
                       vacuum_state::State )

    # vector of creation operators for the slater state
    creation_operator_chain = []

    # translate slater string to chain of creation operators
    for (idx,occupation) in enumerate(slater_string)
        if occupation=='u'
            push!( creation_operator_chain , impurity_creation_operators[2*idx-1] ) 
        elseif occupation=='d'
            push!( creation_operator_chain , impurity_creation_operators[2*idx] ) 
        elseif occupation=='2'
            push!( creation_operator_chain , impurity_creation_operators[2*idx-1] , impurity_creation_operators[2*idx] ) 
        end
    end

    # no operators -> vacuum
    length(creation_operator_chain)==0 && (return vacuum_state)

    # return slater state
    return reduce( * , creation_operator_chain )*vacuum_state
    
end
