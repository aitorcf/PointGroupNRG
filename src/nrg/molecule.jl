function nrg_molecule( 
            label::String ,
            calculation::String ,
            L::Float64 ,
            iterations::Int64 , cutoff_type::String ,
            cutoff_magnitude::R ,
            cg_o_dir::String ,
            multiplet_dir::String ,
            input_file::String ,
            chemical_potential::Float64 ,
            hop_symparams::Dict{ String , Matrix{ComplexF64} } ;
            z::Float64=0.0 ,
            max_spin2::Int64=10 ,
            channels_dos::Dict{ String , Vector{Function} }=Dict{ String , Vector{Function} }() ,
            discretization::String=discretization_default ,
            tridiagonalization::String=tridiagonalization_default ,
            enforce_particle_hole_symmetry::Bool=true,
            minmult::Int64=0 ,
            mine::Float64=0.0 ,
            betabar::Float64=1.0 ,
            spectral::Bool=false ,
            spectral_temperature::Float64=0.0 ,
            extra_iterations::Int64=0 ,
            operator::String="particle", # particle / spin_excitation
            K_factor::Float64=2.0 ,
            spectral_broadening::Float64=1.0 ,
            broadening_distribution::String="gaussian",
            orbitalresolved::Bool=false,
            only_j_diagonalization::Bool=false ,
            dmnrg::Bool=false ,
            compute_impmults::Bool=false ,
            band_width::Float64=1.0 ,
            scale_asymptotic::Bool=true ,
            compute_selfenergy::Bool=false ) where {R<:Real}

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
    #number_of_shell_orbitals = read_shell_input(input_file)
    number_of_shell_orbitals = sum( size(m,1) for (G,m) in hop_symparams )

    # read orbital rotation input
    orbital_rotation = read_rotation_input( input_file,
                                            number_of_impurity_orbitals )
    do_rotation = orbital_rotation!==diagm(map(x->ComplexF64(1.0),axes(orbital_rotation,1)))

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
    impurity_creation_operators_up_rotated = do_rotation ? map( 
        operator->transform_creation_operator(operator,impurity_creation_operators_up,orbital_rotation) ,
        impurity_creation_operators_up
    ) : impurity_creation_operators_up
    impurity_creation_operators_do_rotated = do_rotation ? map( 
        operator->transform_creation_operator(operator,impurity_creation_operators_do,orbital_rotation) ,
        impurity_creation_operators_do
    ) : impurity_creation_operators_do
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
    channels_tridiagonal = discretize_bands( 
        channels_dos ,
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
    scale::Float64 = band_width
    if scale_asymptotic
        codiagonals_first = channels_codiagonals[10][1][1]
        scale *= codiagonals_first
        channels_codiagonals = [
            # n (iterations) loop
            Dict(# I_o loop
                I_o => [# r_o loop
                    r_o_multiplet_codiagonals[n]/codiagonals_first
                    for (r_o,(_,r_o_multiplet_codiagonals)) in enumerate(I_o_multiplets_couplings)
                ]
                for (I_o,I_o_multiplets_couplings) in channels_tridiagonal_int
            )
            for n in 1:(iterations-1)
        ]
    end
    @show scale

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
            #irrEU_imp[G] = ([impurity_spectrum[G...,r]-E_0 for r in 1:rr],U)
            irrEU_imp[G] = ([impurity_spectrum[G...,r] for r in 1:rr],U)
        elseif N==N_0+1
            #irrEU_imp[G] = ([impurity_spectrum[G...,r]-E_0+chemical_potential for r in 1:rr],U)
            irrEU_imp[G] = ([impurity_spectrum[G...,r]+chemical_potential for r in 1:rr],U)
        elseif N==N_0-1
            #irrEU_imp[G] = ([impurity_spectrum[G...,r]-E_0-chemical_potential for r in 1:rr],U)
            irrEU_imp[G] = ([impurity_spectrum[G...,r]-chemical_potential for r in 1:rr],U)
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
    irrEU_notrescaled = copy(irrEU)
    irrEU = Dict{ClearIrrep,Tuple{Vector{Float64},Matrix{ComplexF64}}}(
        G=>( @.rescale(E,L,z,scale) , U ) for (G,(E,U)) in irrEU
    )

    println()
    println( "RESCALED IMPURITY HAMILTONIAN" )
    print_spectrum( irrEU )
    println( "HOPPING PARAMETERS" )
    print_dict(hop_symparams_int)


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
    print_dict(pcg_block)


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
    qq_a_shell = collect(filter(x->x[1]==1,keys(symstates_shell)))
    multiplets_a_shell = Set( (q[1:3]...,q[end]) for q in qq_a_shell )
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
    
    #   ------------------------------------------------------   #
    #%% shell hamiltonian for particle-hole asymmetric systems %%#
    #   ------------------------------------------------------   #
    #
    # occupation operator for the various shell orbitals
    #   { orbital_irrep_int => [{ G_mu=>[N_I_mu_r_mu] }] }
    orbitals2diagonalcounter = Dict{ Int64 , Vector{Dict{ IntIrrep , Vector{Float64} }} }()
    for (orbital_irrep,orbital_multiplets_diagonals) in channels_tridiagonal

        # orbital irrep in int format
        orbital_irrep_int = oirreps2indices[orbital_irrep]

        # occupations (irrEU format) for each orbital multiplet
        # belonging to the orbital_irrep
        orbital_multiplets_occupations = []

        # iterate through orbital multiplets
        for (r_o,_) in enumerate(orbital_multiplets_diagonals) 

            # one-electron orbital multiplet for which to compute the occupations
            orbital = (orbital_irrep,r_o)

            # operator for the chosen one-electron orbital
            counter_operator = electron_counter_sym( symstates_shell ,
                                                     orbital )

            # diagonalize operator
            orbital_count_diagonalization = irrEU2int(
                get_irrEU_initial( symstates_shell, counter_operator ),
                oirreps2indices
            )

            # introduce eigenvalues (number of particles) into dictionary
            push!( orbital_multiplets_occupations , 
                   Dict( G_mu=>N_mu 
                   for (G_mu,(N_mu,_)) in orbital_count_diagonalization) 
            )

        end

        # store multiplet occupations (irrEU format) in main dictionary
        orbitals2diagonalcounter[orbital_irrep_int] = copy(orbital_multiplets_occupations)

    end
    # occupations for the shell symstates
    #   { G_mu => [ r_mu -> [ I_o => [ r_o -> number_of_particles] ] ] }
    G2R_mu = Dict( 
        G_mu => length([m for m in multiplets_shell if m[1:3]==G_mu]) 
        for G_mu in Set( m[1:3] for m in multiplets_shell )
    )
    shell_sym_occupations = Dict{ IntIrrep , Vector{Dict{Int64,Vector{Float64}}} }(
        # (G_mu,R_mu) iteration
        G_mu => [# r_mu in 1:R_mu iteration
                    Dict(# (I_o,I_o_multiplets_occupations) iteration
                        I_o => [# r_o_occupations iteration
                            r_o_occupations[G_mu][r_mu]
                            for r_o_occupations in I_o_multiplets_occupations
                        ] 
                        for (I_o,I_o_multiplets_occupations) in orbitals2diagonalcounter
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
                            sum(# I_o loop
                                sum(# r_o loop
                                    r_o_couplings[n]*r_mu_multiplet_occupations[I_o][r_o] 
                                    for (r_o,(r_o_couplings,_)) in enumerate(I_o_multiplets_couplings)
                                )
                                for (I_o,I_o_multiplets_couplings) in channels_tridiagonal_int
                            )
                            for r_mu_multiplet_occupations in G_mu_multiplets_occupations
                        ]
                for (G_mu,G_mu_multiplets_occupations) in shell_sym_occupations
            )
            for n in 1:iterations
        ]
    if scale_asymptotic
        channels_diagonals = [
            # n (iteration) loop 
            Dict(# G_mu loop
                G_mu => [# r_mu loop
                            sum(# I_o loop
                                sum(# r_o loop
                                    r_o_couplings[n]*r_mu_multiplet_occupations[I_o][r_o]/codiagonals_first
                                    for (r_o,(r_o_couplings,_)) in enumerate(I_o_multiplets_couplings)
                                )
                                for (I_o,I_o_multiplets_couplings) in channels_tridiagonal_int
                            )
                            for r_mu_multiplet_occupations in G_mu_multiplets_occupations
                        ]
                for (G_mu,G_mu_multiplets_occupations) in shell_sym_occupations
            )
            for n in 1:iterations
        ]
    end

    #   ------------------------ #
    ##% conversion to int format #
    #   ------------------------ #
    multiplets_block = Set([ convert_to_int(m,oirreps2indices) 
                             for m in multiplets_block ])
    multiplets_imp = Set([ convert_to_int(m,oirreps2indices) 
                             for m in multiplets_imp ])
    irrEU_notrescaled = Dict( convert_to_int(G,oirreps2indices)=>x 
                             for (G,x) in irrEU_notrescaled )
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
                    collect(multiplets_block),
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

    # impurity (block) lehmann amplitudes for shell one-electrons
    multiplets_a_shell_int = Set( 
        convert_to_int(m,oirreps2indices) 
        for m in multiplets_a_shell
    )

    # all impurity (block) lehmann amplitudes (amplify for spectral later)
    multiplets_a_imp_int = multiplets_a_shell_int

    # compute impurity pcgred
    pcgred_block = get_redmat3(
                    pcg_block ,
                    multiplets_block ,
                    multiplets_a_shell_int ,
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

    impurity_operators = Dict{String,Dict{IntTripleG,Array{ComplexF64,3}}}()
    spectral_functions = Dict{String,Dict{IntMultiplet,Matrix{Float64}}}()
    if spectral
        print( "Computing excitation matrix... ") 

        atom_creops = symcreops_imp
        pcg_atomic_excitations = Dict()
        pcg_triple_excitations = Dict()
        if compute_selfenergy
            pcg_triple_excitations = Dict()
            for (qbra,sbra) in symstates_imp,
                (qket,sket) in symstates_imp,
                (qa,oa)     in atom_creops 

                qa_inversespin = (qa[1:4]...,-qa[5],qa[6])
                ca = oa
                ca_inversespin = atom_creops[qa_inversespin]
                aa_inversespin = adjoint(ca_inversespin)
                c = sbra*(ca_inversespin*aa_inversespin*ca)*sket 
                if !isapprox( abs2(c) , 0.0 )
                    pcg_triple_excitations[(qbra,qa,qket)] = c 
                end

            end
        end
        if operator=="particle"
            for (qbra,sbra) in symstates_imp,
                (qket,sket) in symstates_imp,
                (qa,oa)     in atom_creops 

                c = sbra*oa*sket 
                if !isapprox( abs2(c) , 0.0 )
                    pcg_atomic_excitations[(qbra,qa,qket)] = c 
                end

            end
        elseif operator=="spin_excitation"
            qtriplet = (2,"A1g",1.0,1,1.0,1)
            qa = (0,"A1g",1.0,1,1.0,1)
            qsinglet = (2,"A1g",0.0,1,0.0,1)
            pcg_atomic_excitations[(qtriplet,qa,qsinglet)] = 1.0
        elseif operator=="spin_excitation_complete" # not working
            qa_register = Set()
            for (qa1,oa1) in atom_creops,
                (qa2,oa2) in atom_creops

                # quantum numbers
                m1 = qa1[5]
                m2 = qa2[5]
                m1==m2 && continue
                r1 = qa1[end]
                r2 = qa2[end]
                r2>r1 && continue
                m = m1-m2
                co1 = oa1
                ao2 = adjoint(oa2)
                # same orbital, opposite spin, adjoint
                ao1_spininverted = adjoint(atom_creops[qa2[1:4]...,-qa2[5],qa2[end]])
                co2_spininverted = atom_creops[qa1[1:4]...,-qa1[5],qa2[end]]
                oa = co1*ao2+co2_spininverted*ao1_spininverted

                # outer multiplicity
                qa = (0,"A1g",1.0,1,m,1)
                if qa in qa_register
                    rmax = maximum([
                        qa_reg[end] for qa_reg in qa_register 
                        if qa_reg[1:5]==qa[1:5]
                    ])
                    qa = (qa[1:5]...,rmax+1)
                end
                push!( qa_register , qa )

                for (qbra,sbra) in symstates_imp,
                    (qket,sket) in symstates_imp

                    c = sbra*oa*sket 
                    if !isapprox( abs2(c) , 0.0 )
                        pcg_atomic_excitations[(qbra,qa,qket)] = c 
                    end

                end
            end
        end
        pcg_atomic_excitations = Dict( 
            (convert_to_int(k[1],oirreps2indices),
             convert_to_int(k[2],oirreps2indices),
             convert_to_int(k[3],oirreps2indices))=>v
             for (k,v) in pcg_atomic_excitations )
        pcg_triple_excitations = Dict( 
            (convert_to_int(k[1],oirreps2indices),
             convert_to_int(k[2],oirreps2indices),
             convert_to_int(k[3],oirreps2indices))=>v
             for (k,v) in pcg_triple_excitations )
        multiplets_a_imp = Set()
        if operator=="particle"
            multiplets_a_imp = Set(
                (convert_to_int(q,oirreps2indices)[1:3]...,q[6]) 
                for q in keys(atom_creops)
           )
        elseif operator=="spin_excitation"
            multiplets_a_imp = Set([ 
                (q_a[1:3]...,q_a[end]) for (qbra,q_a,qket) in keys(pcg_atomic_excitations)
            ])
        end

        impurity_operators["particle"] = get_redmat3( 
            pcg_atomic_excitations ,
            multiplets_imp ,
            multiplets_a_imp ,
            cg_o_fullmatint ,
            cg_s_fullmatint 
        )
        GG_a  = Set(G_a for (_,G_a,_) in keys(impurity_operators["particle"]))
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

        #if orbitalresolved 
        #    Mred, AA = setup_redmat_AA_orbitalresolved(
        #                pcg_atomic_excitations ,
        #                multiplets_block ,
        #                multiplets_a_imp ,
        #                cg_o_fullmatint ,
        #                cg_s_fullmatint ,
        #                irrEU ,
        #                oindex2dimensions ;
        #                verbose=false )
        #else 
        #    Mred, AA = setup_redmat_AA(
        #                pcg_atomic_excitations ,
        #                multiplets_block ,
        #                multiplets_a_imp ,
        #                cg_o_fullmatint ,
        #                cg_s_fullmatint ,
        #                irrEU ,
        #                oindex2dimensions ;
        #                verbose=false )
        #    if compute_selfenergy 
        #        Mred_se, AA_se = setup_selfenergy(
        #                    Mred ,
        #                    pcg_triple_excitations ,
        #                    multiplets_block ,
        #                    multiplets_a_imp ,
        #                    cg_o_fullmatint ,
        #                    cg_s_fullmatint ,
        #                    irrEU ,
        #                    oindex2dimensions ;
        #                    verbose=false )
        #    end
        #end
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

        if (operator=="particle" && only_j_diagonalization)
            @show irrEU
            ground_irrep = filter( 
                x->x[2][1][1]==E_0 ,
                collect(irrEU_notrescaled)
            )[1][1]
            ground_multiplet = (ground_irrep...,1)
            @show ground_multiplet
            irrEU_notrescaled_normalized = Dict( G=>(E.-E_0,U) for (G,(E,U)) in irrEU_notrescaled )
            J = compute_J_matrix( impurity_operators["particle"] , 
                                  irrEU_notrescaled_normalized , 
                                  ground_multiplet , 
                                  number_of_impurity_orbitals )
            jj,U = diagonalize_J( J )
            V2 = abs2.(diag(collect(hop_symparams)[1][2]))
            jjV2 = jj[1:number_of_shell_orbitals].*V2
            positivejjV2 = filter( x->x>0 , jjV2 )
            rhoJJ = jjV2./(2*band_width)
            rhoJ = (sum(filter(x->x>0,jjV2)))/(2*band_width)
            first_excited_energy = minimum(filter( 
                x->x>E_0 ,
                [E[1] for (E,U) in values(irrEU_notrescaled_normalized)]
            ))
            kondotemps = 0.4 * first_excited_energy .* sqrt.(abs.(0.5.*positivejjV2)) .* exp.(-ground_multiplet[3]./(2.0 .* positivejjV2))
            kondotemp = 0.4*
                        first_excited_energy*
                        sqrt(abs(rhoJ))*
                        exp(-1.0/rhoJ)
            println( "-----------------------------" )
            println( "Magnetic coupling constants:" )
            print( "j = ")
            println( jj )
            print( "ρJ = " )
            println( rhoJJ )
            print( "V2*j = ")
            println( jjV2 )
            print( "rhoJ = ")
            println( rhoJ )
            print( "estimated Kondo temperature: ")
            println( kondotemp )
            println( "separate Kondo temperatures: ")
            @show kondotemps
            println( "-----------------------------" )
            println( "Unitary matrix U" )
            for i in axes(U,1)

                for j in axes(U,2)
                    print( " " )
                    print( real(U[i,j]) )
                end
                println()

            end
            println( "-----------------------------" )
            println()
            println( "Input orbital rotation" )
            for i in axes(U,1)

                for j in 1:size(U,2)
                    print( " " )
                    print( real(orbital_rotation[i,j]) )
                end
                println()

            end
            println( "-----------------------------" )
            println()
            println( "Combined rotation" )
            c = orbital_rotation*U
            for i in axes(c,1)

                for j in axes(c,2)
                    print( " " )
                    print( real(c[i,j]) )
                end
                println()

            end
            println( "-----------------------------" )
            return
        end

        Mo_tot = length(oirreps2indices) 
        II_a = collect(Set([G[2] for G in get_irreps( multiplets_a_imp )]))
        Ms_atomspin = maximum([m[3] for m in multiplets_block])
        Ms_shellspin = maximum([m[3] for m in multiplets_shell]) 
        Ms_tot = Int64(maximum((max_spin2,Ms_atomspin,Ms_shellspin)))
        Ms_shell = Int64(maximum((Ms_atomspin,Ms_shellspin)))
        SS_a = collect(Set([G[3] for G in get_irreps( multiplets_a_imp )]))
        Karray_orbital,Karray_spin = 
                    compute_Ksum_arrays(
                        oindex2dimensions,
                        cg_o_fullmatint,
                        cg_s_fullmatint,
                        Mo_tot ,
                        II_a ,
                        Ms_tot ,
                        Ms_shell ;
                        SS_a=SS_a)

        println()
    end

    #   -------------------   #
    #%% clebsch-gordan sums %%#
    #   -------------------   #
    print( "Precomputing Clebsch-Gordan sums..." )
    @time begin
    Bsum_o_dict,Bsum_s_dict,Csum_o_dict,Csum_s_dict =
        precompute_CGsums(
                oirreps ,
                union(multiplets_a_imp_int,multiplets_a_shell_int) ,
                union(multiplets_block,multiplets_shell) ,
                max_spin2 ,
                oindex2dimensions ,
                cg_o_fullmatint ,
                cg_s_fullmatint )
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
    @time (irrEU,combinations_uprima) = matdiag_redmat(
                    multiplets_block , 
                    multiplets_shell ,
                    irrEU , 
                    hop_symparams_int , 
                    keys_as_dict_o ,
                    keys_as_dict_s ,
                    Csum_o_array ,
                    Csum_s_array ,
                    Bsum_o_array ,
                    Bsum_s_array ,
                    pcgred_block ,
                    pcgred_shell ,
                    collect(multiplets_a_imp_int) , 
                    collect(multiplets_a_shell_int) ,
                    combinations_uprima ;
                    conduction_diagonals=channels_diagonals[1],
                    verbose=false )
    println()
    println( "SPECTRUM AFTER ADDING FIRST SHELL" )
    print_spectrum(irrEU)

    # impurity thermodynamics 
    if (compute_impmults && calculation=="IMP")
        mm_i,m_imp = update_impmultinfo( 
                        mm_i ,
                        irrEU ,
                        betabar ,
                        oindex2dimensions ,
                        combinations_uprima )
    end

    if spectral
        impurity_operators["particle"] = update_operator( 
            impurity_operators["particle"], 
            collect(multiplets_a_imp) ,
            Karray_orbital ,
            Karray_spin ,
            combinations_uprima ,
            irrEU 
        )
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
        #if orbitalresolved
        #    Mred, AA = update_redmat_AA_CGsummethod_orbitalresolved(
        #                Mred ,
        #                irrEU ,
        #                combinations_uprima ,
        #                collect(multiplets_a_imp) ,
        #                cg_o_fullmatint ,
        #                cg_s_fullmatint ,
        #                Karray_orbital ,
        #                Karray_spin ,
        #                AA ,
        #                oindex2dimensions ;
        #                verbose=false )
        #    print_A( AA[end] )

        #else
        #    Mred, AA = update_redmat_AA_CGsummethod(
        #                Mred ,
        #                irrEU ,
        #                combinations_uprima ,
        #                collect(multiplets_a_imp) ,
        #                cg_o_fullmatint ,
        #                cg_s_fullmatint ,
        #                Karray_orbital ,
        #                Karray_spin ,
        #                AA ,
        #                oindex2dimensions ;
        #                verbose=false )
        #    if compute_selfenergy
        #        Mred_se, AA_se = update_selfenergy_CGsummethod(
        #                    Mred ,
        #                    Mred_se ,
        #                    irrEU ,
        #                    combinations_uprima ,
        #                    collect(multiplets_a_imp) ,
        #                    cg_o_fullmatint ,
        #                    cg_s_fullmatint ,
        #                    Karray_orbital ,
        #                    Karray_spin ,
        #                    AA_se ,
        #                    oindex2dimensions ;
        #                    verbose=false )

        #    end
        #end
    end

    #   =============== #
    #%% NRG CALCULATION #
    #   =============== #
    println()
    println( ":::::::::::::::::::::" )
    println( "--- NRG PROCEDURE ---" )
    println( ":::::::::::::::::::::" )
    println()
    if spectral
        if !compute_selfenergy
            nrg = NRG( label ,
                       calculation ,
                       iterations ,
                       cutoff_type ,
                       cutoff_magnitude ,
                       L ,
                       hop_symparams_int ,
                       irrEU ,
                       multiplets_shell ,
                       cg_o_fullmatint ,
                       cg_s_fullmatint ,
                       keys_as_dict_o ,
                       keys_as_dict_s ,
                       Csum_o_array ,
                       Csum_s_array ,
                       Bsum_o_array ,
                       Bsum_s_array ,
                       pcgred_shell ,
                       collect(multiplets_a_shell_int) , 
                       combinations_uprima ,
                       betabar ,
                       oindex2dimensions ,
                       channels_codiagonals ,
                       max_spin2 ;
                       mine=mine ,
                       z=z ,
                       spectral=true ,
                       spectral_functions=spectral_functions ,
                       impurity_operators=impurity_operators ,
                       spectral_temperature=spectral_temperature ,
                       spectral_broadening=spectral_broadening ,
                       broadening_distribution=broadening_distribution,
                       K_factor=K_factor ,
                       orbitalresolved=orbitalresolved ,
                       extra_iterations=extra_iterations ,
                       #M=Mred ,
                       #AA=AA , 
                       Karray_orbital=Karray_orbital ,
                       Karray_spin=Karray_spin ,
                       multiplets_atomhop=collect(multiplets_a_imp) ,
                       compute_impmults=compute_impmults,
                       mm_i=mm_i,
                       mult2index=mult2index,
                       orbital_multiplets=omults,
                       channels_diagonals=channels_diagonals ,
                       scale=scale )
        else
            nrg = NRG( label ,
                       calculation ,
                       iterations ,
                       cutoff_type ,
                       cutoff_magnitude ,
                       L ,
                       hop_symparams_int ,
                       irrEU ,
                       multiplets_shell ,
                       cg_o_fullmatint ,
                       cg_s_fullmatint ,
                       keys_as_dict_o ,
                       keys_as_dict_s ,
                       Csum_o_array ,
                       Csum_s_array ,
                       Bsum_o_array ,
                       Bsum_s_array ,
                       pcgred_shell ,
                       collect(multiplets_a_shell_int) , 
                       combinations_uprima ,
                       betabar ,
                       oindex2dimensions ,
                       channels_codiagonals ,
                       max_spin2 ;
                       mine=mine ,
                       z=z ,
                       spectral=true ,
                       spectral_broadening=etafac ,
                       broadening_distribution=broadening_distribution,
                       K_factor=K_factor ,
                       orbitalresolved=orbitalresolved,
                       M=Mred ,
                       AA=AA , 
                       Karray_orbital=Karray_orbital ,
                       Karray_spin=Karray_spin ,
                       multiplets_atomhop=collect(multiplets_a_imp) ,
                       compute_impmults=compute_impmults,
                       mm_i=mm_i,
                       mult2index=mult2index,
                       orbital_multiplets=omults,
                       channels_diagonals=channels_diagonals ,
                       scale=scale ,
                       compute_selfenergy=true ,
                       Mred_se=Mred_se ,
                       AA_se=AA_se )
        end
    else
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
                   pcgred_shell,
                   collect(multiplets_a_shell_int), 
                   combinations_uprima,
                   betabar,
                   oindex2dimensions,
                   channels_codiagonals,
                   max_spin2;
                   z=z ,
                   verbose=false ,
                   spectral=false ,
                   compute_impmults=compute_impmults,
                   mm_i=mm_i,
                   mult2index=mult2index,
                   orbital_multiplets=omults,
                   channels_diagonals=channels_diagonals,
                   scale=scale)
    end

    println()
    println( "END OF FULL NRG CALCULATION WITH z=$(z)" )
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
                              number_of_impurity_orbitals::Int64 )::Matrix{ComplexF64}


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

    # if there is no rotation input, return identity
    all(iszero.(rotation_matrix)) && return ComplexF64.(diagm(map(x->1,1:number_of_impurity_orbitals)))

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

# ------------------------------------ #
# J diagonalization and transformation #
# ------------------------------------ #

function read_J( filename::String )::Matrix{Float64}
    return readdlm( filename )
end

function diagonalize_J( J::Matrix{Float64} )::Tuple{Vector{Float64},Matrix{ComplexF64}}

    F = eigen( J , sortby=x->-x)
    jj = real.(F.values)
    U = F.vectors

    return jj,U

end

function transform_creation_operator( creation_operator::Operator , 
                                      all_creation_operators_spin::Vector{<:Operator} ,
                                      U::Matrix{ComplexF64} )::Operator

    # find the index of the operator in the operator vector
    operator_idx = findfirst( x->x==creation_operator , all_creation_operators_spin )

    # transformed operator vector
    transformed_operator_vector = conj.(U[operator_idx,:])

    # compute transformed operator
    transformed_operator = sum( transformed_operator_vector.*all_creation_operators_spin )
    println()

    return transformed_operator

end

function compute_spin_factor( S_c::Float64 , S_0::Float64 , N_c::Int64 , N_0::Int64 )

    if N_c<N_0

        S_c<S_0 && (return S_0)
        S_c>S_0 && (return -(S_0+1))

    elseif N_c>N_0

        S_c<S_0 && (return S_0+0.5)
        S_c>S_0 && (return -(S_0+0.5))

    end
    error( "error in spin factor calculation")
end

function compute_J_matrix( lehmann_reduced::Dict{NTuple{3,IntIrrep},Array{ComplexF64,3}} ,
                           irrEU_imp::Dict{IntIrrep,Tuple{Vector{Float64},Matrix{ComplexF64}}} ,
                           ground_multiplet::IntMultiplet ,
                           number_of_orbitals::Int64 )

    # initialize J matrix
    J = zeros(Float64,number_of_orbitals,number_of_orbitals)

    # ground irrep quantum numbers
    N_0 = ground_multiplet[1]
    S_0 = ground_multiplet[3]/2.0

    # iterate through lehmann amplitudes and impurity multiplets
    for ((G_1,G_a,G_2),mat) in lehmann_reduced,
        r_1 in axes(mat,1),
        r_2 in axes(mat,3)

        # multiplet quantum numbers
        M_1 = (G_1...,r_1)
        M_2 = (G_2...,r_2)

        # if no ground state, continue
        (M_1==ground_multiplet || M_2==ground_multiplet) || continue

        # select excited multiplet and gather related quantities
        M_c = M_1==ground_multiplet ? M_2 : M_1
        N_c = M_c[1]
        S_c = M_c[3]/2.0
        E_c = irrEU_imp[M_c[1:3]...][1][M_c[end]]

        # skip non-contributing excited multiplets
        abs(S_c-S_0)==0.5 || continue

        # compute spin factor
        spin_factor = compute_spin_factor( S_c , S_0 , N_c , N_0 )

        @show M_c, E_c, spin_factor

        # MO submatrix
        @views submat = mat[r_1,:,r_2]

        # tensor product
        for a1 in 1:number_of_orbitals,
            a2 in 1:number_of_orbitals

            J[a1,a2] += submat[a1]*submat[a2]/(spin_factor*E_c)

        end

    end

    return J

end
