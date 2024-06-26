# get symmetry type
ispointspin(symmetry::String) = symmetry=="PS" || symmetry=="pointspin"
isdoublegroup(symmetry::String) = symmetry=="D" || symmetry=="doublegroup"
istotalangularmomentum(symmetry::String) = symmetry=="J" || symmetry=="totalangularmomentum" || symmetry=="S" || symmetry=="spin"
isorbital(symmetry::String) = ispointspin(symmetry) || isdoublegroup(symmetry)

# ======================= #
# TWO-PARTICLE MULTIPLETS #
# ======================= #
function multiplets_2particles_allsymmetries( 
            symmetry::String ,
            multiplets_dir::String ,
            impurity_config::Dict{SF,Int64} ;
            identityrep::String="",
            cg_o_dir::String="",
            max_SJ2::Int64=10) where {SF<:Union{String,Float64}}

    if isorbital(symmetry)
        if cg_o_dir==""
            error("""Clebsch-Gordan directory <cg_o_dir> required for orbital symmetries.
                     Add keyword argument:
                        cg_o_dir=<clebsch_gordan_directory>
            """)
        elseif identityrep==""
            error("""Identity irrep <identityrep> required for orbital symmetries.
                     Add keyword argument:
                        identityrep=<identity_representation>
            """)
        end
    end

    atom_orbital_irreps::Vector{SF} = collect(keys(impurity_config))

    # orbital symmetry
    oirreps::Vector{String},
    oirreps2indices::Dict{String,Int64},
    oirreps2dimensions::Dict{String,Int64},
    oindex2dimensions::Vector{Int64},
    cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}} = begin
        if isorbital(symmetry)

            oirreps,
            oirreps2indices,
            oirreps2dimensions,
            oindex2dimensions,
            cg_o_fullmatint = get_cg_o_info_nonsimple( 
                cg_o_dir , 
                atom_orbital_irreps 
            )

        else # total angular momentum

            identityrep::String = "A"
            oirreps = [identityrep]
            oirreps2indices = Dict(identityrep=>1)
            oirreps2dimensions = Dict(identityrep=>1)
            oindex2dimensions = [1]
            cg_o_fullmatint = Dict(
                (1,1,1) => ones(ComplexF64,1,1,1,1)
            )

        end
    end

    # spin symmetry
    cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} = begin
        if ispointspin(symmetry)

            get_cg_s_fullmatint(max_SJ2)

        elseif isdoublegroup(symmetry)

            Dict( (0,0,0)=>[1.0;;;] )

        end
    end

    #   ------------------------------- #
    #%% symstates, basis and multiplets #
    #   ------------------------------- #
    symstates_atom_noint,
    basis_atom,
    multiplets_atom_noint,
    multiplets_a_atom_noint = begin
        if ispointspin(symmetry)

            get_symstates_basis_multiplets_pointspin_nonsimple( 
                impurity_config,
                oirreps2dimensions,
                identityrep,
                multiplets_dir,
                cg_o_fullmatint ,
                cg_s_fullmatint ,
                oirreps2indices ;
                verbose=true 
            )

        elseif isdoublegroup(symmetry)

            get_symstates_basis_multiplets_doublegroups_nonsimple( 
                    impurity_config,
                    oirreps2dimensions,
                    identityrep,
                    multiplets_dir,
                    cg_o_fullmatint ,
                    cg_s_fullmatint ,
                    oirreps2indices ;
                    verbose=true )

        else # totalangularmomentum

            get_symstates_basis_multiplets_totalangularmomentum( 
                    impurity_config,
                    multiplets_dir;
                    verbose=true )

        end
    end

    multiplets_atom_twopart = filter(
        m->m[1]==2,
        multiplets_atom_noint
    )

    println( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" )
    println( "TWO-PARTICLE ATOMIC MULTIPLETS (with representative)" )
    println( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" )
    for ma in multiplets_atom_twopart 
        println( ma )
        println( symstates_atom_noint[(ma[1:3]...,1,ma[3],ma[4])] )
        println()
    end
end

# ================= #
# IMPURITY SPECTRUM #
# ================= #
function impurity_spectrum_allsymmetries(
            symmetry::String , # pointspin PS, double group D, total angular momentum J
            multiplets_dir::String ,
            impurity_config::Dict{SF,Int64} ,
            shell_config::Dict{SF,Int64} ,
            epsilon_symparams::Dict{ SF , Vector{ComplexF64} } ,
            u_symparams::Dict{ Tuple{String,Float64} , Matrix{ComplexF64} } ,
            hop_symparams::Dict{ SF , Matrix{ComplexF64} } ;
            cg_o_dir::String="",
            identityrep="",
            max_SJ2::Int64=10 ,
    ) where {SF<:Union{String,Float64}}

    println( "*****************************" )
    println( "Impurity spectrum calculation" )
    println( "*****************************" )
    println()

    if isorbital(symmetry)
        if cg_o_dir==""
            error("Clebsch-Gordan directory <cg_o_dir> required for orbital symmetries.")
        elseif identityrep==""
            error("Identity irrep <identityrep> required for orbital symmetries.")
        end
    end

    # orbital irreps present in the atom
    atom_orbital_irreps::Vector{SF} = collect(keys(impurity_config))

    println( "====================" )
    println( "SETUP AND PARAMETERS" )
    println( "====================" )
    println( "OCCUPATION ENERGIES" )
    print_dict( epsilon_symparams ) 
    println( "COULOMB PARAMETERS" )
    print_dict( u_symparams ) 
    println( "HYBRIDIZATION PARAMETERS" )
    print_dict( hop_symparams )
    println()


    # hiztegia
    hiztegia = Dict{String,Any}()
    if ispointspin(symmetry)

        hiztegia = Dict{String,Any}( o=>o for (o,_) in impurity_config )
        merge!( hiztegia , Dict( "u"=>0.5 , "d"=>-0.5 ) )

    elseif isdoublegroup(symmetry)

        hiztegia = Dict{String,Any}( o=>o for (o,_) in impurity_config )
        merge!( hiztegia , Dict( "-"=>0.0 ) )

    end # no hiztegia for total angular momentum symmetry

    #   ==========================   #
    #%% SYMMETRY-RELATED VARIABLES %%#
    #   ==========================   #

    # orbital symmetry
    oirreps::Vector{String} = []
    oirreps2indices::Dict{String,Int64} = Dict()
    oirreps2dimensions::Dict{String,Int64} = Dict()
    oindex2dimensions::Vector{Int64} = []
    cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}} = Dict()
    if isorbital(symmetry)

        (oirreps,
         oirreps2indices,
         oirreps2dimensions,
         oindex2dimensions,
         cg_o_fullmatint) = get_cg_o_info_nonsimple( cg_o_dir , atom_orbital_irreps )

    else # total angular momentum

        identityrep::String = "A"
        oirreps = [identityrep]
        oirreps2indices = Dict(identityrep=>1)
        oirreps2dimensions = Dict(identityrep=>1)
        oindex2dimensions = [1]
        cg_o_fullmatint = Dict(
            (1,1,1) => ones(ComplexF64,1,1,1,1)
        )

    end

    # for dmnrg
    shell_dimension::Int64 = 0
    if ispointspin(symmetry)

        shell_dimension = reduce( * , [4^(oirreps2dimensions[I]*R) for (I,R) in shell_config] )

    elseif isdoublegroup(symmetry)

        shell_dimension = reduce( * , [(oirreps2dimensions[I]*R) for (I,R) in shell_config] )

    elseif istotalangularmomentum(symmetry)

        shell_dimension = reduce( * , [jdim(J)*R for (J,R) in shell_config] )

    end

    # spin symmetry
    cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} = Dict{NTuple{3,Int64},Array{ComplexF64,3}}()
    if ispointspin(symmetry)

        cg_s_fullmatint = get_cg_s_fullmatint(max_SJ2)

    elseif isdoublegroup(symmetry)

        cg_s_fullmatint = Dict( (0,0,0)=>[1.0;;;] )

    end

    #   ===========   #
    #%% ATOMIC PART %%#
    #   ===========   #
    println()
    println( ":::::::::::::::::::" )
    println( "--- ATOMIC PART ---" )
    println( ":::::::::::::::::::" )
    println()

    #   ------------------------------- #
    #%% symstates, basis and multiplets #
    #   ------------------------------- #
    symstates_atom_noint, basis_atom, _, _ = begin
        if ispointspin(symmetry)

            get_symstates_basis_multiplets_pointspin_nonsimple( 
                impurity_config,
                oirreps2dimensions,
                identityrep,
                multiplets_dir,
                cg_o_fullmatint ,
                cg_s_fullmatint ,
                oirreps2indices ;
                verbose=true 
                )

        elseif isdoublegroup(symmetry)

            get_symstates_basis_multiplets_doublegroups_nonsimple( 
                    impurity_config,
                    oirreps2dimensions,
                    identityrep,
                    multiplets_dir,
                    cg_o_fullmatint ,
                    cg_s_fullmatint ,
                    oirreps2indices ;
                    verbose=true )

        else # totalangularmomentum

            get_symstates_basis_multiplets_totalangularmomentum( 
                    impurity_config,
                    multiplets_dir;
                    verbose=true )

        end
    end

    #   ----- #
    #%% model #
    #   ----- #
    # operators
    if isorbital(symmetry)

        epsilon = epsilon_sym( symstates_atom_noint , epsilon_symparams ; verbose=false )
        coulomb = u_sym( symstates_atom_noint , u_symparams ; verbose=false )

    else

        epsilon = epsilon_sym_totalangularmomentum( symstates_atom_noint , epsilon_symparams ; verbose=false )
        coulomb = u_sym_totalangularmomentum( symstates_atom_noint , u_symparams ; verbose=false )

    end

    # hamiltonian 
    H::Operator{typeof(basis_atom)} = epsilon + coulomb 

    #   --------   #
    #%% spectrum %%#
    #   --------   #
    irrEU_clear::Dict{ Tuple{Int64,String,Float64} , Tuple{Vector{Float64},Matrix{ComplexF64}} } = get_irrEU_initial(
        symstates_atom_noint,
        H;
        verbose=true
    )::Dict{ Tuple{Int64,String,Float64} , Tuple{Vector{Float64},Matrix{ComplexF64}} }
    println( "------------------------------------" )
    println( "IMPURITY SPECTRUM" )
    println()
    print_spectrum( irrEU_clear )
    println( "------------------------------------" )
    println()
end

# ==================== #
# FULL NRG CALCULATION #
# ==================== #

# impurity part
# - stardard anderson method
function construct_impurity(
        symmetry::String ,
        calculation::String ,
        impurity_config::Dict{SF,Int64} ,
        oirreps2dimensions::Dict{String,Int64} ,
        oirreps2indices::Dict{String,Int64} ,
        identityrep::String ,
        multiplets_dir::String ,
        cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}} ,
        cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} ,
        hiztegia ,
        epsilon_symparams::Dict{SF,Vector{ComplexF64}} ,
        u_symparams::Dict{Tuple{String,Float64},Matrix{ComplexF64}} ;
        until::String=""
    ) where {SF<:Union{String,Float64}}

    #   ------------------------------- #
    #%% symstates, basis and multiplets #
    #   ------------------------------- #
    symstates_atom_noint::Dict{Tuple{Int64,String,Float64,Int64,Float64,Int64},State} = Dict()
    multiplets_atom_noint::Set{Tuple{Int64,String,Float64,Int64}} = Set()
    multiplets_a_atom_noint::Set{Tuple{Int64,String,Float64,Int64}} = Set()
    mult2index::Dict{ClearMultiplet,Int64} = Dict()
    orbital_multiplets::Vector{ClearMultiplet} = []
    if calculation=="IMP"
        symstates_atom_noint,
        basis_atom,
        multiplets_atom_noint,
        multiplets_a_atom_noint = begin
            if ispointspin(symmetry)

                get_symstates_basis_multiplets_pointspin_nonsimple( 
                    impurity_config,
                    oirreps2dimensions,
                    identityrep,
                    multiplets_dir,
                    cg_o_fullmatint ,
                    cg_s_fullmatint ,
                    oirreps2indices ;
                    verbose=true 
                    )

            elseif isdoublegroup(symmetry)

                get_symstates_basis_multiplets_doublegroups_nonsimple( 
                        impurity_config,
                        oirreps2dimensions,
                        identityrep,
                        multiplets_dir,
                        cg_o_fullmatint ,
                        cg_s_fullmatint ,
                        oirreps2indices ;
                        verbose=true )

            else # totalangularmomentum

                get_symstates_basis_multiplets_totalangularmomentum(
                        impurity_config,
                        multiplets_dir;
                        verbose=true )

            end
        end

        if until=="2-particle multiplets"

            multiplets_atom_twopart = filter(
                m->m[1]==2,
                multiplets_atom_noint
            )

            println( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" )
            println( "TWO-PARTICLE ATOMIC MULTIPLETS (with representative)" )
            println( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" )
            for ma in multiplets_atom_twopart 
                println( ma )
                println( symstates_atom_noint[(ma[1:3]...,1,ma[3],ma[4])] )
                println()
            end

            return Dict{IntIrrep,Tuple{Vector{Float64},Matrix{ComplexF64}}}(
                (0,0,0)=>(zeros(Float64,0),zeros(ComplexF64,0,0))
            )
        end

        orbital_multiplets = ordered_multiplets(multiplets_atom_noint)
        mult2index = Dict( m=>i for (i,m) in enumerate(orbital_multiplets))
        multiplets_atom::Set{NTuple{4,Int64}} = multiplets2int( multiplets_atom_noint , 
                                                                oirreps2indices )
        multiplets_a_atom::Set{NTuple{4,Int64}} = multiplets2int( multiplets_a_atom_noint , 
                                                                  oirreps2indices )
    else 
        multiplets_atom_noint = Set([(0,identityrep,0.0,1)]) 
        multiplets_atom = multiplets2int(multiplets_atom_noint,oirreps2indices)
        multiplets_a_atom = multiplets_atom
    end

    #   -------------------------- #
    #%% reduced lehmann amplitudes #
    #   -------------------------- #
    lehmann_iaj::IntIrrepPCGNS = begin
        if calculation=="CLEAN"

            IntIrrepPCGNS()

        else
            if isorbital(symmetry)

                get_lehmann_reduced_orbitalsym( 
                    basis_atom ,
                    symstates_atom_noint::ClearSymstateDict ,
                    multiplets_atom::IntMultipletSet ,
                    hiztegia ,
                    oirreps2indices::Dict{String,Int64} ,
                    cg_o_fullmatint::IntCGNS ,
                    cg_s_fullmatint::IntCG ;
                    verbose=false
                )::IntIrrepPCGNS

            else

                get_lehmann_reduced_totalangularmomentum( 
                    basis_atom ,
                    symstates_atom_noint::ClearSymstateDict ,
                    multiplets_atom::IntMultipletSet ;
                    verbose=false 
                )::IntIrrepPCGNS

            end
        end
    end

    #   -------------------- #
    #%% impurity hamiltonian #
    #   -------------------- #
    if calculation=="IMP"

        # operators
        if isorbital(symmetry)

            epsilon = epsilon_sym( symstates_atom_noint , epsilon_symparams ; verbose=false )
            coulomb = u_sym( symstates_atom_noint , u_symparams ; verbose=false )

        else

            epsilon = epsilon_sym_totalangularmomentum( symstates_atom_noint , epsilon_symparams ; verbose=false )
            coulomb = u_sym_totalangularmomentum( symstates_atom_noint , u_symparams ; verbose=false )

        end

        # hamiltonian 
        H::Operator{typeof(basis_atom)} = epsilon + coulomb 

    end

    #   --------   #
    #%% spectrum %%#
    #   --------   #
    irrEU_clear::Dict{ Tuple{Int64,String,Float64} , Tuple{Vector{Float64},Matrix{ComplexF64}} } = begin
        if calculation=="IMP"

            get_irrEU_initial(symstates_atom_noint,H;verbose=true)::Dict{ Tuple{Int64,String,Float64} , Tuple{Vector{Float64},Matrix{ComplexF64}} }

        elseif calculation=="CLEAN" 

            get_irrEU_initial(identityrep,oirreps2indices)::Dict{ Tuple{Int64,String,Float64} , Tuple{Vector{Float64},Matrix{ComplexF64}} }

        end
    end
    println( "------------------------------------" )
    println( "IMPURITY SPECTRUM" )
    println()
    print_spectrum( irrEU_clear )
    println( "------------------------------------" )
    println()
    irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } = 
        irrEU2int( irrEU_clear , oirreps2indices ) 

    return (
        multiplets_atom,
        multiplets_a_atom,
        mult2index,
        orbital_multiplets,
        lehmann_iaj,
        irrEU
    )
end
# - ionic anderson method
function construct_impurity(
        spectrum::Dict{ClearIrrep,Vector{Float64}} ,
        lehmann_iaj::Dict{ClearTripleG,Array{ComplexF64,4}} ,
        calculation::String ,
        oirreps2indices::Dict{SF,Int64} ,
        identityrep::String ,
    ) where {SF<:Union{String,Float64}}

    #   ------------------------------- #
    #%% symstates, basis and multiplets #
    #   ------------------------------- #
    multiplets_atom_noint::Set{Tuple{Int64,String,Float64,Int64}} = Set()
    multiplets_a_atom_noint::Set{Tuple{Int64,String,Float64,Int64}} = Set()
    mult2index::Dict{ClearMultiplet,Int64} = Dict()
    orbital_multiplets::Vector{ClearMultiplet} = []
    if calculation=="IMP"
        multiplets_atom_noint = Set((G...,r) for (G,E) in spectrum for r in 1:length(E))
        multiplets_a_atom_noint = Set((G_a...,r_a) for ((_,G_a,_),arr) in lehmann_iaj for r_a in 1:size(arr,3))   
        orbital_multiplets = ordered_multiplets(multiplets_atom_noint)
        mult2index = Dict( m=>i for (i,m) in enumerate(orbital_multiplets))
        multiplets_atom::Set{NTuple{4,Int64}} = multiplets2int( multiplets_atom_noint , 
                                                                oirreps2indices )
        multiplets_a_atom::Set{NTuple{4,Int64}} = multiplets2int( multiplets_a_atom_noint , 
                                                                  oirreps2indices )
    else 
        multiplets_atom_noint = Set([(0,identityrep,0.0,1)]) 
        multiplets_atom = multiplets2int(multiplets_atom_noint,oirreps2indices)
        multiplets_a_atom = multiplets_atom
    end

    #   -------------------------- #
    #%% reduced lehmann amplitudes #
    #   -------------------------- #
    lehmann_iaj::IntIrrepPCGNS = begin
        if calculation=="CLEAN"

            IntIrrepPCGNS()

        else

            Dict(
                (convert_to_int(G_i,oirreps2indices),
                 convert_to_int(G_a,oirreps2indices),
                 convert_to_int(G_j,oirreps2indices)) => arr
                 for ((G_i,G_a,G_j),arr) in lehmann_iaj
            )

        end
    end

    #   --------   #
    #%% spectrum %%#
    #   --------   #
    irrEU_clear::Dict{ClearIrrep,Tuple{Vector{Float64},Matrix{ComplexF64}}} = begin
        if calculation=="IMP"

            Dict( (G,(E,Matrix{ComplexF64}(I,length(E),length(E)))) for (G,E) in spectrum )

        elseif calculation=="CLEAN" 

            get_irrEU_initial(identityrep,oirreps2indices)::Dict{ Tuple{Int64,String,Float64} , Tuple{Vector{Float64},Matrix{ComplexF64}} }

        end
    end
    println( "------------------------------------" )
    println( "IMPURITY SPECTRUM" )
    println()
    print_spectrum( irrEU_clear )
    println()
    println( "------------------------------------" )
    println()
    irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} } = 
        irrEU2int( irrEU_clear , oirreps2indices ) 

    return (
        multiplets_atom,
        multiplets_a_atom,
        mult2index,
        orbital_multiplets,
        lehmann_iaj,
        irrEU
    )
end
# input impurity
function construct_impurity(
        input_file::String
    )

    #   -------   #
    #%% reading %%#
    #   -------   #

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
    #
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
end

#function nrg_full_allsymmetries( 
function nrgfull(
            symmetry::String , # pointspin PS, double group D, total angular momentum J
            label::String ,
            L::Float64 ,
            iterations::Int64 ,
            cutoff::IF ,
            shell_config::Dict{SF,Int64} ,
            tunneling::Dict{ SF , Matrix{ComplexF64} } ;
            # directory containing the multiplets calculated by compute_multiplets
            multiplets_path::String="multiplets" ,
            # calculation with impurity ("IMP") or without ("CLEAN")
            calculation::String="IMP" ,
            # impurity input
            # - standard anderson
            impurity_config::Dict{SF,Int64}=Dict{SF,Int64}() ,
            onsite::Dict{SF,Vector{ComplexF64}}=Dict{SF,Vector{ComplexF64}}() ,
            interaction::Dict{Tuple{String,Float64},Matrix{ComplexF64}}=Dict{ Tuple{String,Float64} , Matrix{ComplexF64} }() ,
            # - ionic anderson
            spectrum::Dict{ClearIrrep,Vector{Float64}}=Dict{ClearIrrep,Vector{Float64}}() ,
            lehmann_iaj::Dict{ClearTripleG,Array{ComplexF64,4}}=Dict{ClearTripleG,Array{ComplexF64,4}}() ,
            # - ionic from input file
            input_file::String="" ,
            # run only until a certain point for partial information
            # - 2-particle multiplets
            # - impurity spectrum
            # - impurity-shell spectrum
            until::String = "",
            # other
            clebschgordan_path::String="" ,
            identityrep::String="" ,
            distributed::Bool=false ,
            z::Float64=0.0 ,
            max_SJ2::Int64=10 ,
            channels_dos::Dict{ SF , Vector{Function} }=Dict{ SF , Vector{Function} }() ,
            discretization::String=discretization_default ,
            tridiagonalization::String=tridiagonalization_default ,
            enforce_particle_hole_symmetry::Bool=true,
            mine::Float64=0.0 ,
            betabar::Float64=1.0 ,
            spectral::Bool=false ,
            spectral_broadening::Float64=0.5 ,
            broadening_distribution::String="loggaussian" ,
            K_factor::Float64=2.0 ,
            spectral_temperature::Float64=0.0 ,
            extra_iterations::Int64=0 ,
            dmnrg::Bool=false ,
            compute_impurity_projections::Bool=false ,
            scale_asymptotic::Bool=true ,
            band_width::Float64=1.0 ,
            print_spectrum_levels::Int64=0 ,
    ) where {SF<:Union{String,Float64},IF<:Union{Int64,Float64}}


    println( "********************************" )
    println( "Full NRG calculation with z=$(z)" )
    println( "********************************" )
    println()

    # impurity model
    model = begin
        if length(spectrum)!==0 && length(interaction)==length(impurity_config)==length(onsite)==0 && input_file==""
            "ionic"
        elseif length(impurity_config!==0) && length(spectrum)==0 && input_file==""
            "standard"
        elseif input_file!=="" && length(spectrum)==length(onsite)==length(interaction)==length(impurity_config)==0
            "input"
        else
            error(
                """Give either of these as input:
                     - spectrum and (optionally) lehmann amplitudes.
                     - impurity configuration, on-site energies and interaction.
                     - input file (experimental). """
            )
        end
    end

    # interaction format for total angular momentum
    if isorbital(symmetry)
        interaction = Dict(
            (k[1],0.0)=>v
            for (k,v) in interaction
        )
    else
        interaction = Dict(
            ("A",k[2])=>v
            for (k,v) in interaction
        )
    end

    # =================== #
    # ERRORS AND WARNINGS #
    # =================== #

    # cutoff_type from cutoff type
    cutoff_type = IF==Int64 ? "multiplet" : "energy"

    # spectral only
    if (spectral && calculation=="CLEAN") 
        println( "CLEAN calculation → spectral function will not be computed" )
    end

    # impurity projections only with IMP
    if compute_impurity_projections && calculation=="CLEAN"
        println("CLEAN calculation → impurity multiplet weights will not be coputed.")
    end
    compute_impurity_projections = compute_impurity_projections && (calculation=="IMP")

    if isorbital(symmetry)
        if clebschgordan_path==""
            error("Clebsch-Gordan directory <clebschgordan_path> required for orbital symmetries.")
        elseif identityrep==""
            error("Identity irrep <identityrep> required for orbital symmetries.")
        end
    end

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
    @show cutoff
    @show max_SJ2
    @show spectral
    println()
    if model=="standard"
        println("IMPURITY CONFIGURATION")
        print_dict(impurity_config)
        println("ON-SITE ENERGIES")
        print_dict(onsite) 
        println("INTERACTION PARAMETERS")
        print_dict(interaction) 
    elseif model=="ionic"
        println("SPECTRUM")
        print_dict(spectrum)
        println("LEHMANN AMPLITUDES")
        print_dict(lehmann_iaj)
    elseif model=="input"
        println("INPUT FILE: $input_file")
    end
    println("SHELL CONFIGURATION")
    print_dict(shell_config)
    println( "HYBRIDIZATION PARAMETERS" )
    print_dict( tunneling )
    println()

    #   ========   #
    #%% SYMMETRY %%#
    #   ========   #

    # orbital irreps present in the atom
    atom_orbital_irreps::Vector{SF} = begin 

        if model=="standard"

            collect(keys(impurity_config))

        elseif model=="ionic" # ionic

            imp_orbital_irreps::Vector{SF} = collect(Set(
                I for ((_,I,_),_) in spectrum
            ))
            shell_orbital_irreps::Vector{SF} = collect(keys(shell_config))
            collect(Set(vcat(imp_orbital_irreps,shell_orbital_irreps)))

        elseif model=="input"

            # just the identity
            collect(keys(shell_config)) 

        end
    end

    # hiztegia
    hiztegia = Dict{String,Any}()
    if ispointspin(symmetry)

        hiztegia = Dict{String,Any}( o=>o for (o,_) in impurity_config )
        merge!( hiztegia , Dict{String,Any}(o=>o for (o,_) in shell_config) )
        merge!( hiztegia , Dict( "u"=>0.5 , "d"=>-0.5 ) )

    elseif isdoublegroup(symmetry)

        hiztegia = Dict{String,Any}( o=>o for (o,_) in impurity_config )
        merge!( hiztegia , Dict{String,Any}(o=>o for (o,_) in shell_config) )
        merge!( hiztegia , Dict( "-"=>0.0 ) )

    end # no hiztegia for total angular momentum symmetry

    # orbital symmetry
    oirreps::Vector{String} = []
    oirreps2indices::Dict{String,Int64} = Dict()
    oirreps2dimensions::Dict{String,Int64} = Dict()
    oindex2dimensions::Vector{Int64} = []
    cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}} = Dict()
    if isorbital(symmetry)

        (oirreps,
         oirreps2indices,
         oirreps2dimensions,
         oindex2dimensions,
         cg_o_fullmatint) = get_cg_o_info_nonsimple( clebschgordan_path , atom_orbital_irreps )

    else # total angular momentum

        identityrep::String = "A"
        oirreps = [identityrep]
        oirreps2indices = Dict(identityrep=>1)
        oirreps2dimensions = Dict(identityrep=>1)
        oindex2dimensions = [1]
        cg_o_fullmatint = Dict(
            (1,1,1) => ones(ComplexF64,1,1,1,1)
        )

    end
    println("Orbital irreps and indices:")
    number_of_orbital_irreps = length(oirreps)
    for i in 1:number_of_orbital_irreps
        @printf "%5s " i
    end
    println()
    for i in 1:number_of_orbital_irreps
        @printf "%5s " oirreps[i]
    end
    println(); println()
    print("Performing sanity check on Clebsch-Gordan coefficients... ")
    clebsch_gordan_sanity_check(cg_o_fullmatint)
    println("Done\n\n")

    # for dmnrg
    shell_dimension::Int64 = 0
    if ispointspin(symmetry)

        shell_dimension = reduce( * , [4^(oirreps2dimensions[I]*R) for (I,R) in shell_config] )

    elseif isdoublegroup(symmetry)

        shell_dimension = reduce( * , [(oirreps2dimensions[I]*R) for (I,R) in shell_config] )

    elseif istotalangularmomentum(symmetry)

        shell_dimension = reduce( * , [jdim(J)*R for (J,R) in shell_config] )

    end

    # spin symmetry
    cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} = Dict{NTuple{3,Int64},Array{ComplexF64,3}}()
    if ispointspin(symmetry)

        cg_s_fullmatint = get_cg_s_fullmatint(max_SJ2)

    elseif isdoublegroup(symmetry)

        cg_s_fullmatint = Dict( (0,0,0)=>[1.0;;;] )

    end

    # G1 x G2 = n1×G1' + n2×G2' + ...
    # (G1,G2) => [G1',G2',...]
    cg_o_fullmatint_keys = keys(cg_o_fullmatint)
    cg_s_fullmatint_keys = keys(cg_s_fullmatint)
    keys_as_dict_o::Dict{NTuple{2,Int64},Vector{NTuple{2,Int64}}} = Dict(
        (I1,I2)=>collect(Set((I_3,size(arr,1)) for ((I_1,I_2,I_3),arr) in cg_o_fullmatint if (I_1,I_2)==(I1,I2)))
        for (I1,I2,_) in cg_o_fullmatint_keys 
    )
    cg_o_comb2A::Dict{ NTuple{3,Int64} , Int64 } = Dict(
        (I1,I2,I3)=>size(arr,1)
        for ((I1,I2,I3),arr) in cg_o_fullmatint
    )
    keys_as_dict_s::Dict{NTuple{2,Int64},Vector{Int64}} = Dict(
        (S1,S2)=>collect(Set(x[3] for x in cg_s_fullmatint_keys if x[1:2]==(S1,S2)))
        for (S1,S2,_) in cg_s_fullmatint_keys 
    )

    #   =============   #
    #%% IMPURITY PART %%#
    #   =============   #
    println()
    println( ":::::::::::::::::::::" )
    println( "--- IMPURITY PART ---" )
    println( ":::::::::::::::::::::" )
    println()
    multiplets_atom,
    multiplets_a_atom,
    mult2index,
    orbital_multiplets,
    lehmann_iaj,
    irrEU = begin

        # ionic model
        if model=="ionic"

            construct_impurity(
                spectrum::Dict{ClearIrrep,Vector{Float64}} ,
                lehmann_iaj::Dict{ClearTripleG,Array{ComplexF64,4}} ,
                calculation::String ,
                oirreps2indices::Dict{String,Int64} ,
                identityrep::String
            )


        # standard anderson model
        elseif model=="standard"

            ret = construct_impurity(
                symmetry::String ,
                calculation::String ,
                impurity_config::Dict{SF,Int64} ,
                oirreps2dimensions::Dict{String,Int64} ,
                oirreps2indices::Dict{String,Int64} ,
                identityrep::String ,
                multiplets_path::String ,
                cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}} ,
                cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} ,
                hiztegia ,
                onsite::Dict{SF,Vector{ComplexF64}} ,
                interaction::Dict{Tuple{String,Float64},Matrix{ComplexF64}} ;
                until=until
            )
            until=="2-particle multiplets" && return
            ret

        elseif model=="input"

        number_of_impurity_orbitals,
        impurity_spectrum,
        impurity_symstates = read_impurity_input(input_file)

        end
    end
    if until=="impurity spectrum" || until=="2-particle multiplets"
        return
    end


    #   ==============   #
    #%% DISCRETIZATION %%#
    #   ==============   #
    println()
    println( "::::::::::::::::::::::" )
    println( "--- DISCRETIZATION ---" )
    println( "::::::::::::::::::::::" )
    println()

    # default behavior: eta=x->1/2
    if length(channels_dos)==0 
        channels_dos = Dict{SF,Vector{Function}}( 
            orbital_irrep=>Function[x->0.5 for i in 1:size(hop_matrix,1)]
            for (orbital_irrep,hop_matrix) in tunneling 
        )
    end

    # channel coupling parameters
    print("Discretizing conduction bands... ")
    channels_tridiagonal = discretize_bands( 
        channels_dos ,
        L ,
        z , 
        iterations ;
        discretization=discretization ,
        tridiagonalization=tridiagonalization ,
        enforce_particle_hole_symmetry=enforce_particle_hole_symmetry 
    )
    println("Done.\n")
    channels_tridiagonal_int::Dict{Int64,Vector{Tuple{Vector{Float64},Vector{Float64}}}} = Dict()
    channels_codiagonals::Vector{Dict{Int64,Vector{Float64}}} = []
    if isorbital(symmetry)

        channels_tridiagonal_int = Dict(
            oirreps2indices[oirrep]=>v 
            for (oirrep,v) in channels_tridiagonal
        )

        # [ n -> { I_o => [ r_o -> diagonal_element ] } ]
        channels_codiagonals = [
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

    else # total angular momentum

        channels_tridiagonal_int = Dict(
            doublej(J)=>v 
            for (J,v) in channels_tridiagonal
        )

        # [ n -> { J2_o => [ J2_o -> diagonal_element ] } ]
        channels_codiagonals = [
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

    end

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
    # 
    if scale_asymptotic
        codiagonals_first = collect(values(channels_codiagonals[10]))[1][1]
        scale *= codiagonals_first
        channels_codiagonals = [
            # n (iterations) loop
            Dict(# IJ2_o loop
                IJ2_o => [# r_o loop
                    r_o_multiplet_codiagonals[n]/codiagonals_first
                    for (r_o,(_,r_o_multiplet_codiagonals)) in enumerate(IJ2_o_multiplets_couplings)
                ]
                for (IJ2_o,IJ2_o_multiplets_couplings) in channels_tridiagonal_int
            )
            for n in 1:(iterations-1)
        ]
    end
    println()
    println( "Scale: band-width × asymptotic hopping = $scale")
    println()

    # STILL NOT IMPLEMENTED
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

    #   -------------------- #
    #%% rescaled hamiltonian #
    #   -------------------- #
    println( "RESCALED HAMILTONIAN" )
    irrEU = Dict( k=>(@.rescale(E,L,z,scale),U) for (k,(E,U)) in irrEU )
    println("- Spectrum")
    print_spectrum(irrEU)
    println()
    oindices2irreps = Dict( v=>k for (k,v) in oirreps2indices )
    A_L = 0.5*(L+1)/(L-1)*log(L)  # correction factor
    hopscale = discretization=="yoshida1990" ? scale/sqrt(A_L) : scale
    tunneling_int::Dict{Int64,Matrix{ComplexF64}} = begin
        if isorbital(symmetry)

            Dict{Int64,Matrix{ComplexF64}}( oirreps2indices[k]=>(@.rescale(v,L,z,hopscale)) for (k,v) in tunneling )

        else

            Dict{Int64,Matrix{ComplexF64}}( doublej(J)=>(@.rescale(v,L,z,hopscale)) for (J,v) in tunneling )

        end
    end
    println("- hopping")
    print_dict(tunneling)
    println()


    #   ------------------------   #
    #%% impurity quantum numbers %%#
    #   ------------------------   #
    mm_i::Dict{IntMultiplet,Vector{Float64}} = Dict()
    m_imp::Vector{Float64} = []
    if compute_impurity_projections
        mm_i, m_imp = setup_impmultinfo(
            orbital_multiplets,
            irrEU,
            betabar,
            oindex2dimensions,
            oirreps2indices
        )
        println( "IMPURITY COMPOSITION" )
        println()
        @printf "%20s   %-10s\n" "multiplets" "weights" 
        for (multiplet,weight) in zip(orbital_multiplets,m_imp)
            @printf "%20s   %-.3f\n" multiplet weight
        end
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
    multiplets_a_shell_noint::Set{Tuple{Int64,String,Float64,Int64}} = begin
        if ispointspin(symmetry)

            get_symstates_basis_multiplets_pointspin_nonsimple(
                shell_config,
                oirreps2dimensions,
                identityrep,
                multiplets_path,
                cg_o_fullmatint ,
                cg_s_fullmatint ,
                oirreps2indices ;
                verbose=true 
            )

        elseif isdoublegroup(symmetry)

            get_symstates_basis_multiplets_doublegroups_nonsimple(
                shell_config,
                oirreps2dimensions,
                identityrep,
                multiplets_path,
                cg_o_fullmatint ,
                cg_s_fullmatint ,
                oirreps2indices ;
                verbose=true 
            )

        else # totalangularmomentum

            get_symstates_basis_multiplets_totalangularmomentum(
                shell_config,
                multiplets_path;
                verbose=true 
            )

        end
    end
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
    for (IJ,orbital_multiplets_diagonals) in channels_tridiagonal

        # orbital irrep in int format
        orbital_irrep_int = isorbital(symmetry) ? oirreps2indices[IJ] : doublej(IJ)

        # occupations (irrEU format) for each orbital multiplet
        # belonging to the orbital_irrep
        orbital_multiplets_occupations = []

        # iterate through orbital multiplets
        for (r_o,_) in enumerate(orbital_multiplets_diagonals) 

            # one-electron orbital multiplet for which to compute the occupations
            orbital = (IJ,r_o)

            # operator for the chosen one-electron orbital
            counter_operator = (
                if isorbital(symmetry)

                    electron_counter_sym( 
                        symstates_shell_noint ,
                        orbital 
                    )

                else

                    electron_counter_sym_totalangularmomentum( 
                        symstates_shell_noint ,
                        orbital 
                    )

                end
            )


            # diagonalize operator
            orbital_count_diagonalization = irrEU2int(get_irrEU_initial( symstates_shell_noint , counter_operator ),oirreps2indices)

            # introduce eigenvalues (number of particles) into dictionary
            push!( orbital_multiplets_occupations , Dict( G_mu=>N_mu for (G_mu,(N_mu,_)) in orbital_count_diagonalization) )

        end

        # store multiplet occupations (irrEU format) in main dictionary
        orbitals2diagonalcounter[orbital_irrep_int] = copy(orbital_multiplets_occupations)

    end
    # occupations for the shell symstates
    #   { G_mu => [ r_mu -> [ IJ2_o => [ r_o -> number_of_particles] ] ] }
    G2R_mu = Dict( 
        G_mu => length([m for m in multiplets_shell if m[1:3]==G_mu]) 
        for G_mu in Set( m[1:3] for m in multiplets_shell )
    )
    shell_sym_occupations = Dict{ IntIrrep , Vector{Dict{Int64,Vector{Float64}}} }(
        # (G_mu,R_mu) iteration
        G_mu => [# r_mu in 1:R_mu iteration
                    Dict(# (IJ2_o,IJ2_o_multiplets_occupations) iteration
                        IJ2_o => [# r_o_occupations iteration
                            r_o_occupations[G_mu][r_mu]
                            for r_o_occupations in IJ2_o_multiplets_occupations
                        ] 
                        for (IJ2_o,IJ2_o_multiplets_occupations) in orbitals2diagonalcounter
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
                        sum(# IJ2_o loop
                            sum(# r_o loop
                                r_o_couplings[n]*r_mu_multiplet_occupations[IJ2_o][r_o] 
                                for (r_o,(r_o_couplings,_)) in enumerate(IJ2_o_multiplets_couplings)
                            )
                            for (IJ2_o,IJ2_o_multiplets_couplings) in channels_tridiagonal_int
                        )
                        for r_mu_multiplet_occupations in G_mu_multiplets_occupations
                    ]
            for (G_mu,G_mu_multiplets_occupations) in shell_sym_occupations
        )
        for n in 1:iterations
    ]

    #   -------------------------- #
    #%% reduced lehmann amplitudes #
    #   -------------------------- #
    lehmann_muanu::Dict{IntTripleG,Array{ComplexF64,4}} = begin
        if isorbital(symmetry)

            get_lehmann_reduced_orbitalsym( 
                basis_shell,
                symstates_shell_noint,
                multiplets_shell,
                hiztegia,
                oirreps2indices,
                cg_o_fullmatint,
                cg_s_fullmatint;
                verbose=false
            )

        else

            get_lehmann_reduced_totalangularmomentum( 
                basis_shell,
                symstates_shell_noint,
                multiplets_shell;
                verbose=false
            )

        end
    end


    #   ====================================   #
    #%% COUPLING IMPURITY TO INNERMOST SHELL %%#
    #   ====================================   #
    println()
    println( "::::::::::::::::::::::::::::::::::::::::::::" )
    println( "--- COUPLING IMPURITY TO INNERMOST SHELL ---" )
    println( "::::::::::::::::::::::::::::::::::::::::::::" )
    println()

    #   ------------------------------ #
    #%% precompute clebsch-gordan sums #
    #   ------------------------------ #
    #
    print("Computing Clebsch-Gordan sums... ")
    dsum,ksum,fsum = begin
        if isorbital(symmetry)

            ClebschGordanSums(
                length(oindices2irreps),
                collect(union(Set(x[2] for x in multiplets_atom),Set(x[2] for x in multiplets_shell))),
                collect(Set(x[2] for x in multiplets_a_shell if x[1]==1)),
                cg_o_fullmatint,
                isdoublegroup(symmetry) ? 0 : max_SJ2,
                maximum(vcat([ x[3] for x in multiplets_atom ],[x[3] for x in multiplets_shell])),
                cg_s_fullmatint;
                compute_fsum=spectral,
                verbose=false 
            )

        else # totalangularmomentum

            ClebschGordanSums(
                    maximum([m[3] for m in vcat(collect(multiplets_a_shell),collect(multiplets_a_atom))]), # maximum_J2_onelectron
                    max_SJ2,
                    maximum([m[3] for m in vcat(collect(multiplets_shell),collect(multiplets_atom))]); # maximum_J2_onelectron
                    compute_fsum=spectral,
                    verbose=false 
                )

        end
    end
    println("Done.\n")

    #   -------- #
    #%% spectral #
    #   -------- #
    impurity_operators = Dict{String,Dict{IntTripleG,Array{ComplexF64,4}}}()
    spectral_functions = Dict{String,Dict{IntMultiplet,Matrix{Float64}}}()
    if spectral

        print("Setup for spectral function calculation... ")

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

        add_correlation_contribution_nonsimple!(
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

        println("Done.\n")
    end

    #   ------------------   #
    #%% combinations μ,i→u %%#
    #   ------------------   #
    combinations_Gu_muiualpha = get_combinations_Gu_muiualpha_initial(
        identityrep ,
        calculation ,
        multiplets_atom ,
        oirreps2indices ;
        verbose=true 
    )

    #   ---------------------------------------   #
    #%% matrix construction and diagonalization %%#
    #   ---------------------------------------   #
    print("Diagonalizing impurity + innermost shell... ")
    (irrEU,combinations_Gu_muiualpha) = begin
        if isorbital(symmetry)

            construct_and_diagonalize_uHv_orbitalsym( 
                multiplets_atom , 
                multiplets_shell ,
                irrEU , 
                tunneling_int , 
                keys_as_dict_o ,
                keys_as_dict_s ,
                dsum ,
                lehmann_iaj ,
                lehmann_muanu ,
                collect(multiplets_a_atom) , 
                collect(multiplets_a_shell) ;
                conduction_diagonals=channels_diagonals[1],
                verbose=false
            )

        else # totalangularmomentum

            construct_and_diagonalize_uHv_totalangularmomentum( 
                multiplets_atom , 
                multiplets_shell ,
                irrEU , 
                tunneling_int , 
                dsum ,
                lehmann_iaj ,
                lehmann_muanu ,
                collect(multiplets_a_atom) , 
                collect(multiplets_a_shell) ;
                conduction_diagonals=channels_diagonals[1],
                verbose=false 
            )

        end
    end
    println("Done.\n")

    println( "-----------------------------------------------" )
    println( "SPECTRUM OF ATOM + INNERMOST SHELL (NORMALIZED)" )
    println()
    print_spectrum( irrEU )
    println()
    println( "-----------------------------------------------" )
    println()
    if until=="impurity-shell spectrum"
        return
    end

    #   --------------------------- #
    #%% update impurity information # 
    #   --------------------------- #
    if compute_impurity_projections
        mm_i,m_imp = update_impmultinfo_allsymmetries(
            mm_i ,
            irrEU ,
            betabar ,
            oindex2dimensions ,
            combinations_Gu_muiualpha )
    end

    #   --------------------------- #
    #%% update spectral information # 
    #   --------------------------- #
    print("Updating spectral function variables... ")
    if spectral 

        impurity_operators["particle"] = update_operator_nonsimple( 
            impurity_operators["particle"],
            collect(multiplets_a_atom) ,
            fsum ,
            combinations_Gu_muiualpha ,
            irrEU ,
            cg_o_comb2A
        )
        add_correlation_contribution_nonsimple!(
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

    end
    println("Done.")


    #   =============   #
    #%% NRG PROCEDURE %%# 
    #   =============   #
    println()
    println( ":::::::::::::::::::::" )
    println( "--- NRG PROCEDURE ---" )
    println( ":::::::::::::::::::::" )
    println()
    if !spectral

        nrg = NRG_allsymmetries(
            symmetry,
            label,
            calculation,
            iterations,
            cutoff_type,
            cutoff,
            L,
            irrEU,
            multiplets_shell,
            keys_as_dict_o,
            keys_as_dict_s,
            cg_o_comb2A,
            dsum,
            ksum,
            lehmann_muanu,
            collect(multiplets_a_shell),
            combinations_Gu_muiualpha,
            betabar,
            oindex2dimensions,
            channels_codiagonals,
            max_SJ2;
            mine=mine,
            z=z,
            compute_impmults=compute_impurity_projections,
            mult2index=mult2index,
            orbital_multiplets=orbital_multiplets,
            mm_i=mm_i,
            channels_diagonals=channels_diagonals,
            scale=scale,
            print_spectrum_levels=print_spectrum_levels
        )

    elseif dmnrg # not updated

        nrg = NRG( label ,
                   calculation ,
                   iterations,
                   cutoff_type,
                   cutoff,
                   L,
                   tunneling_int,
                   copy(irrEU),
                   multiplets_shell,
                   cg_o_fullmatint,
                   cg_s_fullmatint,
                   keys_as_dict_o , keys_as_dict_s ,
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
             cutoff,
             L,
             tunneling_int,
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

        nrg = NRG_allsymmetries(
            symmetry,
            label,
            calculation,
            iterations,
            cutoff_type,
            cutoff,
            L,
            irrEU,
            multiplets_shell,
            keys_as_dict_o,
            keys_as_dict_s,
            cg_o_comb2A,
            dsum,
            ksum,
            lehmann_muanu,
            collect(multiplets_a_shell),
            combinations_Gu_muiualpha,
            betabar,
            oindex2dimensions,
            channels_codiagonals,
            max_SJ2;
            mine=mine,
            z=z,
            compute_impmults=compute_impurity_projections,
            mult2index=mult2index,
            orbital_multiplets=orbital_multiplets,
            mm_i=mm_i,
            channels_diagonals=channels_diagonals,
            spectral=true,
            spectral_functions=spectral_functions,
            spectral_broadening=spectral_broadening,
            broadening_distribution=broadening_distribution,
            K_factor=K_factor,
            impurity_operators=impurity_operators,
            spectral_temperature=spectral_temperature,
            extra_iterations=extra_iterations,
            multiplets_atomhop=collect(multiplets_a_atom) ,
            fsum=fsum,
            scale=scale,
            print_spectrum_levels=print_spectrum_levels
        )

    end

    println()
    println( "END OF FULL NRG CALCULATION WITH z=$(z)" )
end

# ========= #
# NRG STEPS #
# ========= #

function NRG_allsymmetries(
              symmetry::String ,
              label::String ,
              calculation::String ,
              iterations::Int64, 
              cutoff_type::String, 
              cutoff_magnitude::Number,
              L::Float64,
              irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
              multiplets_shell::Set{NTuple{4,Int64}}, 
              keys_as_dict_o::Dict{NTuple{2,Int64},Vector{NTuple{2,Int64}}} ,
              keys_as_dict_s::Dict{NTuple{2,Int64},Vector{Int64}} ,
              cg_o_comb2A::Dict{ NTuple{3,Int64} , Int64 } ,
              dsum::DSum ,
              ksum::KSum ,
              lehmann_muanu::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,4} },
              multiplets_a::Vector{NTuple{4,Int64}} , 
              combinations_Gu_muiualpha::Dict{NTuple{3,Int64}, Vector{Tuple{IntMultiplet,IntMultiplet,IntMultiplet,Int64}}},
              betabar::Float64 ,
              oindex2dimensions::Vector{Int64} ,
              channels_codiagonals::Vector{Dict{Int64,Vector{Float64}}},
              max_SJ2::Int64;
              verbose::Bool=false ,
              minmult::Int64=0 , 
              mine::Float64=0.0 ,
              z::Float64=0.0 ,
              spectral::Bool=false ,
              spectral_functions::Dict{String,Dict{IntMultiplet,Matrix{Float64}}}=Dict{String,Dict{IntMultiplet,Matrix{Float64}}}() ,
              K_factor::Float64=2.0 ,
              spectral_broadening::Float64=1.0 ,
              broadening_distribution::String="loggaussian",
              dmnrg::Bool=false ,
              dmnrg_run::Int64=1 ,
              shell_dimension::Int64=0 ,
              density_matrices::Vector{Dict{IntIrrep,Matrix{Float64}}}=Dict{IntIrrep,Matrix{Float64}}[] ,
              half_weight_idx::Int64=0 ,
              half_weight_energy::Float64=0.0 ,
              fsum::FSum=EmptyClebschGordanSum(FSum) ,
              multiplets_atomhop::Vector{NTuple{4,Int64}}=NTuple{4,Int64}[] ,
              scale::Float64=1.0 ,
              compute_impmults::Bool=false ,
              mult2index::Dict{ClearMultiplet,Int64}=Dict{ClearMultiplet,Int64}() ,
              orbital_multiplets::Vector{ClearMultiplet}=ClearMultiplet[] ,
              mm_i::Dict{NTuple{4,Int64},Vector{Float64}}=Dict{NTuple{4,Int64},Vector{Float64}}() ,
              write_spectrum::Bool=false ,
              channels_diagonals::Vector{Dict{IntIrrep,Vector{Float64}}}=[] ,
              impurity_operators::Dict{String,Dict{IntTripleG,Array{ComplexF64,4}}}=Dict{String,Dict{IntTripleG,Array{ComplexF64,4}}}() ,
              spectral_temperature::Float64=0.0 ,
              extra_iterations::Int64=0 ,
              compute_selfenergy::Bool=false ,
              print_spectrum_levels::Int64=0
    )

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
    spectral && (isdir("spectral") || mkdir("spectral"))

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
        (multiplets_block, discarded) = cut_off_nonsimple!( 
            irrEU;
            type=cutoff_type,
            cutoff=cutoff_magnitude,
            safeguard=true,
            minmult=minmult,
            mine=mine,
            verbose=false,
            impurity_operator=!iszero(length(impurity_operators)) ? impurity_operators["particle"] : Dict{IntTripleG,Array{ComplexF64,4}}()
        )

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
            IJ2_o => diagm(ComplexF64.(channels_codiagonals[n-1][IJ2_o])) 
            for IJ2_o in keys(channels_codiagonals[1])
        )
        println( "CONDUCTION COUPLING TERMS" )
        println( "Codiagonals (hoppings)")
        @printf "  %-8s  %-10s  %-s\n" "orbital" "multiplet" "amplitude"
        for (IJ2,hop_matrix) in hop_symparams,
            cartesian_index in CartesianIndices(hop_matrix)

            i,j = Tuple(cartesian_index)
            isapprox(hop_matrix[i,j],0.0) && continue

            @printf "  %-8s  %-3i => %-3i  %-.3f\n" IJ2 i j hop_matrix[i,j]

        end
        println( "Diagonals (occupation energies)")
        @printf "  %-8s  %-10s  %-s\n" "irrep" "multiplet" "amplitude"
        for (IJ2_o,IJ2_o_multiplets_diagonals) in channels_diagonals[n],
            (r_o,r_o_multiplet_diagonal) in enumerate(IJ2_o_multiplets_diagonals)

            @printf "  %-8s  %-3i => %-.3f\n" IJ2_o r_o r_o_multiplet_diagonal

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
        diagonalization_performance = @timed (irrEU,combinations_Gu_muiualpha) = (
            if isorbital(symmetry)

                construct_and_diagonalize_uHv_orbitalsym(
                    multiplets_block,
                    multiplets_shell,
                    irrEU,
                    hop_symparams,
                    keys_as_dict_o,
                    keys_as_dict_s,
                    dsum,
                    lehmann_iaj,
                    lehmann_muanu,
                    multiplets_a,
                    multiplets_a;
                    conduction_diagonals=channels_diagonals[1]
                )

            else # totalangularmomentum

                construct_and_diagonalize_uHv_totalangularmomentum( 
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

            end
        )

        if !iszero(print_spectrum_levels)
            print_spectrum(irrEU;multiplet_cutoff=print_spectrum_levels)
        end

        # information
        maximum_SJ2 = maximum(collect( G[3] for (G,(E,U)) in irrEU ))
        if maximum_SJ2>max_SJ2
            error( """
                ERROR: max_SJ2 too low!
                max_SJ2 = $max_SJ2
                maximum 2S (2J) found in calculation = $maximum_SJ2
            """)
        end
        maximum_spin2_all_iterations = maximum([ maximum_SJ2 , maximum_spin2_all_iterations ])
        maximum_eigenenergy_after_diagonalization = maximum([ maximum(E) for (_,(E,_)) in irrEU ])
        println( "DIAGONALIZATION" )
        println( "  time: $(diagonalization_performance.time)s" )
        println( "  memory: $(diagonalization_performance.bytes*10^-6)Mb" )
        println( "  garbage collection time: $(diagonalization_performance.gctime)s" )
        println( "  maximum 2S (2J) needed: $maximum_SJ2")
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
            @show spectrum_even
            (n%2==0)  && save_spectrum!( spectrum_even , irrEU )
            (n%2!==0) && save_spectrum!( spectrum_odd  , irrEU )
        end

        # impurity info
        if compute_impmults
            mm_i,m_imp = update_impmultinfo_allsymmetries(
                mm_i ,
                irrEU ,
                betabar ,
                oindex2dimensions ,
                combinations_Gu_muiualpha 
            )
            push!(impmults,m_imp)
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
                impurity_operators["particle"] = update_operator_nonsimple( 
                    impurity_operators["particle"], 
                    collect(multiplets_atomhop) ,
                    fsum ,
                    combinations_Gu_muiualpha ,
                    irrEU ,
                    cg_o_comb2A
                )
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
                add_correlation_contribution_nonsimple!(
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
            if any(any(isnan.(v)) for v in values(impurity_operators["particle"]))
                error("NaN in impurity operators")
            end
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
        @show spectrum_even
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
