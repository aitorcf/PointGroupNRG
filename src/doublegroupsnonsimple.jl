function cg_orbital_nonsimple( I_1::String , I_2::String , path ; verbose=false )

    # STRING version. 
    # Given two orbital irreps I_1 and I_2, it searches in path 
    # for the file containing CG information and returns it in 
    # the form of a dictionary:
    #
    #           cg[I_1,i_1,I_2,i_2,I_3,i_3,r_3] = ( I_1 , i_1 ; I_2 , i_2 | I_3 , i_3 , r_3 )
    #
    #
    cg::Dict{ Tuple{String,Int64,String,Int64,String,Int64,Int64} , ComplexF64 } = 
        Dict{ Tuple{String,Int64,String,Int64,String,Int64,Int64} , ComplexF64 }()

    file = [ x for x in readdir("$(path)/") 
               if (occursin("$(I_1)x$(I_2)",x) || occursin("$(I_2)x$(I_1)",x)) ][1]
    verbose && @show file 

    inverted = false
    occursin("$(I_2)x$(I_1)",file) && (inverted=true)

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

function get_cg_o_fullmatint_nonsimple( 
            cg_o_fulldict::Dict{Tuple{String,Int64,String,Int64,String,Int64,Int64},ComplexF64} , 
            oirreps::Vector{String} )::Dict{ NTuple{3,Int64} , Array{ComplexF64,4} } 
    # 
    # Given a CG coefficient dictionary in the form 
    #
    #           cg[I_1,i_1,I_2,i_2,I_3,i_3] = ( I_1 , i_1 ; I_2 , i_2 | I_3 , i_3 ),
    #
    # where I_1, I_2 and I_3 are STRINGS, it transforms it into a dictionary of matrices 
    #
    #           cg[I_1,I_2,I_3] = M(I_1,I_2,I_3), where
    #           M(I_1,I_2;I_3)_{i_1,i_2,i_3}=(I_1,i_1;I_2,i_2|I_3,i_3),
    #
    # where all indices are INTEGERS.
    #
    oirreps2indices = Dict( o=>i for (i,o) in enumerate(oirreps) )
    combs = Set( (k[1],k[3],k[5]) for k in keys(cg_o_fulldict) )

    function find_max_r( I1::String , I2::String , I3::String , 
                         cgodict::Dict{Tuple{String,Int64,String,Int64,String,Int64,Int64},ComplexF64} )
        matching_keys = Set( k for k in keys(cgodict) if (k[1],k[3],k[5])==(I1,I2,I3) )
        return maximum([k[end] for k in matching_keys])
    end


    cg_o_fullmat::Dict{Tuple{Int64,Int64,Int64},Array{ComplexF64,4}} = Dict() 
    for (Is1,Is2,Is3) in combs
        D1 = maximum(Set( k[2] for k in keys(cg_o_fulldict) if k[1]==Is1 ))
        D2 = maximum(Set( k[4] for k in keys(cg_o_fulldict) if k[3]==Is2 ))
        D3 = maximum(Set( k[6] for k in keys(cg_o_fulldict) if k[5]==Is3 ))
        R3 = find_max_r(Is1,Is2,Is3,cg_o_fulldict)
        I1 = oirreps2indices[Is1]
        I2 = oirreps2indices[Is2]
        I3 = oirreps2indices[Is3]
        push!( cg_o_fullmat , 
               (I1,I2,I3) => 
               ComplexF64[ get(cg_o_fulldict,(Is1,i1,Is2,i2,Is3,i3),zero(ComplexF64))
                           for r3=1:R3, i1=1:D1, i2=1:D2, i3=1:D3 ]
         )
    end
    return cg_o_fullmat 
end

function get_cg_o_info_nonsimple( 
            cg_o_dir::String , 
            atom_orbital_irreps::Vector{String} ;
            verbose=false )

    if verbose
        println( "GETTING ORBITAL CG INFORMATION" )
        println()
    end

    # all needed orbital irreps, their indices and dimensions.
    oirreps::Vector{String} = cg_shortcircuit_nonsimple( cg_o_dir , 
                                                         atom_orbital_irreps )::Vector{String}
    oirreps2indices::Dict{String,Int64} = Dict( o=>i for (i,o) in enumerate(oirreps) )
    if verbose 
        @show oirreps 
        @show oirreps2indices
    end

    # clebsch-gordan matrix
    cg_o_full::Dict{Tuple{String,Int64,String,Int64,String,Int64,Int64},ComplexF64} = get_cg_o_fulldict_nonsimple( oirreps , cg_o_dir )
    cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} = get_cg_o_fullmatint_nonsimple( cg_o_full , oirreps )
    if verbose 
        println( "Orbital CG matrix:" )
        print_dict( cg_o_fullmatint ) 
        println()
    end


    # orbital dimensions
    oirreps2dimensions::Dict{String,Int64} = Dict()
    for (ostring,oindex) in oirreps2indices
        for ((I1,I2,I3),IIImat) in cg_o_fullmatint 
            valid = false
            I1==oindex && (valid=true; i=1) 
            I2==oindex && (valid=true; i=2) 
            I3==oindex && (valid=true; i=3) 
            valid || continue
            oirreps2dimensions[ostring] = size(IIImat)[i+1]
        end
    end

    #oirreps2dimensions = Dict( "Eg" => 2 ,
    #                           "A1g"=> 1 ,
    #                           "A2g"=> 1 )
    oindex2dimensions::Vector{Int64} = collect( oirreps2dimensions[I] for I in oirreps )

    return (oirreps,
            oirreps2indices,
            oirreps2dimensions,
            oindex2dimensions,
            cg_o_fullmatint)

end

function nrg_full_doublegroups_nonsimple( 
            label::String ,
            calculation::String ,
            L::Float64 ,
            iterations::Int64 ,
            cutoff_type::String ,
            cutoff_magnitude ,
            cg_o_dir::String ,
            multiplets_dir::String ,
            impurity_config::Dict{String,Int64} ,
            shell_config::Dict{String,Int64} ,
            identityrep::String ,
            epsilon_symparams::Dict{ String , Vector{ComplexF64} } ,
            u_symparams::Dict{ Tuple{String,Float64} , Matrix{ComplexF64} } ,
            hop_symparams::Dict{ String , Matrix{ComplexF64} } ;
            distributed::Bool=false ,
            z::Float64=0.0 ,
            max_spin2::Int64=0 ,
            channels_dos::Dict{ String , Vector{Function} }=Dict{ String , Vector{Function} }() ,
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

    # defaults
    precompute_iaj = true
    spectral_method = "sakai1989"

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
    atom_orbital_irreps::Vector{String} = collect(keys(impurity_config))

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
    @show max_spin2
    @show spectral
    println( "OCCUPATION ENERGIES" )
    print_dict( epsilon_symparams ) 
    println( "COULOMB PARAMETERS" )
    print_dict( u_symparams ) 
    println( "HYBRIDIZATION PARAMETERS" )
    print_dict( hop_symparams )
    println()


    # hiztegia
    hiztegia = Dict{String,Any}( o=>o for (o,_) in impurity_config )
    merge!( hiztegia , Dict( "-"=>0.0 ) )

    #   ==========================   #
    #%% SYMMETRY-RELATED VARIABLES %%#
    #   ==========================   #

    # orbital symmetry
    (oirreps::Vector{String},
     oirreps2indices::Dict{String,Int64},
     oirreps2dimensions::Dict{String,Int64},
     oindex2dimensions::Vector{Int64},
     cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}}) = get_cg_o_info( cg_o_dir , atom_orbital_irreps )
    x=y

    # for dmnrg
    shell_dimension = reduce( * , [4^(oirreps2dimensions[I]*R) for (I,R) in shell_config] )

    # spin symmetry
    cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} = Dict( (0,0,0)=>[1.0;;;] )

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
    scale::Float64 = band_width
    if scale_asymptotic
        codiagonals_first = collect(values(channels_codiagonals[10]))[1][1]
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
    hop_symparams_int = Dict{Int64,Matrix{ComplexF64}}( oirreps2indices[k]=>(@.rescale(v,L,z,hopscale)) for (k,v) in hop_symparams )
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
            get_symstates_basis_multiplets_doublegroups( 
                    impurity_config,
                    oirreps2dimensions,
                    identityrep,
                    multiplets_dir,
                    cg_o_dir ;
                    verbose=true )
        orbital_multiplets = ordered_multiplets(multiplets_atom_noint)
        mult2index = Dict( m=>i for (i,m) in enumerate(orbital_multiplets))
        multiplets_atom::Set{NTuple{4,Int64}} = multiplets2int( multiplets_atom_noint , 
                                                                oirreps2indices )
        multiplets_a_atom::Set{NTuple{4,Int64}} = multiplets2int( multiplets_a_atom_noint , 
                                                                  oirreps2indices )
    else 
        multiplets_atom_noint = Set([(0,identityrep,0.0,1)]) 
        multiplets_atom = multiplets2int(multiplets_atom_noint,
                                                oirreps2indices)
        multiplets_a_atom = multiplets_atom
    end

    #   ------------------------ #
    #%% reduced pcg coefficients #
    #   ------------------------ #
    pcgred_atom::IntIrrepPCG = 
        calculation=="CLEAN" ? 
        IntIrrepPCG() :
        get_pcgred( basis_atom ,
                    symstates_atom_noint::ClearSymstateDict ,
                    multiplets_atom::IntMultipletSet ,
                    hiztegia ,
                    oirreps2indices::Dict{String,Int64} ,
                    cg_o_fullmatint::IntCG ,
                    cg_s_fullmatint::IntCG ;
                    verbose=false )::IntIrrepPCG

    #   ------------- #
    #%% impurity atom #
    #   ------------- #
    if calculation=="IMP"

        # operators
        epsilon::Operator{typeof(basis_atom)} = epsilon_sym( symstates_atom_noint , epsilon_symparams ; verbose=false )
        coulomb::Operator{typeof(basis_atom)} = u_sym( symstates_atom_noint , u_symparams ; verbose=false )

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
    println( "ATOMIC SPECTRUM" )
    println()
    print_spectrum( irrEU_clear )
    println()
    println( "------------------------------------" )
    println()
    irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} } = 
        irrEU2int( irrEU_clear , oirreps2indices ) 


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
        get_symstates_basis_multiplets_doublegroups(
                shell_config,
                oirreps2dimensions,
                identityrep,
                multiplets_dir,
                cg_o_dir ;
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
            counter_operator = electron_counter_sym( symstates_shell_noint ,
                                                     orbital )

            # diagonalize operator
            orbital_count_diagonalization = irrEU2int(get_irrEU_initial( symstates_shell_noint , counter_operator ),oirreps2indices)

            # introduce eigenvalues (number of particles) into dictionary
            push!( orbital_multiplets_occupations , Dict( G_mu=>N_mu for (G_mu,(N_mu,_)) in orbital_count_diagonalization) )

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

    #   ------------------------ #
    #%% reduced pcg coefficients #
    #   ------------------------ #
    pcgred_shell::Dict{NTuple{3,NTuple{3,Int64}},Array{ComplexF64,3}} = get_pcgred( 
                basis_shell ,
                symstates_shell_noint ,
                multiplets_shell ,
                hiztegia ,
                oirreps2indices ,
                cg_o_fullmatint ,
                cg_s_fullmatint ;
                verbose=true )


    #   ================================   #
    #%% COUPLING ATOM TO INNERMOST SHELL %%#
    #   ================================   #
    println()
    println( "::::::::::::::::::::::::::::::::::::::::" )
    println( "--- COUPLING ATOM TO INNERMOST SHELL ---" )
    println( "::::::::::::::::::::::::::::::::::::::::" )
    println()

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
    end

    #   ------------------------------ #
    #%% precompute clebsch-gordan sums #
    #   ------------------------------ #
    Bsum_o_dict,Bsum_s_dict,Csum_o_dict,Csum_s_dict =
        precompute_CGsums(
                oirreps ,
                union(multiplets_a_atom,multiplets_a_shell) ,
                union(multiplets_atom,multiplets_shell) ,
                max_spin2 ,
                oindex2dimensions ,
                cg_o_fullmatint ,
                cg_s_fullmatint ;
                verbose=true)
    Bsum_o_array,Bsum_s_array,Csum_o_array,Csum_s_array = 
        CGsums_dict2array( Bsum_o_dict,
                           Bsum_s_dict,
                           Csum_o_dict,
                           Csum_s_dict ) 

    #   -------- #
    #%% spectral #
    #   -------- #
    impurity_operators = Dict{String,Dict{IntTripleG,Array{ComplexF64,3}}}()
    spectral_functions = Dict{String,Dict{IntMultiplet,Matrix{Float64}}}()
    if spectral

        impurity_operators["particle"] = pcgred_atom
        GG_a  = Set(G_a for (_,G_a,_) in keys(pcgred_atom))
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
        
    #   ---------------   #
    #%% combinations u' %%#
    #   ---------------   #
    combinations_uprima = get_combinations_uprima_initial(
                identityrep ,
                calculation ,
                multiplets_atom ,
                oirreps2indices ;
                verbose=true )

    #   ---------------------------------------   #
    #%% matrix construction and diagonalization %%#
    #   ---------------------------------------   #
    (irrEU,combinations_uprima) = matdiag_redmat( 
                    multiplets_atom , 
                    multiplets_shell ,
                    irrEU , 
                    hop_symparams_int , 
                    keys_as_dict_o ,
                    keys_as_dict_s ,
                    Csum_o_array ,
                    Csum_s_array ,
                    Bsum_o_array ,
                    Bsum_s_array ,
                    pcgred_atom ,
                    pcgred_shell ,
                    collect(multiplets_a_atom) , 
                    collect(multiplets_a_shell) ,
                    combinations_uprima ;
                    conduction_diagonals=channels_diagonals[1],
                    verbose=false ,
                    distributed=distributed ,
                    precompute_iaj=precompute_iaj );
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
                        combinations_uprima )
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
                   precompute_iaj=precompute_iaj ,
                   compute_impmults=compute_impmults ,
                   mult2index=mult2index ,
                   orbital_multiplets=orbital_multiplets ,
                   mm_i=mm_i ,
                   channels_diagonals=channels_diagonals )
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
