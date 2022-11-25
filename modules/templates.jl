using Glob

        
function nrg_full_thermo( 
            label::String ,
            calculation::String ,
            L::Float64 ,
            z::Float64 ,
            distributed::Bool ,
            iterations::Int64 ,
            cutoff_type::String ,
            cutoff_magnitude::R ,
            max_spin2::Int64 ,
            cg_o_dir::String ,
            asym_dir::String ,
            atom_config::Dict{String,Int64} ,
            shell_config::Dict{String,Int64} ,
            identityrep::String ,
            epsilon_symparams::Dict{ Tuple{String,Int64} , ComplexF64 } ,
            u_symparams::Dict{ Tuple{String,Int64} , Matrix{ComplexF64} } ,
            hop_symparams::Dict{ String , Matrix{ComplexF64} } ;
            discretization="standard" ,
            distworkers::Int64=0 ,
            method::String="" ,
            minmult::Int64=0 ,
            mine::Float64=0.0 ,
            betabar::Float64=1.0 ,
            spectral::Bool=false ,
            etafac::Float64=1.0 ) where {R<:Real}

    if (spectral && calculation=="CLEAN") 
        println( "ERROR: calculation must be IMP for computing the spectral function" )
        return nothing 
    end
    
    # orbital irreps present in the atom
    atom_orbital_irreps::Vector{String} = collect(keys(atom_config))

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
    distributed && @show distworkers
    distributed && @show method
    println( "OCCUPATION ENERGIES" )
    print_dict( epsilon_symparams ) 
    println( "COULOMB PARAMETERS" )
    print_dict( u_symparams ) 
    println( "HYBRIDIZATION PARAMETERS" )
    print_dict( hop_symparams )
    println()
    
            
    # hiztegia
    hiztegia = Dict{String,Any}( o=>o for (o,_) in atom_config )
    merge!( hiztegia , Dict( "u"=>0.5, "d"=>-0.5 ) )

    #   ==========================   #
    #%% SYMMETRY-RELATED VARIABLES %%#
    #   ==========================   #

    # orbital symmetry
    (oirreps,
     oirreps2indices,
     oirreps2dimensions,
     oindex2dimensions,
     cg_o_fullmatint) = get_cg_o_info( cg_o_dir , atom_orbital_irreps )

    # spin symmetry
    cg_s_fullmatint = get_cg_s_fullmatint( max_spin2 );


    #   ===========   #
    #%% ATOMIC PART %%#
    #   ===========   #
    println()
    println( ":::::::::::::::::::" )
    println( "--- ATOMIC PART ---" )
    println( ":::::::::::::::::::" )
    println()
    
    #   ------------------- #
    #%% rescaled parameters #
    #   ------------------- #
    hop_symparams     = Dict( oirreps2indices[k]=>@.rescale(v,L,z,discretization) for (k,v) in hop_symparams )
    epsilon_symparams = Dict( k=>@.rescale(v,L,z,discretization) for (k,v) in epsilon_symparams )
    u_symparams       = Dict( k=>@.rescale(v,L,z,discretization) for (k,v) in u_symparams )
    println( "RESCALED PARAMETERS FOR H0" )
    @show epsilon_symparams 
    @show u_symparams 
    @show hop_symparams
    println()

    #   ------------------------------- #
    #%% symstates, basis and multiplets #
    #   ------------------------------- #
    if calculation=="IMP"
        symstates_atom_noint,basis_atom,multiplets_atom_noint,multiplets_a_atom_noint = 
            get_symstates_basis_multiplets( 
                    atom_config,
                    oirreps2dimensions,
                    identityrep,
                    asym_dir,
                    cg_o_dir ;
                    verbose=true )
        omults = ordered_multiplets(multiplets_atom_noint)
        mult2index = Dict( m=>i for (i,m) in 
                           enumerate(omults))
        global multiplets_atom = multiplets2int( multiplets_atom_noint , 
                                                 oirreps2indices )
        global multiplets_a_atom = multiplets2int( multiplets_a_atom_noint , 
                                                   oirreps2indices )
    else 
        multiplets_atom_noint = Set([(0,identityrep,0.0,1)]) 
        global multiplets_atom = multiplets2int(multiplets_atom_noint,
                                                oirreps2indices)
        global multiplets_a_atom = multiplets_atom
    end

    #   ------------------------ #
    #%% reduced pcg coefficients #
    #   ------------------------ #
    pcgred_atom = 
        calculation=="CLEAN" ? 
        Dict{NTuple{3,NTuple{3,Int64}},Array{ComplexF64,3}}() :
        get_pcgred( basis_atom ,
                    symstates_atom_noint ,
                    multiplets_atom ,
                    hiztegia ,
                    oirreps2indices ,
                    cg_o_fullmatint ,
                    cg_s_fullmatint ;
                    verbose=false )
        
    #   ------------- #
    #%% impurity atom #
    #   ------------- #
    if calculation=="IMP"

        # operators
        epsilon = epsilon_sym( symstates_atom_noint , epsilon_symparams ; verbose=false )
        coulomb = u_sym( symstates_atom_noint , u_symparams ; verbose=false )

        # hamiltonian 
        global H = epsilon + coulomb 

    end

    #   -----   #
    #%% irreu %%#
    #   -----   #
    if calculation=="IMP"
        global irrEU = get_irrEU_initial(symstates_atom_noint,H;verbose=true)
    elseif calculation=="CLEAN" 
        global irrEU = get_irrEU_initial(identityrep,oirreps2indices)
    end
    print_spectrum( irrEU )
    irrEU = irrEU2int( irrEU , oirreps2indices ) 


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
    symstates_shell_noint,
    basis_shell,
    multiplets_shell_noint,
    multiplets_a_shell_noint = 
        get_symstates_basis_multiplets( 
                shell_config,
                oirreps2dimensions,
                identityrep,
                asym_dir,
                cg_o_dir ;
                verbose=true )
    multiplets_shell = multiplets2int( multiplets_shell_noint , 
                                       oirreps2indices )
    multiplets_a_shell = multiplets2int( multiplets_a_shell_noint , 
                                         oirreps2indices )

    #   ------------------------ #
    #%% reduced pcg coefficients #
    #   ------------------------ #
    pcgred_shell = get_pcgred( 
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
    mm_i,m_imp = setup_impmultinfo( 
                    multiplets_atom ,
                    irrEU ,
                    betabar ,
                    oindex2dimensions )
    println( "IMPURITY COMPOSITION" )
    @show m_imp
    println()

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
                cg_s_fullmatint )
    Bsum_o_array,Bsum_s_array,Csum_o_array,Csum_s_array = 
        CGsums_dict2array( Bsum_o_dict,
                           Bsum_s_dict,
                           Csum_o_dict,
                           Csum_s_dict ) 

    #   -------- #
    #%% spectral #
    #   -------- #
    if spectral 

        global M = pcgred_atom 
        part0 = get_partition0(irrEU,oindex2dimensions)
        A = redM2A( M,
                    collect(multiplets_a_atom),
                    cg_o_fullmatint,
                    cg_s_fullmatint,
                    irrEU,
                    part0)
        global AA = [A]

        Mo_tot = length(oirreps2indices) 
        II_a = collect(Set([G[2] for G in get_irreps( multiplets_a_atom )]))
        Ms_atomspin = maximum([m[3] for m in multiplets_atom])
        Ms_shellspin = maximum([m[3] for m in multiplets_shell]) 
        Ms_tot = maximum((Ms_atomspin,Ms_shellspin))
        global Karray_orbital,Karray_spin = 
                    compute_Ksum_arrays(
                        oindex2dimensions,
                        cg_o_fullmatint,
                        cg_s_fullmatint,
                        Mo_tot ,
                        II_a ,
                        max_spin2 ,
                        Ms_tot )
        global alpha = compute_ebar0_z( z , L ; discretization=discretization )
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
                    hop_symparams , 
                    cg_o_fullmatint , 
                    cg_s_fullmatint ,
                    Csum_o_array ,
                    Csum_s_array ,
                    Bsum_o_array ,
                    Bsum_s_array ,
                    pcgred_atom ,
                    pcgred_shell ,
                    collect(multiplets_a_atom) , 
                    collect(multiplets_a_shell) ,
                    combinations_uprima , 
                    oindex2dimensions ;
                    verbose=false ,
                    distributed=distributed );
    print( "AFTER ADDING INNERMOST SHELL, " )
    print_spectrum( irrEU )

    #   --------------------------- #
    #%% update impurity information # 
    #   --------------------------- #
    mm_i,m_imp = update_impmultinfo( 
                    mm_i ,
                    irrEU ,
                    betabar ,
                    oindex2dimensions ,
                    combinations_uprima )
    if spectral 
        global M, AA = update_redmat_AA_CGsummethod(
                M,
                irrEU ,
                combinations_uprima ,
                collect(multiplets_a_atom) ,
                cg_o_fullmatint ,
                cg_s_fullmatint ,
                Karray_orbital ,
                Karray_spin ,
                AA ,
                oindex2dimensions ;
                verbose=false )
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
        nrg = NRG( iterations,
                   cutoff_type,
                   cutoff_magnitude,
                   L,
                   hop_symparams,
                   irrEU,
                   multiplets_shell,
                   cg_o_fullmatint,
                   cg_s_fullmatint,
                   Csum_o_array ,
                   Csum_s_array ,
                   Bsum_o_array ,
                   Bsum_s_array ,
                   pcgred_shell,
                   collect(multiplets_a_shell), 
                   combinations_uprima,
                   betabar,
                   oindex2dimensions,
                   mm_i ;
                   mine=mine ,
                   distributed=distributed ,
                   method=method ,
                   z=z ,
                   discretization=discretization ,
                   verbose=false )
    else 
        nrg = NRG( iterations,
                   cutoff_type,
                   cutoff_magnitude,
                   L,
                   hop_symparams,
                   irrEU,
                   multiplets_shell,
                   cg_o_fullmatint,
                   cg_s_fullmatint,
                   Csum_o_array ,
                   Csum_s_array ,
                   Bsum_o_array ,
                   Bsum_s_array ,
                   pcgred_shell,
                   collect(multiplets_a_shell), 
                   combinations_uprima,
                   betabar,
                   oindex2dimensions,
                   mm_i ;
                   mine=mine ,
                   distributed=distributed ,
                   method=method ,
                   z=z ,
                   discretization=discretization ,
                   verbose=false ,
                   spectral=true ,
                   M=M,
                   AA=AA , 
                   etafac=etafac ,
                   Karray_orbital=Karray_orbital ,
                   Karray_spin=Karray_spin ,
                   multiplets_atomhop=collect(multiplets_a_atom) ,
                   alpha=Float64(alpha) )
    end

    println()

    #   ===========   #
    #%% PERFORMANCE %%#
    #   ===========   #

    println( "===========" )
    println( "PERFORMANCE" )
    println( "===========" )
    print_performance_onestep(nrg)
    println()

    #   ==============   #
    #%% SAVING TO FILE %%#
    #   ==============   #

    # impurity properties 
    if calculation=="IMP" 
        write_impurity_info( nrg , omults , mult2index , label , z )
    end

    if spectral 
        writedlm( "spectral/spectral.dat" , nrg.specfunc )
    end

    # thermodata for this given value of z
    write_thermodata_onez( nrg , calculation , label , z )

    # thermo diff
    if calculation=="IMP"
        if length(glob("thermodata/thermo_clean_$(label)_z$z.dat"))!==0 
            write_thermodiff( label , z )
        end
    end

    spectral && writedlm( "spectral/spectral.dat" , nrg.specfunc )

end

function multiplets_2part( 
            cg_o_dir ,
            asym_dir ,
            atom_orbital_irreps ,
            atom_config ,
            identityrep )

    # orbital symmetry
    (oirreps,
     oirreps2indices,
     oirreps2dimensions,
     oindex2dimensions,
     cg_o_fullmatint) = get_cg_o_info( cg_o_dir , atom_orbital_irreps )

    # spin symmetry
    cg_s_fullmatint = get_cg_s_fullmatint( max_spin2 );

    
    symstates_atom_noint,basis_atom,multiplets_atom_noint,multiplets_a_atom_noint = 
        get_symstates_basis_multiplets( 
                atom_config,
                oirreps2dimensions,
                identityrep,
                asym_dir,
                cg_o_dir ;
                verbose=true )
    multiplets_atom_twopart = filter( 
            m->m[1]==2 ,
            multiplets_atom_noint )

    println( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" )
    println( "TWO-PARTICLE ATOMIC MULTIPLETS (with representative)" )
    println( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" )
    for ma in multiplets_atom_twopart 
        println( ma )
        println( symstates_atom_noint[(ma[1:3]...,1,ma[3],ma[4])] )
        println()
    end
end
