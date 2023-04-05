using Glob

include( "shell.jl" )
include( "lanczos.jl" )

function get_cg_o_info( 
            cg_o_dir::String , 
            atom_orbital_irreps::Vector{String} ;
            verbose=false )

    if verbose
        println( "GETTING ORBITAL CG INFORMATION" )
        println()
    end

    # all needed orbital irreps, their indices and dimensions.
    oirreps = cg_shortcircuit( cg_o_dir , 
                               atom_orbital_irreps )::Vector{String}
    oirreps2indices = Dict( o=>i for (i,o) in enumerate(oirreps) )
    if verbose 
        @show oirreps 
        @show oirreps2indices
    end

    # clebsch-gordan matrix
    cg_o_full = get_cg_o_fulldict( oirreps , cg_o_dir )
    cg_o_fullmatint = get_cg_o_fullmatint( cg_o_full , oirreps )
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
            oirreps2dimensions[ostring] = size(IIImat)[i]
        end
    end

    #oirreps2dimensions = Dict( "Eg" => 2 ,
    #                           "A1g"=> 1 ,
    #                           "A2g"=> 1 )
    oindex2dimensions = collect( oirreps2dimensions[I] for I in oirreps )

    return (oirreps,oirreps2indices,oirreps2dimensions,oindex2dimensions,cg_o_fullmatint)

end

function get_symstates_basis_multiplets( 
            atom_config::Dict{String,Int64},
            oirreps2dimensions::Dict{String,Int64} ,
            identityrep::String ,
            asym_dir::String ,
            cg_o_dir::String ;
            verbose=true )

    # hiztegia
    hiztegia::Dict{String,Any} = Dict( o=>o for (o,_) in atom_config )
    merge!( hiztegia , Dict( "u"=>0.5, "d"=>-0.5 ) )

    # no impurity (clean)
    if length(atom_config)==0 
        hilbert = HilbertSpace( ) 
    end

    # symstates
    separate_symstates = []
    for (oirrep,multiplicity) in atom_config
        for m in 1:multiplicity
            onemult_states = get_atomic_states(m,oirrep,oirreps2dimensions)
            onemult_hilbert = HilbertSpace( onemult_states )
            onemult_symstates = oneirrep_symstates( 
                                    onemult_hilbert ,
                                    hiztegia ,
                                    identityrep ,
                                    "$(asym_dir)$(oirrep)_julia/" )
            push!( separate_symstates , onemult_symstates )
        end
    end
    symstates = separate_symstates[1]
    for symstates_new in separate_symstates[2:end]
        symstates = cg_reduce_product_states(
                        symstates , 
                        symstates_new ,
                        cg_o_dir )
    end
    
    # basis 
    basis = collect(values(symstates))[1].basis 

    # multiplets 
    multiplets = get_multiplets( symstates )
    multiplets_a = Set( m for m in multiplets if m[1]==1 )

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
function get_symstates_basis_multiplets( 
            shell_config::Dict{String,Vector{Int64}},
            oirreps2dimensions::Dict{String,Int64} ,
            identityrep::String ,
            asym_dir::String ,
            cg_o_dir::String ;
            verbose=true )

    # hiztegia
    hiztegia::Dict{String,Any} = Dict( o=>o for (o,_) in atom_config )
    merge!( hiztegia , Dict( "u"=>0.5, "d"=>-0.5 ) )

    # symstates
    separate_symstates = []
    for (oirrep,multiplicities) in atom_config
        for m in multiplicities
            onemult_states = get_atomic_states(m,oirrep,oirreps2dimensions)
            onemult_hilbert = HilbertSpace( onemult_states )
            onemult_symstates = oneirrep_symstates( 
                                    onemult_hilbert ,
                                    hiztegia ,
                                    identityrep ,
                                    "$(asym_dir)$(oirrep)_julia/" )
            push!( separate_symstates , onemult_symstates )
        end
    end
    symstates = separate_symstates[1]
    for symstates_new in separate_symstates[2:end]
        symstates = cg_reduce_product_states(
                        symstates , 
                        symstates_new ,
                        cg_o_dir )
    end
    
    # basis 
    basis = collect(values(symstates))[1].basis 

    # multiplets 
    multiplets = get_multiplets( symstates )
    multiplets_a = Set( m for m in multiplets if m[1]==1 )

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

function rescale( 
            num::N ,
            L::Float64 ,
            z::Float64 ,
            discretization::String ;
            iterations::Int64=30 ,
            eta::Function=x->1.0 ) where {N<:Number}
    # apply initial rescaling to parameter A:
    #
    #   A → A/( √Λ ϵ^z_0(Λ) ),
    #
    # where ϵ^z_0(Λ) is the hopping to the first
    # shell (not the innermost one, but the next one).
    
    # α = ϵ_0^z (scaling factor)
    if discretization!=="lanczos" 
        a::Float64 = compute_ebar0_z(
                         z,
                         L;
                         discretization=discretization)
    elseif discretization=="lanczos" 
        a = get_hoppings( iterations , L , z , eta )[2][1]
    end

    return num/(a*sqrt(L))
end
function rescale( 
            num::N ,
            L::Float64 ,
            z::Float64 ,
            scale::Float64 ) where {N<:Number}
    # apply initial rescaling to parameter A:
    #
    #   A → A/( √Λ ϵ^z_0(Λ) ),
    #
    # where ϵ^z_0(Λ) is the hopping to the first
    # shell (not the innermost one, but the next one).
    #
    return num/(scale*sqrt(L))
end

# IMP method
function get_irrEU_initial( 
            symstates_0::Dict{T,S} , 
            H::O ;
            verbose=false ) where {T<:Tuple,S<:State,O<:Operator}

    irreps_0 = get_irreps( symstates_0 )
    irrEU_imp = symdiag( irreps_0 , symstates_0 , H )
    if verbose
        println( "irrEU (for IMP)" )
        print_dict( irrEU_imp )
        println()
    end
    irrEU_imp = normalize_irrEU( irrEU_imp )
    return irrEU_imp
end
# CLEAN method
function get_irrEU_initial( 
            identityrep::String ,
            oirreps2indices::Dict{String,Int64} ;
            verbose=false )

    irrEU_clean = get_irrEU_clean( identityrep )

    if verbose 
        println( "irrEU (for CLEAN)" )
        print_dict( irrEU_clean )
        println()
    end

    return irrEU_clean
end

function get_combinations_uprima_initial( 
            identityrep::String ,
            calculation::String ,
            multiplets_0 ,
            oirreps2indices::Dict{String,Int64} ;
            verbose=true )

    combinations_uprima = 
        Dict{ Tuple{Int64,Int64,Int64,Int64} , NTuple{2,Tuple{Int64,Int64,Int64,Int64}} }()

    m_vac = (0,oirreps2indices[identityrep],0,1)
    if calculation=="IMP"
        for m_mu in multiplets_0
            push!( combinations_uprima , m_mu=>(m_mu,m_vac) )
        end 
    elseif calculation=="CLEAN" 
        push!( combinations_uprima , m_vac=>(m_vac,m_vac) )
    end

    irreps_uprima = Set( k[1:3] for k in keys(combinations_uprima) )
    combinations_uprima = 
            Dict{ NTuple{3,Int64} , Vector{NTuple{3,NTuple{4,Int64}}} }(
                G => NTuple{3,NTuple{4,Int64}}[
                    (m_u,m_mu,m_i)
                    for (m_u,(m_mu,m_i)) in combinations_uprima 
                    if m_u[1:3]==G
                    ]
                for G in irreps_uprima
            ) 
    if verbose 
        println( "COMBINATIONS U' FOR N=0" )
        for (G,combs) in combinations_uprima
            println( "$G => $combs" )
        end
        println()
    end

    return combinations_uprima 
end

function irrEU2int( irrEU , oirreps2indices )
    return Dict( (convert_to_int(G,oirreps2indices),(E,U)) 
                 for (G,(E,U)) in irrEU )
end
function multiplets2int( multiplets , oirreps2indices )
    return Set( convert_to_int(m,oirreps2indices) for m in multiplets )
end

function get_multiplets_block( 
            calculation ,
            multiplets_0 ,
            identityrep ,
            oirreps2indices )
    if calculation=="IMP" 
        return multiplets_0 
    elseif calculation=="CLEAN" 
        return Set([(0,oirreps2indices[identityrep],0,1)])
    end
end

function setup_impmultinfo( 
            multiplets_block ,
            irrEU ,
            betabar ,
            oindex2dimensions )
    
    omults = ordered_multiplets(multiplets_block)
    mult2index = Dict( m=>i for (i,m) in 
                       enumerate(omults))
    mm_i = Dict( 
                m=>[(i==mult2index[m] ? 1.0 : 0.0)
                    for i in 1:length(multiplets_block)] 
                    for m in omults
               )
    m_imp = mult_thermo( irrEU ,
                         betabar ,
                         oindex2dimensions ,
                         mm_i )

    return mm_i,m_imp
end
function update_impmultinfo( 
                mm_i ,
                irrEU ,
                betabar ,
                oindex2dimensions ,
                combinations_uprima )

    mm_i = imp_mults( irrEU ,
                      oindex2dimensions ,
                      combinations_uprima ,
                      mm_i )
    m_imp = mult_thermo( irrEU ,
                         betabar ,
                         oindex2dimensions ,
                         mm_i )
    return mm_i,m_imp 
end

function check_pcgred( 
            pcgred ,
            pcgdict ,
            cg_o_fullmatint ,
            cg_s_fullmatint ,
            symstates ,
            shell_cops )

    for ((G_i,G_a,G_j),redmat) in pcgred 

        for r_i in 1:size(redmat,1),
            r_a in 1:size(redmat,2),
            r_j in 1:size(redmat,3)

            (N_i,I_i,S_i) = G_i[1:3]
            (N_a,I_a,S_a) = G_a[1:3]
            (N_j,I_j,S_j) = G_j[1:3]

            cgomat = cg_o_fullmatint[(I_a,I_j,I_i)]
            cgsmat = cg_s_fullmatint[(S_a,S_j,S_i)]

            (dI_a,dI_j,dI_i) = size(cgomat)
            (dS_a,dS_j,dS_i) = size(cgsmat)

            for i_i  in 1:dI_i,
                (si_i,s_i) in enumerate(-S_i:2:S_i),
                i_j  in 1:dI_j,
                (si_j,s_j) in enumerate(-S_j:2:S_j),
                i_a  in 1:dI_a,
                (si_a,s_a) in enumerate(-S_a:2:S_a)

                q_i = (G_i...,i_i,s_i,r_i)
                q_a = (G_a...,i_a,s_a,r_a)
                q_j = (G_j...,i_j,s_j,r_j)

                dictmatel = get( pcgdict ,(q_i,q_a,q_j), zero(ComplexF64) )
                cg = conj(cgomat[i_a, i_j, i_i ]*cgsmat[si_a,si_j,si_i])
                redmatel = redmat[r_i, r_a, r_j ]*cg

                if !(dictmatel≈redmatel)

                    println( "''''''''''''''''''''" )
                    println( "not matching" )
                    @show G_i,G_a,G_j
                    @show dictmatel
                    @show redmat[r_i,r_a,r_j]
                    @show cg
                    @show redmatel 
                    println( "pcgdict states" )
                    println()
                    print( "| $q_j ) = " ) 
                    println( symstates[q_j] )
                    println()
                    print( "f^†( $q_a ) * | $q_j ) = " )
                    println( shell_cops[q_a] * symstates[q_j] ) 
                    println()
                    print( "( $q_i | = " ) 
                    println( symstates[q_i] )
                    println( "''''''''''''''''''''" )

                end
                
            end
            end
    end

end
function prepare_pcgred( 
            symstates_0 ,
            basis_0 ,
            hiztegia ,
            oirreps2indices ,
            cg_o_fullmatint ,
            cg_s_fullmatint ,
            multiplets_block ,
            multiplets_shell ;
            verbose=false )

    # pcg in dict format
    pcg = get_pseudoCG( symstates_0 , 
                        basis_0 , 
                        hiztegia , 
                        oirreps2indices )
    #if verbose 
    #    println( "PCG DICT" )
    #    print_dict( pcg ) 
    #    println()
    #end

    # pcg in reduced matrix format
    multiplets_a = collect(filter( x->x[1]==1 , multiplets_shell ))
    pcgred_atom = get_redmat2( pcg ,
                              multiplets_shell ,
                              multiplets_a ,
                              cg_o_fullmatint ,
                              cg_s_fullmatint ;
                              verbose=false )
    pcgred_shell = get_redmat2( pcg ,
                               multiplets_shell ,
                               multiplets_a ,
                               cg_o_fullmatint ,
                               cg_s_fullmatint ;
                               verbose=false )
    if verbose
        println( "pcgred atom" )
        print_dict( pcgred_atom )
        println()
        println( "pcgred shell" )
        print_dict( pcgred_shell )
        println()
    end

    return (pcgred_atom,pcgred_shell)
end

function get_pcgred( 
            basis ,
            symstates_noint ,
            multiplets,
            hiztegia ,
            oirreps2indices ,
            cg_o_fullmatint ,
            cg_s_fullmatint ;
            verbose=false ) 

    # pcg in dict format
    pcg = get_pseudoCG( symstates_noint , 
                        basis , 
                        hiztegia , 
                        oirreps2indices )
    #if verbose 
    #    println( "PCG DICT" )
    #    print_dict( pcg ) 
    #    println()
    #end

    # pcg in reduced matrix format
    multiplets_a = collect(filter( x->x[1]==1 , multiplets ))
    pcgred = get_redmat2( pcg ,
                          multiplets ,
                          multiplets_a ,
                          cg_o_fullmatint ,
                          cg_s_fullmatint ;
                          verbose=false )
    if verbose
        println( "PCGRED" )
        print_dict( pcgred )
        println()
    end

    return pcgred
end

function precompute_CGsums(
            oirreps::Vector{String} ,
            multiplets_a ,
            multiplets_all ,
            max_spin2 ,
            oindex2dimensions ,
            cg_o_fullmatint ,
            cg_s_fullmatint ;
            verbose=false )
            
    Mo_tot    = length( oirreps ) 
    II_a      = collect(Set( m[2] for m in multiplets_a ))
    Ms_tot    = max_spin2
    Ms_shell  = maximum(collect( m[3] for m in multiplets_all ))

    Csum_o_dict = compute_CG_Csumdict_o( 
                        oindex2dimensions ,
                        cg_o_fullmatint ,
                        Mo_tot ,
                        II_a )
    Csum_s_dict = compute_CG_Csumdict_s(
                        cg_s_fullmatint ,
                        Ms_tot ,
                        Ms_shell )
    Bsum_o_dict = compute_CG_Bsumdict_o(
                        oindex2dimensions ,
                        cg_o_fullmatint ,
                        Mo_tot ,
                        II_a )
    Bsum_s_dict = compute_CG_Bsumdict_s( 
                        cg_s_fullmatint ,
                        Ms_tot ,
                        Ms_shell )

    if verbose 
        println( "B sum dict | orbital" )
        print_dict( Bsum_o_dict )
        println( "B sum dict | spin" )
        print_dict( Bsum_s_dict )
        println( "C sum dict | orbital" )
        print_dict( Csum_o_dict )
        println( "C sum dict | spin" )
        print_dict( Csum_s_dict )
    end

    return (Bsum_o_dict,Bsum_s_dict,Csum_o_dict,Csum_s_dict)
end
function CGsums_dict2array(
            Bsum_o_dict,
            Bsum_s_dict,
            Csum_o_dict,
            Csum_s_dict ;
            verbose=false )
    
    # B orbital 
    Bodims1 = maximum([ x[1] for x in keys(Bsum_o_dict) ])
    Bodims2 = maximum([ x[2] for x in keys(Bsum_o_dict) ])
    Bodims3 = maximum([ x[3] for x in keys(Bsum_o_dict) ])
    Bodims4 = maximum([ x[4] for x in keys(Bsum_o_dict) ])
    Bodims5 = maximum([ x[5] for x in keys(Bsum_o_dict) ])
    Bodims6 = maximum([ x[6] for x in keys(Bsum_o_dict) ])
    # C orbital 
    Codims1 = maximum([ x[1] for x in keys(Csum_o_dict) ])
    Codims2 = maximum([ x[2] for x in keys(Csum_o_dict) ])
    Codims3 = maximum([ x[3] for x in keys(Csum_o_dict) ])
    Codims4 = maximum([ x[4] for x in keys(Csum_o_dict) ])
    Codims5 = maximum([ x[5] for x in keys(Csum_o_dict) ])
    Codims6 = maximum([ x[6] for x in keys(Csum_o_dict) ])
    # B spin 
    Bsdims1 = maximum([ (x[1]+1) for x in keys(Bsum_s_dict) ])
    Bsdims2 = maximum([ (x[2]+1) for x in keys(Bsum_s_dict) ])
    Bsdims3 = maximum([ (x[3]+1) for x in keys(Bsum_s_dict) ])
    Bsdims4 = maximum([ (x[4]+1) for x in keys(Bsum_s_dict) ])
    Bsdims5 = maximum([ (x[5]+1) for x in keys(Bsum_s_dict) ])
    Bsdims6 = maximum([ (x[6]+1) for x in keys(Bsum_s_dict) ])
    # C spin
    Csdims1 = maximum([ (x[1]+1) for x in keys(Csum_s_dict) ])
    Csdims2 = maximum([ (x[2]+1) for x in keys(Csum_s_dict) ])
    Csdims3 = maximum([ (x[3]+1) for x in keys(Csum_s_dict) ])
    Csdims4 = maximum([ (x[4]+1) for x in keys(Csum_s_dict) ])
    Csdims5 = maximum([ (x[5]+1) for x in keys(Csum_s_dict) ])
    Csdims6 = maximum([ (x[6]+1) for x in keys(Csum_s_dict) ])
    Bsum_o_array = zeros( ComplexF64 , Bodims1 , Bodims2 , Bodims3 , Bodims4 , Bodims5 , Bodims6 )
    Bsum_s_array = zeros( ComplexF64 , Bsdims1 , Bsdims2 , Bsdims3 , Bsdims4 , Bsdims5 , Bsdims6 )
    Csum_o_array = zeros( ComplexF64 , Codims1 , Codims2 , Codims3 , Codims4 , Codims5 , Codims6 )
    Csum_s_array = zeros( ComplexF64 , Csdims1 , Csdims2 , Csdims3 , Csdims4 , Csdims5 , Csdims6 )
    for (k,v) in Bsum_o_dict 
        Bsum_o_array[k...] = v
    end
    for (k,v) in Bsum_s_dict 
        Bsum_s_array[(k.+1)...] = v
    end
    for (k,v) in Csum_o_dict 
        Csum_o_array[k...] = v
    end
    for (k,v) in Csum_s_dict 
        Csum_s_array[(k.+1)...] = v
    end

    return (Bsum_o_array,Bsum_s_array,Csum_o_array,Csum_s_array)
end

function print_spectrum( 
            irrEU ;
            savefile=false ,
            label="" )

    spectrum = [ ((G...,r),e) 
                 for (G,(E,U)) in irrEU 
                 for (r,e) in enumerate(E) ]
    sort!( spectrum , by=x->x[2] )

    println( "SPECTRUM" )
    for (m,e) in spectrum
        @show m, e
    end
    println()

    if savefile 
        open( "spectrum_$label.dat" , "w" ) do f
            for (m,e) in spectrum
                write( f , "$m  |  $e\n" )
            end
        end
    end

end

function print_irrep_products( cg_o_fullmatint , oirreps2indices ) 
    oindices2irreps = Dict( v=>k for (k,v) in oirreps2indices )
    combs = Set((x[1],x[2]) for x in keys(cg_o_fullmatint) )
    products = Dict( c=>[] for c in combs )
    for k in keys(cg_o_fullmatint),
        c in combs
        k[1:2]==c && push!(products[c],k[3]) 
    end
    println( "ORBITAL IRREP PRODUCTS" )
    for ((G1,G2),GG3) in products 
        I1 = oindices2irreps[G1]
        I2 = oindices2irreps[G2]
        II3 = [oindices2irreps[G] for G in GG3]
        print( "$I1 ⊗ $I2 = " )
        for I3 in II3[1:(end-1)]
            print( "$I3 ⊕ " )
        end
        println( "$(II3[end])" )
    end
    println()

end

function construct_custom_orbital_space( 
            external_label::Int64,
            orbital_names::Vector{String},
            oirreps2dimensions::Dict{String,Int64},
            identityrep::String,
            hiztegia::D,
            CG_PATH::String,
            ASYM_PATH::String
            ) where {D<:Dict}

    seed_tuples = [(external_label,o) for o in orbital_names]

    symstates_orbitals = []
    for (el,o) in seed_tuples

        oirrep = hiztegia[o]
        state_tuples = [(el,o,i,m) for i in 1:oirreps2dimensions[oirrep]
                                   for m in ["u","d"]]
        hilbert_orbital = HilbertSpace(state_tuples)

        irrep = hiztegia[o]
        symstates_nor = oneirrep_symstates( 
                            hilbert_orbital , 
                            hiztegia ,
                            identityrep ,
                            "$(ASYM_PATH)$(irrep)_julia/" ) 

        push!(
            symstates_orbitals ,
            Dict( (q[1:5]...,1)=>s 
                  for (q,s) in symstates_nor )
        )
    end

    symstates_total = symstates_orbitals[1]
    for i in 2:length(symstates_orbitals)
        symstates_total = cg_reduce_product_states(
                            symstates_total ,
                            symstates_orbitals[i] ,
                            CG_PATH )
    end
    
    # basis 
    basis = collect(values(symstates_total))[1].basis 
    
    # multiplets 
    multiplets = get_multiplets( symstates_total )

    return (basis,multiplets,symstates_total)

end
        
function nrg_full( 
            label::String ,
            calculation::String ,
            L::Float64 ,
            iterations::Int64 ,
            cutoff_type::String ,
            cutoff_magnitude::R ,
            cg_o_dir::String ,
            asym_dir::String ,
            atom_config::Dict{String,Int64} ,
            shell_config::Dict{String,Int64} ,
            identityrep::String ,
            epsilon_symparams::Dict{ String , Vector{ComplexF64} } ,
            u_symparams::Dict{ Tuple{String,Int64} , Matrix{ComplexF64} } ,
            hop_symparams::Dict{ String , Matrix{ComplexF64} } ;
            distributed::Bool=false,
            z::Float64=0.0 ,
            max_spin2::Int64=10 ,
            channel_etas::Dict{ String , Vector{Function} }=Dict{ String , Vector{Function} }() ,
            discretization="standard" ,
            distworkers::Int64=0 ,
            method::String="" ,
            minmult::Int64=0 ,
            mine::Float64=0.0 ,
            betabar::Float64=1.0 ,
            spectral::Bool=false ,
            etafac::Float64=1.0 ,
            Nz=1 ) where {R<:Real}

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

    # ------------------------------------ #
    # improve dictionary lookup for pcgred #
    # ------------------------------------ #
    #max_population_atom  = sum( 2*O*oirreps2dimensions[I] for (I,O) in atom_config )
    #max_population_shell = sum( 2*O*oirreps2dimensions[I] for (I,O) in shell_config )
    #Ndim = max_population_atom + iterations*max_population_shell
    #Idim = length(oirreps)
    #Sdim = max_spin2
    #maxISdim = maximum((Idim,Sdim))
    #ISdim = Idim*Sdim
    #Gdim = Ndim*Idim*Sdim
    #dimtup = ( ISdim , Sdim , 1 ) 
    #begin global
    #    @eval function Base.hash( x::NTuple{3,NTuple{3,Int64}} )
    #        UInt64(sum( sum.((x[1].*$dimtup,x[2].*$dimtup,x[3].*$dimtup)).*($(Gdim^2),$Gdim,1) ))
    #    end
    #end

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

    # default behavior: eta=x->1
    if length(channel_etas)==0 
        channel_etas = Dict( k=>Function[x->0.5 for i in 1:size(v)[1]]
                             for (k,v) in hop_symparams 
        )
    end
    
    # symmetry-wise channel structure 
    #       irrep => [ ϵ , ̄ϵ , ξ ]
    channel_symstructure = Dict( 
            k=>[get_hoppings(iterations,L,z,coupling)
                for coupling in couplings]
            for (k,couplings) in channel_etas
    )
    
    # symmetry-wise etabar 
    etabar_sym = Dict( 
            k=>[s[4] for s in v]
            for (k,v) in channel_symstructure 
    )
    @show etabar_sym
    println()

    # scale parameter: first asymptotic hopping element 
    #       irrep => [ ̄ϵ[1] ]
    scale_symparams = Dict( 
            k=>[h[2][1] for h in s] 
            for (k,s) in channel_symstructure 
    )
    @show scale_symparams
    println()
    scale = maximum([v for (k,V) in scale_symparams for v in V])
    factor_symparams = Dict( k=>v./scale for (k,v) in scale_symparams )
    @show factor_symparams
    println()

    # hopping parameters: ξ = ϵ / ̄ϵ
    #       irrep => [ ξ ]
    xi_symparams::Dict{Int64,Vector{Vector{ComplexF64}}} = Dict( 
            oirreps2indices[k]=>[ComplexF64.(h[3].*factor_symparams[k][i]) for (i,h) in enumerate(s)] 
            for (k,s) in channel_symstructure
    )
    
    #   ------------------- #
    #%% rescaled parameters #
    #   ------------------- #
    oindices2irreps = Dict( v=>k for (k,v) in oirreps2indices )
    println()
    if discretization=="lanczos"
        hop_symparams     = Dict( oirreps2indices[k]=>(@.rescale(v,L,z,scale)) for (k,v) in hop_symparams )
        for (k,v) in hop_symparams 
            for i in 1:size(v)[2]
                @show hop_symparams[k][:,i]
                @show etabar_sym[oindices2irreps[k]][i]
                hop_symparams[k][:,i] .*= etabar_sym[oindices2irreps[k]][i]
                @show hop_symparams[k][:,i]
                println()
            end
        end
        epsilon_symparams = Dict( k=>@.rescale(v,L,z,scale) for (k,v) in epsilon_symparams )
        u_symparams       = Dict( k=>@.rescale(v,L,z,scale) for (k,v) in u_symparams )
    else
        hop_symparams     = Dict( oirreps2indices[k]=>@.rescale(v,L,z,discretization;iterations=iterations,eta=x->1.0) for (k,v) in hop_symparams )
        epsilon_symparams = Dict( k=>@.rescale(v,L,z,discretization;iterations=iterations,eta=x->1.0) for (k,v) in epsilon_symparams )
        u_symparams       = Dict( k=>@.rescale(v,L,z,discretization;iterations=iterations,eta=x->1.0) for (k,v) in u_symparams )
    end
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
                   xi_symparams ,
                   mm_i ;
                   mine=mine ,
                   distributed=distributed ,
                   method=method ,
                   z=z ,
                   discretization=discretization ,
                   verbose=false ,
                   Nz=Nz)
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
            atom_config ,
            identityrep )

    atom_orbital_irreps = collect(keys(atom_config))

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

function atomic_spectrum( 
            label::String ,
            z::Float64 ,
            cg_o_dir::String ,
            asym_dir::String ,
            atom_config::Dict{String,Int64} ,
            identityrep::String ,
            epsilon_symparams::Dict{ Tuple{String,Int64} , ComplexF64 } ,
            u_symparams::Dict{ Tuple{String,Int64} , Matrix{ComplexF64} } ;
            max_spin2::Int64=10 ,
            distributed::Bool=false ,
            method::String="distfor" ,
            calculation::String="IMP" )

    # orbital irreps present in the atom
    atom_orbital_irreps::Vector{String} = collect(keys(atom_config))

    println( "====================" )
    println( "SETUP AND PARAMETERS" )
    println( "====================" )
    @show distributed 
    @show z
    distributed && @show distworkers
    distributed && @show method
    println( "OCCUPATION ENERGIES" )
    print_dict( epsilon_symparams ) 
    println( "COULOMB PARAMETERS" )
    print_dict( u_symparams ) 
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
    epsilon_symparams = Dict( k=>@.rescale(v,L,z,discretization) for (k,v) in epsilon_symparams )
    u_symparams       = Dict( k=>@.rescale(v,L,z,discretization) for (k,v) in u_symparams )
    println( "RESCALED PARAMETERS FOR H0" )
    @show epsilon_symparams 
    @show u_symparams 
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

end
