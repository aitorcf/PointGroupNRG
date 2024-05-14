# utilities 
spin_doubleint2index( s2::Float64 , S2::Float64 ) = Int64(0.5*(s2+S2)+1)
spin_index2doubleint( si::Int64 , S2::Int64 ) = Int64(2*(si-(S2+2)/2.0))

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
               if (occursin("_$(I_1)x$(I_2).",x) || occursin("_$(I_2)x$(I_1).",x)) ][1]
    verbose && @show file 

    inverted = (I_2!==I_1 && occursin("$(I_2)x$(I_1)",file))

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

function cg_shortcircuit_nonsimple( CG_PATH , oirreps ; verbose=false )

    seeds::Vector{String} = oirreps
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
    if verbose 
        @show produced 
        println()
    end

    if Set(produced)==Set(seeds)
        return produced 
    else 
        return cg_shortcircuit_nonsimple( CG_PATH , produced )
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
               ComplexF64[ get(cg_o_fulldict,(Is1,i1,Is2,i2,Is3,i3,r3),zero(ComplexF64))
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
    cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}} = get_cg_o_fullmatint_nonsimple( cg_o_full , oirreps )
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

function cg_reduce_product_states_nonsimple( symstates_1::Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , S } , 
                                             symstates_2::Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , S } ,
                                             cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,4} } ,
                                             cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,4} } ,
                                             oirreps2indices::Dict{String,Int64} )::Dict{ Tuple{Int64,String,Float64,Int64,Float64,Int64} , S } where {S<:State}
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

    oindices2irreps = [oirreps2indices[i] for i in 1:length(oirreps2indices)]

    spin_index2float( si::Int64 , S2::Int64 ) = si-(S2+2)/2.0

    irrep_combinations = Iterators.product( get_multiplets(symstates_1) ,
                                            get_multiplets(symstates_2) )

    symstates = Dict()
    for ircomb in irrep_combinations

        (N_1::Int64,I_1::String,S_1::Float64,r_1::Int64) = ircomb[1]
        (N_2::Int64,I_2::String,S_2::Float64,r_2::Int64) = ircomb[2]

        #redcomb::NTuple{2,Tuple{Int64,String,Float64}} = ( (N_1,I_1,S_1) , (N_2,I_2,S_2) )
        N_3::Int64 = N_1 + N_2
        #r_3::NTuple{2,Int64} = (r_1,r_2)

        #cg_comb::Tuple{Dict{ Tuple{String,Int64,String,Int64,String,Int64} , ComplexF64 },Dict{ NTuple{6,Float64} , ComplexF64 }} = 
        #        cg_combination( redcomb , oh_path )

        #cg_o::Dict{ Tuple{String,Int64,String,Int64,String,Int64} , ComplexF64 } = cg_comb[1]
        #cg_s::Dict{ NTuple{6,Float64} , ComplexF64 } = cg_comb[2]

        for ((I1,I2,I3),arr_orbital) in cg_o_fullmatint,
            ((S21,S22,S23),arr_spin) in cg_s_fullmatint

            (I1==oirreps2indices(I_1) && I2==oirreps2indices(I_2)) || continue
            (S21==S_1/2.0 && S22==S_2/2.0)                         || continue

            I_3 = oindices2irreps[I3]

            S_1 = S21/2.0
            S_2 = S22/2.0
            S_3 = S23/2.0

            @show I_1,I_2,I_3
            @show S_1,S_2,S_3

            for ro_3 in axes(arr_orbital,1),
                i_1 in axes(arr_orbital,2), 
                i_2 in axes(arr_orbital,3),
                i_3 in axes(arr_orbital,4),
                si_1 in axes(arr_spin,1),
                si_2 in axes(arr_spin,2),
                si_3 in axes(arr_spin,3)

                r_3 = (r_1,r_2,ro_3)

                s_1 = spin_index2float(si_1,S21)
                s_2 = spin_index2float(si_2,S22)
                s_3 = spin_index2float(si_3,S23)

                q_o = (I_1,i_1,I_2,i_2,I_3,i_3)
                q_s = (S_1,s_1,S_2,s_2,S_3,s_3)

                c_o = arr_orbital[ro_3,i_1,i_2,i_3]
                c_s = arr_spin[si_1,si_2,si_3]

                q_3::Tuple{ NTuple{3,Int64} , NTuple{3,String} , NTuple{3,Float64} , Int64 , Float64 , NTuple{3,Int64} } = 
                    ((N_1,N_2,N_3),(I_1,I_2,I_3),(S_1,S_2,S_3),i_3,s_3,r_3)

                s1::S = symstates_1[( N_1 , I_1 , S_1 , q_o[2] , q_s[2] , r_1 )]
                s2::S = symstates_2[( N_2 , I_2 , S_2 , q_o[4] , q_s[4] , r_2 )]
                merge!( + , symstates , Dict(q_3 => c_o*c_s*x(s1,s2)) )

            end


        end

        #for (q_o::Tuple{String,Int64,String,Int64,String,Int64},c_o::ComplexF64) in cg_o, 
        #    (q_s::NTuple{6,Float64},c_s::ComplexF64) in cg_s 

        #    I_3::String = q_o[end-1]
        #    S_3::Float64 = q_s[end-1]
        #    mu_3::Int64 = q_o[end]
        #    m_3::Float64 = q_s[end]

        #    q_3::Tuple{ NTuple{3,Int64} , NTuple{3,String} , NTuple{3,Float64} , Int64 , Float64 , NTuple{2,Int64} } = 
        #        ((N_1,N_2,N_3),(I_1,I_2,I_3),(S_1,S_2,S_3),mu_3,m_3,r_3)

        #    s1::S = symstates_1[( N_1 , I_1 , S_1 , q_o[2] , q_s[2] , r_1 )]
        #    s2::S = symstates_2[( N_2 , I_2 , S_2 , q_o[4] , q_s[4] , r_2 )]
        #    merge!( + , symstates , Dict(q_3 => c_o*c_s*x(s1,s2)) )
        #end
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
        S_3::Float64 = qnums[3][3]
        mu_3::Int64  = qnums[4]
        m_3::Float64 = qnums[5]

        p::Tuple{Int64,String,Float64,Int64,Float64} = ( N_3 , I_3 , S_3 , mu_3 , m_3 )
        merge!( + , principals , Dict(p => 1) )

        q_3::Tuple{Int64,String,Float64,Int64,Float64,Int64} = ( p... , principals[p] )

        get!( symstates_new , q_3 , state )
    end

    return symstates_new
end

# ---------------- #
# BASIS AND STATES #
# ---------------- #

function get_atomic_states_doublegroups( 
            external_label::Int64 ,
            ostring::String , 
            oirreps2dimensions::Dict{String,Int64} )::Vector{Tuple{Int64,String,Int64,String}}

    root = ( external_label , ostring )
    statetuples = Tuple{Int64,String,Int64,String}[]

    for orbital in 1:oirreps2dimensions[ostring]

        push!( statetuples , (root...,orbital,"-") )

    end

    return statetuples
end


function get_symstates_basis_multiplets_doublegroups_nonsimple( 
            atom_config::Dict{String,Int64},
            oirreps2dimensions::Dict{String,Int64} ,
            identityrep::String ,
            asym_dir::String ,
            cg_o_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,4} } ,
            cg_s_fullmatint::Dict{ NTuple{3,Int64} , Array{ComplexF64,3} } ,
            oirreps2indices::Dict{String,Int64} ;
            verbose=true )

    # hiztegia
    hiztegia::Dict{String,Any} = Dict( o=>o for (o,_) in atom_config )
    merge!( hiztegia , Dict( "-"=>0.0 ) )

    # symstates
    separate_symstates = []
    for (oirrep,multiplicity) in atom_config
        for m in 1:multiplicity
            onemult_states = get_atomic_states_doublegroups(m,oirrep,oirreps2dimensions)
            onemult_hilbert = HilbertSpace( onemult_states )
            onemult_symstates = oneirrep_symstates( 
                                    onemult_hilbert ,
                                    hiztegia ,
                                    identityrep ,
                                    "$(asym_dir)/$(oirrep)/" )
            push!( separate_symstates , onemult_symstates )
        end
    end
    symstates = separate_symstates[1]
    for symstates_new in separate_symstates[2:end]
        symstates = cg_reduce_product_states_nonsimple(
                        symstates , 
                        symstates_new ,
                        cg_o_fullmatint ,
                        cg_s_fullmatint ,
                        oirreps2indices )
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

function get_redmat_nonsimple( dictmat , 
                               multiplets_atom , 
                               multiplets_operator ,
                               cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}} , 
                               cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} ;
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

    # valid irrep combinations
    cg_o_combs = keys(cg_o_fullmatint)
    cg_s_combs = keys(cg_s_fullmatint)

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
        Dict( (G1,G2,G3)=>zeros(ComplexF64,size(cg_o_fullmatint[G2[2],G3[2],G1[2]],1),irr2mult_ij[G1],irr2mult_a[G2],irr2mult_ij[G3]) 
              for (G1,G2,G3) in Gcombs )

    # iterate over irreps (with multiplicities)
    @inbounds for (G_i,R_i) in irrmult_ij,
                  (G_a,R_a) in irrmult_a,
                  (G_j,R_j) in irrmult_ij

        ((G_i,G_a,G_j) in Gcombs) || continue

        # irrep quantum numbers
        (N_i,I_i,S_i) = G_i
        (N_j,I_j,S_j) = G_j
        (N_a,I_a,S_a) = G_a

        # filter irrep combination
        N_i==(N_j+N_a)                || continue 
        ((I_a,I_j,I_i) in cg_o_combs) || continue
        ((S_a,S_j,S_i) in cg_s_combs) || continue

        if verbose
            println( "IRREPS" )
            @show G_i,G_a,G_j
            println()
        end

        # clebsch-gordan matrices
        cgomat_aji = @view cg_o_fullmatint[I_a,I_j,I_i][:,:,:,:]
        cgsmat_aji = @view cg_s_fullmatint[S_a,S_j,S_i][:,:,:,:]

        # iterate over irrep multiplets 
        @inbounds for r_i in 1:R_i,
                      r_a in 1:R_a,
                      r_j in 1:R_j

            verbose && (@show r_i,r_a,r_j)

            for α in axes(cgomat_aji,1)

                verbose && (@show α)

                @views cgomat_aji_alpha = cgomat_aji[α,:,:,:]

                # Store  < m_i || f^+_{m_a} || m_j > = 0
                # without having to search in array
                matel = zero(ComplexF64)

                for i_a in axes(cgomat_aji_alpha,1),
                    i_j in axes(cgomat_aji_alpha,2),
                    si_a in axes(cgsmat_aji,1),
                    si_j in axes(cgsmat_aji,2)

                    # arbitrary
                    i_i = 1
                    si_i = 1

                    # full quantum numbers
                    q_i = (G_i...,i_i,spin_index2doubleint(si_i,S_i),r_i)
                    q_a = (G_a...,i_a,spin_index2doubleint(si_a,S_a),r_a)
                    q_j = (G_j...,i_j,spin_index2doubleint(si_j,S_j),r_j)

                    # lehmann amplitude
                    lehmann = get(dictmat,(q_i,q_a,q_j),zero(ComplexF64))
                    iszero(lehmann) && continue

                    # clebsch-gordan
                    cgo = cgomat_aji_alpha[i_a,i_j,i_i]
                    cgs = cgsmat_aji[si_a,si_j,si_i]
                    cg  = cgo*cgs

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
    end
    return redmat
end

function pcg_nonsimple_sanity_check( pcgdict::IntQPCG ,
                                     pcgred::IntIrrepPCGNS ,
                                     cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}} ,
                                     cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} ;
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

        Nu,Iu,Su,iu,su,ru = q_u
        Na,Ia,Sa,ia,sa,ra = q_a
        Nv,Iv,Sv,iv,sv,rv = q_v

        cg_o = cg_o_fullmatint[Ia,Iv,Iu]
        cg_s = cg_s_fullmatint[Sa,Sv,Su]

        sia = Int64(0.5*(sa+Sa)+1)
        siv = Int64(0.5*(sv+Sv)+1)
        siu = Int64(0.5*(su+Su)+1)

        Gu = (Nu,Iu,Su)
        Ga = (Na,Ia,Sa)
        Gv = (Nv,Iv,Sv)

        if verbose
            println("Dict matrix element:")
            @show q_u,q_a,q_v
            @show matel
            println()
            println( "------------" )
            println("Orbital Clebsch-Gordan coefficients")
            @show Gu,Ga,Gv
            @show Ia,Iu,Iv
            for α in axes(cg_o,1)

                println( "outer multipliticty α=$α" )

                for ia in axes(cg_o,2),
                    iv in axes(cg_o,3),
                    iu in axes(cg_o,4)

                    c = cg_o[α,ia,iv,iu]
                    println( "( $ia $iv | $iu ) = $c" )

                end
                println()
            end
            println( "------------" )
        end

        wignereckart = 0.0
        for α in axes(cg_o,1)

            red = pcgred[Gu,Ga,Gv][α,ru,ra,rv]
            verbose && println( "red = ",red," || ",ia,iv,iu," → " , cg_o[α,ia,iv,iu] )
            wignereckart += red*conj(cg_o[α,ia,iv,iu]*cg_s[sia,siv,siu])

        end

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

function get_pcgred_nonsimple( 
            basis::CB ,
            symstates_noint::ClearSymstateDict ,
            multiplets::IntMultipletSet,
            hiztegia::D ,
            oirreps2indices::Dict{String,Int64} ,
            cg_o_fullmatint::IntCGNS ,
            cg_s_fullmatint::IntCG ;
            verbose::Bool=false )::IntIrrepPCGNS where {CB<:CanonicalBasis,D<:Dict}

    # pcg in dict format
    pcg::IntQPCG = get_pseudoCG( symstates_noint , 
                                 basis , 
                                 hiztegia , 
                                 oirreps2indices )
    # pcg in reduced matrix format
    multiplets_a::IntMultipletVector = collect(filter( x->x[1]==1 , multiplets ))
    pcgred = get_redmat_nonsimple(
        pcg,
        multiplets ,
        multiplets_a ,
        cg_o_fullmatint ,
        cg_s_fullmatint ;
        verbose=true 
    )
    if verbose
        println( "PCGRED" )
        print_dict( pcgred )
        println()
    end

    pcg_nonsimple_sanity_check( pcg , pcgred , cg_o_fullmatint , cg_s_fullmatint ; verbose=false )

    return pcgred
end

function compute_D_cgsum_nonsimple(
            number_of_orbital_irreps::Int64 ,
            orbital_irreps_impurity_shell::Vector{Int64} ,
            orbital_irreps_one_electron::Vector{Int64} ,
            cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}} ,
            maximum_spin2::Int64 ,
            maximum_spin2_impurity::Int64 ,
            cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} ;
            verbose=false )

    # The D sum appears in the calculation of matrix elements between
    # new multiplets:
    #       ⟨u,αu||H||v,αv⟩ ~ ∑ ⟨i||a||j⟩ × ⟨ν||a||μ⟩^* × D + h.c.

    # ORBITAL SECTOR
    #
    maximum_orbital_multiplicity = maximum([
        maximum(size(cg_o_array)) 
        for cg_o_array in values(cg_o_fullmatint) 
    ])
    noi  = number_of_orbital_irreps
    noiis = maximum(orbital_irreps_impurity_shell)
    noio = maximum(orbital_irreps_one_electron)
    mom  = maximum_orbital_multiplicity
    D_orbital = zeros(ComplexF64,noi,noio,noi,noi,noiis,noiis,mom,mom,mom,mom) # Iu,Ia,Ii,Ij,Imu,Inu,αu,αv,αi,αj (Iv=Iu)
    for Iuv    in 1:number_of_orbital_irreps,
        Ia     in orbital_irreps_one_electron,
        Ii     in 1:number_of_orbital_irreps,
        Ij     in 1:number_of_orbital_irreps,
        Imu    in orbital_irreps_impurity_shell,
        Inu    in orbital_irreps_impurity_shell,
        αmuiuv in 1:maximum_orbital_multiplicity,
        αnujuv in 1:maximum_orbital_multiplicity,
        αaji   in 1:maximum_orbital_multiplicity,
        αamunu in 1:maximum_orbital_multiplicity

        element::ComplexF64 = zero(ComplexF64)

        haskey( cg_o_fullmatint , (Imu,Ii,Iuv) ) || continue
        haskey( cg_o_fullmatint , (Inu,Ij,Iuv) ) || continue
        haskey( cg_o_fullmatint , (Ia,Ij,Ii)   ) || continue
        haskey( cg_o_fullmatint , (Ia,Imu,Inu) ) || continue

        αmuiuv ≤ size(cg_o_fullmatint[Imu,Ii,Iuv],1)  || continue
        αnujuv ≤ size(cg_o_fullmatint[Inu,Ij,Iuv],1)  || continue
        αaji   ≤ size(cg_o_fullmatint[Ia,Ij,Ii],1)    || continue
        αamunu ≤ size( cg_o_fullmatint[Ia,Imu,Inu],1) || continue

        @views begin
            cg_o_muiuv = cg_o_fullmatint[Imu,Ii,Iuv][αmuiuv,:,:,:]
            cg_o_nujuv = cg_o_fullmatint[Inu,Ij,Iuv][αnujuv,:,:,:]
            cg_o_aji   = cg_o_fullmatint[Ia,Ij,Ii][αaji,:,:,:]
            cg_o_amunu = cg_o_fullmatint[Ia,Imu,Inu][αamunu,:,:,:]
        end

        iuv::Int64 = 1
        for ia in axes(cg_o_aji,1),
            ij in axes(cg_o_aji,2),
            ii in axes(cg_o_aji,3),
            imu in axes(cg_o_amunu,2),
            inu in axes(cg_o_amunu,3)

            element += conj(cg_o_muiuv[imu,ii,iuv])*cg_o_nujuv[inu,ij,iuv]*
                       conj(cg_o_aji[ia,ij,ii])*    cg_o_amunu[ia,imu,inu]

        end

        D_orbital[Iuv,Ia,Ii,Ij,Imu,Inu,αmuiuv,αnujuv,αaji,αamunu] = element

    end

    # SPIN SECTOR
    #
    maximum_spin2_onelectron = iszero(maximum_spin2) ? 0 : 1
    ms2 = maximum_spin2
    ms2o = maximum_spin2_onelectron
    ms2i = maximum_spin2_impurity
    D_spin = zeros(ComplexF64,ms2,ms2o,ms2,ms2,ms2i,ms2i) # Su,Sa,Si,Sj,Smu,Snu (Sv=Su)
    for Suv  in 1:number_of_orbital_irreps,
        Sa  in orbital_irreps_one_electron,
        Si  in 1:number_of_orbital_irreps,
        Sj  in 1:number_of_orbital_irreps,
        Smu in orbital_irreps_impurity_shell,
        Snu in orbital_irreps_impurity_shell

        element::ComplexF64 = zero(ComplexF64)

        haskey( cg_s_fullmatint , (Smu,Si,Suv) ) || continue
        haskey( cg_s_fullmatint , (Snu,Sj,Suv) ) || continue
        haskey( cg_s_fullmatint , (Sa,Sj,Si)   ) || continue
        haskey( cg_s_fullmatint , (Sa,Smu,Snu) ) || continue

        @views begin
            cg_s_muiuv = cg_s_fullmatint[Smu,Si,Suv][:,:,:]
            cg_s_nujuv = cg_s_fullmatint[Snu,Sj,Suv][:,:,:]
            cg_s_aji   = cg_s_fullmatint[Sa,Sj,Si][:,:,:]
            cg_s_amunu = cg_s_fullmatint[Sa,Smu,Snu][:,:,:]
        end

        siuv::Int64 = 1
        for sia  in axes(cg_s_aji,1),
            sij  in axes(cg_s_aji,2),
            sii  in axes(cg_s_aji,3),
            simu in axes(cg_s_amunu,2),
            sinu in axes(cg_s_amunu,3)

            element += conj(cg_s_muiuv[simu,sii,siuv])*cg_s_nujuv[sinu,sij,siuv]*
                       conj(cg_s_aji[sia,sij,sii])*    cg_s_amunu[sia,simu,sinu]

        end

        D_spin[Suv,Sa,Si,Sj,Smu,Snu] = element

    end

    return D_orbital,D_spin
end
function compute_F_cgsum_nonsimple(
            number_of_orbital_irreps::Int64 ,
            orbital_irreps_impurity_shell::Vector{Int64} ,
            orbital_irreps_one_electron::Vector{Int64} ,
            cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}} ,
            maximum_spin2::Int64 ,
            maximum_spin2_impurity::Int64 ,
            cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} ;
            verbose=false )

    # The K sum appears in the calculation of matrix elements
    # of block operators:
    #       ⟨u||A||v⟩_β = ∑_α ⟨i||A||j⟩_{α,n-1} × sgn(A,μ) × F(α,β) × δ_μν

    # ORBITAL SECTOR
    #
    maximum_orbital_multiplicity = maximum([
        maximum(size(cg_o_array)) 
        for cg_o_array in values(cg_o_fullmatint) 
    ])
    noi  = number_of_orbital_irreps
    noiis = maximum(orbital_irreps_impurity_shell)
    noio = maximum(orbital_irreps_one_electron)
    mom  = maximum_orbital_multiplicity
    F_orbital = zeros(ComplexF64,noi,noi,noio,noi,noiis,noiis,mom,mom,mom,mom) # Iu,Iv,Ia,Iij,Imu,Inu,αuijmu,αvijnu,αmuanu,αuav (Iv=Iu)
    for Iu     in 1:number_of_orbital_irreps,
        Iv     in 1:number_of_orbital_irreps,
        Ia     in orbital_irreps_one_electron,
        Iij    in 1:number_of_orbital_irreps,
        Imu    in orbital_irreps_impurity_shell,
        Inu    in orbital_irreps_impurity_shell,
        αmuiju in 1:maximum_orbital_multiplicity,
        αnuijv in 1:maximum_orbital_multiplicity,
        αanumu in 1:maximum_orbital_multiplicity,
        αavu   in 1:maximum_orbital_multiplicity

        element::ComplexF64 = zero(ComplexF64)

        haskey( cg_o_fullmatint , (Imu,Iij,Iu) ) || continue
        haskey( cg_o_fullmatint , (Inu,Iij,Iv) ) || continue
        haskey( cg_o_fullmatint , (Ia,Inu,Imu) ) || continue
        haskey( cg_o_fullmatint , (Ia,Iv,Iu)   ) || continue

        αmuiju ≤ size(cg_o_fullmatint[Imu,Iij,Iu],1) || continue
        αnuijv ≤ size(cg_o_fullmatint[Inu,Iij,Iu],1) || continue
        αanumu ≤ size(cg_o_fullmatint[Ia,Inu,Imu],1) || continue
        αavu   ≤ size(cg_o_fullmatint[Ia,Iv,Iu],1)   || continue

        @views begin
            cg_o_muiju = cg_o_fullmatint[Imu,Iij,Iu][αmuiju,:,:,:]
            cg_o_nuijv = cg_o_fullmatint[Inu,Iij,Iu][αnuijv,:,:,:]
            cg_o_anumu = cg_o_fullmatint[Ia,Inu,Imu][αanumu,:,:,:]
            cg_o_avu   = cg_o_fullmatint[Ia,Iv,Iu][αavu,:,:,:]
        end

        iu::Int64 = 1
        for iij in axes(cg_o_muiju,2),
            imu in axes(cg_o_muiju,1),
            inu in axes(cg_o_anumu,2),
            ia  in axes(cg_o_anumu,1),
            iv  in axes(cg_o_avu,2)

            element += conj(cg_o_muiju[imu,iij,iu])*cg_o_nuijv[inu,iij,iv]*
                       conj(cg_o_anumu[ia,inu,imu])*cg_o_avu[ia,iv,iu]

        end

        K_orbital[Iu,Iv,Ia,Iij,Imu,Inu,αmuiju,αnuijv,αanumu,αavu] = element

    end

    # SPIN SECTOR
    #
    maximum_spin2_onelectron = iszero(maximum_spin2) ? 0 : 1
    ms2 = maximum_spin2
    ms2o = maximum_spin2_onelectron
    ms2i = maximum_spin2_impurity
    K_spin = zeros(ComplexF64,ms2,ms2,ms2o,ms2,ms2i,ms2i) # Su,Sv,Sa,Sij,Smu,Snu (Si=Sj)
    for Su  in 1:number_of_orbital_irreps,
        Sv  in 1:number_of_orbital_irreps,
        Sa  in orbital_irreps_one_electron,
        Sij  in 1:number_of_orbital_irreps,
        Smu in orbital_irreps_impurity_shell,
        Snu in orbital_irreps_impurity_shell

        element::ComplexF64 = zero(ComplexF64)

        haskey( cg_s_fullmatint , (Smu,Sij,Su) ) || continue
        haskey( cg_s_fullmatint , (Snu,Sij,Sv) ) || continue
        haskey( cg_s_fullmatint , (Sa,Snu,Smu) ) || continue
        haskey( cg_s_fullmatint , (Sa,Sv,Su)   ) || continue

        @views begin
            cg_s_muiju = cg_s_fullmatint[Smu,Sij,Su][:,:,:]
            cg_s_nuijv = cg_s_fullmatint[Snu,Sij,Sv][:,:,:]
            cg_s_anumu = cg_s_fullmatint[Sa,Snu,Smu][:,:,:]
            cg_s_avu   = cg_s_fullmatint[Sa,Sv,Su][:,:,:]
        end

        siu::Int64 = 1
        for siij  in axes(cg_s_muiju,2),
            simu  in axes(cg_s_muiju,1),
            sinu  in axes(cg_s_nuijv,1),
            sia   in axes(cg_s_avu,1),
            siv   in axes(cg_s_avu,2)

            element += conj(cg_s_muiju[simu,siij,siu])*cg_s_nuijv[sinu,siij,siv]*
                       conj(cg_s_anumu[sia,sinu,simu])*cg_s_avu[sia,siv,siu]

        end

        K_spin[Su,Sv,Sa,Sij,Smu,Snu] = element

    end

    return F_orbital,F_spin
end
function compute_K_cgsum_nonsimple(
            number_of_orbital_irreps::Int64 ,
            orbital_irreps_impurity_shell::Vector{Int64} ,
            orbital_irreps_one_electron::Vector{Int64} ,
            cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}} ,
            maximum_spin2::Int64 ,
            maximum_spin2_impurity::Int64 ,
            cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} ;
            verbose=false )

    # The K sum appears in the calculation of matrix elements
    # of local operators
    #       ⟨u||A||v⟩_{n,β} = ∑_α ⟨μ||A||ν⟩_{n,α} × K(α,β) × δ_ij 

    # ORBITAL SECTOR
    #
    maximum_orbital_multiplicity = maximum([
        maximum(size(cg_o_array)) 
        for cg_o_array in values(cg_o_fullmatint) 
    ])
    noi  = number_of_orbital_irreps
    noiis = maximum(orbital_irreps_impurity_shell)
    noio = maximum(orbital_irreps_one_electron)
    mom  = maximum_orbital_multiplicity
    K_orbital = zeros(ComplexF64,noi,noi,noio,noi,noiis,noiis,mom,mom,mom,mom) # Iu,Iv,Ia,Iij,Imu,Inu,αuijmu,αvijnu,αmuanu,αuav (Iv=Iu)
    for Iu     in 1:number_of_orbital_irreps,
        Iv     in 1:number_of_orbital_irreps,
        Ia     in orbital_irreps_one_electron,
        Iij    in 1:number_of_orbital_irreps,
        Imu    in orbital_irreps_impurity_shell,
        Inu    in orbital_irreps_impurity_shell,
        αmuiju in 1:maximum_orbital_multiplicity,
        αnuijv in 1:maximum_orbital_multiplicity,
        αanumu in 1:maximum_orbital_multiplicity,
        αavu   in 1:maximum_orbital_multiplicity

        element::ComplexF64 = zero(ComplexF64)

        haskey( cg_o_fullmatint , (Imu,Iij,Iu) ) || continue
        haskey( cg_o_fullmatint , (Inu,Iij,Iv) ) || continue
        haskey( cg_o_fullmatint , (Ia,Inu,Imu) ) || continue
        haskey( cg_o_fullmatint , (Ia,Iv,Iu)   ) || continue

        αmuiju ≤ size(cg_o_fullmatint[Imu,Iij,Iu],1) || continue
        αnuijv ≤ size(cg_o_fullmatint[Inu,Iij,Iu],1) || continue
        αanumu ≤ size(cg_o_fullmatint[Ia,Inu,Imu],1) || continue
        αavu   ≤ size(cg_o_fullmatint[Ia,Iv,Iu],1)   || continue

        @views begin
            cg_o_muiju = cg_o_fullmatint[Imu,Iij,Iu][αmuiju,:,:,:]
            cg_o_nuijv = cg_o_fullmatint[Inu,Iij,Iu][αnuijv,:,:,:]
            cg_o_anumu = cg_o_fullmatint[Ia,Inu,Imu][αanumu,:,:,:]
            cg_o_avu   = cg_o_fullmatint[Ia,Iv,Iu][αavu,:,:,:]
        end

        iu::Int64 = 1
        for iij in axes(cg_o_muiju,2),
            imu in axes(cg_o_muiju,1),
            inu in axes(cg_o_anumu,2),
            ia  in axes(cg_o_anumu,1),
            iv  in axes(cg_o_avu,2)

            element += conj(cg_o_muiju[imu,iij,iu])*cg_o_nuijv[inu,iij,iv]*
                       conj(cg_o_anumu[ia,inu,imu])*cg_o_avu[ia,iv,iu]

        end

        K_orbital[Iu,Iv,Ia,Iij,Imu,Inu,αmuiju,αnuijv,αanumu,αavu] = element

    end

    # SPIN SECTOR
    #
    maximum_spin2_onelectron = iszero(maximum_spin2) ? 0 : 1
    ms2 = maximum_spin2
    ms2o = maximum_spin2_onelectron
    ms2i = maximum_spin2_impurity
    K_spin = zeros(ComplexF64,ms2,ms2,ms2o,ms2,ms2i,ms2i) # Su,Sv,Sa,Sij,Smu,Snu (Si=Sj)
    for Su  in 1:number_of_orbital_irreps,
        Sv  in 1:number_of_orbital_irreps,
        Sa  in orbital_irreps_one_electron,
        Sij  in 1:number_of_orbital_irreps,
        Smu in orbital_irreps_impurity_shell,
        Snu in orbital_irreps_impurity_shell

        element::ComplexF64 = zero(ComplexF64)

        haskey( cg_s_fullmatint , (Smu,Sij,Su) ) || continue
        haskey( cg_s_fullmatint , (Snu,Sij,Sv) ) || continue
        haskey( cg_s_fullmatint , (Sa,Snu,Smu) ) || continue
        haskey( cg_s_fullmatint , (Sa,Sv,Su)   ) || continue

        @views begin
            cg_s_muiju = cg_s_fullmatint[Smu,Sij,Su][:,:,:]
            cg_s_nuijv = cg_s_fullmatint[Snu,Sij,Sv][:,:,:]
            cg_s_anumu = cg_s_fullmatint[Sa,Snu,Smu][:,:,:]
            cg_s_avu   = cg_s_fullmatint[Sa,Sv,Su][:,:,:]
        end

        siu::Int64 = 1
        for siij  in axes(cg_s_muiju,2),
            simu  in axes(cg_s_muiju,1),
            sinu  in axes(cg_s_nuijv,1),
            sia   in axes(cg_s_avu,1),
            siv   in axes(cg_s_avu,2)

            element += conj(cg_s_muiju[simu,siij,siu])*cg_s_nuijv[sinu,siij,siv]*
                       conj(cg_s_anumu[sia,sinu,simu])*cg_s_avu[sia,siv,siu]

        end

        K_spin[Su,Sv,Sa,Sij,Smu,Snu] = element

    end

    return K_orbital,K_spin
end

# μ,i→uα : [μ,i,u,α] for n=0
function get_combinations_Gu_muiualpha_initial( 
            identityrep::String ,
            calculation::String ,
            multiplets_0 ,
            oirreps2indices::Dict{String,Int64} ;
            verbose=true )

    # Combinations as a vector
    combinations_muiualpha::Vector{Tuple{IntMultiplet,IntMultiplet,IntMultiplet,Int64}} = []

    m_vac = (0,oirreps2indices[identityrep],0,1)
    if calculation=="IMP"
        for m_mu in multiplets_0
            push!( combinations_muiualpha , (m_mu,m_vac,m_mu,1) )
        end 
    elseif calculation=="CLEAN" 
        push!( combinations_muiualpha , (m_vac,m_vac,m_vac,1) )
    end

    irreps_u = Set( u[1:3] for (_,_,u,_) in combinations_muiualpha )
    combinations_Gu_muiualpha= Dict{ IntIrrep , Vector{Tuple{IntMultiplet,IntMultiplet,IntMultiplet,Int64}} }(
        G => [ 
            (m_mu,m_i,m_u,alpha)
            for (m_mu,m_i,m_u,alpha) in combinations_muiualpha
            if m_u[1:3]==G
        ] for G in irreps_u
    ) 

    if verbose 
        println( "COMBINATIONS μ,i→u FOR N=0" )
        for (G,combs) in combinations_Gu_muiualpha
            println( "$G => $combs" )
        end
        println()
    end

    return combinations_Gu_muiualpha
end
function get_combinations_Gu_muiualpha( 
            multiplets_block::Set{NTuple{4,Int64}} , 
            multiplets_shell::Set{NTuple{4,Int64}} ,
            keys_as_dict_o::Dict{NTuple{2,Int64},Vector{NTuple{2,Int64}}} ,
            keys_as_dict_s::Dict{NTuple{2,Int64},Vector{Int64}} ;
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
        IIAA_u::Vector{NTuple{2,Int64}} = keys_as_dict_o[(I_mu,I_i)]
        SS_u::Vector{Int64} = keys_as_dict_s[(S_mu,S_i)]

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

function construct_and_diagonalize_uHv( 
        multiplets_block::Set{NTuple{4,Int64}} , 
        multiplets_shell::Set{NTuple{4,Int64}} ,
        irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
        hop_symparams::Dict{ Int64 , Matrix{ComplexF64} } ,
        keys_as_dict_o::Dict{NTuple{2,Int64},Vector{NTuple{2,Int64}}} ,
        keys_as_dict_s::Dict{NTuple{2,Int64},Vector{Int64}} ,
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
    combinations_Gu_muiualpha_new = get_combinations_Gu_muiualpha( 
        multiplets_block ,
        multiplets_shell ,
        keys_as_dict_o ,
        keys_as_dict_s
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
                                hoparam_matrix = hop_symparams[G_a[2]]

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

function compute_lehmann_iaj(
        multiplets_a_block::Vector{NTuple{4,Int64}}, 
        ksum::KSum ,
        lehmann_muanu::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,4} },
        combinations_Gu_muiualpha::Dict{ IntIrrep , Vector{Tuple{IntMultiplet,IntMultiplet,IntMultiplet,Int64}}} ,
        irrEU::Dict{ IntIrrep , Tuple{Vector{Float64},Matrix{ComplexF64}} } ,
        cg_o_comb2A::Dict{ NTuple{3,Int64} , Int64 } # μ,i,u => n_u (A)
    )::Dict{NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,4} } 


    # G => multiplicity of G
    #
    # i (n) / u (n-1)
    G2R_uv_precutoff::Dict{IntIrrep,Int64} = Dict( G_u=>size(U_u,1) for (G_u,(_,U_u)) in irrEU )
    G2R_uv_postcutoff::Dict{IntIrrep,Int64} = Dict( G_u=>length(E_u) for (G_u,(E_u,_)) in irrEU )
    # a
    G2R_a::Dict{NTuple{3,Int64},Int64} = 
        Dict( G=>R for (G,R) in get_irreps( Set(multiplets_a_block) ; multiplicity=true ) )

    # initialize result dictionary
    lehmann_iaj::Dict{NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,4} } = Dict()

    # irrep -> decomposition -> multiplet
    Gu2GmuGi2rmuiualpha::Dict{IntIrrep,Dict{NTuple{2,IntIrrep},Vector{NTuple{4,Int64}}}} = Dict(
        G_u => Dict( (G_mu,G_i) => [
                    (m_mu[4],m_i[4],m_u[4],α_u) 
                    for (m_mu,m_i,m_u,α_u) in muiualphas 
                    if (m_i[1:3]==G_i && m_mu[1:3]==G_mu)
               ] for (G_mu,G_i) in Set( (m_mu[1:3],m_i[1:3]) for (m_mu,m_i,_,_) in muiualphas)
        )
        for (G_u,muiualphas) in combinations_Gu_muiualpha if G_u in keys(irrEU)
    )

    # full-sized matrices for avoiding allocations
    R_uv_precutoff_max::Int64 = maximum(values(G2R_uv_precutoff))
    R_uv_postcutoff_max::Int64 = maximum(values(G2R_uv_postcutoff))
    R_a_max::Int64 = maximum(values(G2R_a))
    A_max::Int64 = maximum(values(cg_o_comb2A))
    uav_matrix_full::Array{ComplexF64,4} = zeros(ComplexF64,A_max,R_uv_precutoff_max,R_a_max,R_uv_precutoff_max)
    tmp_full::Matrix{ComplexF64} = zeros(ComplexF64,R_uv_precutoff_max,R_uv_postcutoff_max)

    # G_iu, G_jv iteration
    # i',j' for block states from previous step
    for (G_iu,GmuGip2rmuiualpha) in Gu2GmuGi2rmuiualpha,
        (G_jv,GnuGjp2rnujvalpha) in Gu2GmuGi2rmuiualpha

        # irrep quantum numbers 
        N_iu,I_iu,S_iu = G_iu
        N_jv,I_jv,S_jv = G_jv

        # early discard 
        N_iu==(N_jv+1) || continue

        # multiplicities
        R_u_precutoff = G2R_uv_precutoff[G_iu]
        R_v_precutoff = G2R_uv_precutoff[G_jv]
        R_u_postcutoff = G2R_uv_postcutoff[G_iu]
        R_v_postcutoff = G2R_uv_postcutoff[G_jv]

        # transformation matrices
        @views begin
            U_iu = irrEU[G_iu][2][1:R_u_precutoff,1:R_u_postcutoff] 
            U_jv = irrEU[G_jv][2][1:R_v_precutoff,1:R_v_postcutoff] 
            tmp = tmp_full[1:R_u_precutoff,1:R_v_postcutoff]
        end

        # G_a iteration
        for (G_a,R_a) in G2R_a

            # irrep quantum numbers
            N_a,I_a,S_a = G_a

            # < Γu || f†_Γa || Γv >
            B = get( cg_o_comb2A , (I_a,I_jv,I_iu) , 0 )
            iszero(B) && continue
            @views uav_matrix = uav_matrix_full[1:B,1:R_u_precutoff,1:R_a,1:R_v_precutoff]
            uav_matrix .= zero(ComplexF64)

            # G_i',G_mu,G_j',G_nu iteration
            for ((G_mu,G_ip),rmuiualpha) in GmuGip2rmuiualpha,
                ((G_nu,G_jp),rnujvalpha) in GnuGjp2rnujvalpha

                # irrep quantum numbers 
                G_ipjp = G_ip
                N_ipjp,I_ipjp,S_ipjp = G_ipjp
                N_mu,I_mu,S_mu = G_mu
                N_nu,I_nu,S_nu = G_nu

                # early discard
                G_ip==G_jp     || continue
                N_mu==(N_nu+1) || continue

                # clebsch-gordan sum
                k_spin,k_orbital_array = ksum[(G_iu,G_jv,G_a,G_ipjp,G_mu,G_nu)]
                iszero(k_spin) && continue
                iszero(length(k_orbital_array)) && continue

                # <Γμ||f^†_Γa||Γν>_α , α=1,…,A
                lehmann_muanu_array = get( lehmann_muanu , (G_mu,G_a,G_nu) , zeros(ComplexF64,0,0,0,0) )
                iszero(length(lehmann_muanu_array)) && continue
                A = size(lehmann_muanu_array,1)

                # u,i,mu , v,j,nu  multiplet iteration
                for (r_mu,r_i,r_u,αu) in rmuiualpha,
                    (r_nu,r_j,r_v,αv) in rnujvalpha

                    r_i==r_j || continue

                    for α in 1:A,
                        β in 1:B

                        k = k_spin*k_orbital_array[αu,αv,α,β]

                        for r_a in 1:R_a
                            uav_matrix[β,r_u,r_a,r_v] += k * lehmann_muanu_array[α,r_mu,r_a,r_nu]
                        end # r_a iteration

                    end # α,β iteration
                end # u,i,mu , v,j,nu , a multiplet iteration
            end # G_i,G_mu,G_j,G_nu iteration

            # final discard
            #isapprox(sum(abs2.(uav_matrix)),zero(ComplexF64)) && continue

            # transform matrix
            #lehmann_iaj[G_iu,G_a,G_jv] = zeros(ComplexF64,B,R_u_postcutoff,R_a,R_v_postcutoff)
            lehmann_matrix = zeros(ComplexF64,B,R_u_postcutoff,R_a,R_v_postcutoff)
            for r_a in 1:R_a,
                β in 1:B

                #@views begin
                #    uav_matrix_view = uav_matrix[β,:,r_a,:]
                #    lehmann_matrix_view = lehmann_iaj[G_iu,G_a,G_jv][β,:,r_a,:]
                #end

                #lehmann_matrix_view .= (U_iu' * uav_matrix_view * U_jv)
                #@inbounds begin
                #    mul!( tmp , uav_matrix_view , U_jv )
                #    mul!( lehmann_matrix_view , U_iu' , tmp )
                #end
                @views begin
                    mul!( tmp , uav_matrix[β,:,r_a,:] , U_jv )
                    mul!( lehmann_matrix[β,:,r_a,:] , U_iu' , tmp )
                end

            end # r_a,α,β iteration
            lehmann_iaj[G_iu,G_a,G_jv] = copy(lehmann_matrix)

        end # G_a iteration
    end # G_u, G_v iteration

    return lehmann_iaj
end

function cut_off_nonsimple!( 
            irrEU::Dict{ Tuple{Int64,Int64,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} };
            type::String="multiplet" , 
            cutoff::T=200 , 
            safeguard::Bool=true , 
            safeguard_tol::Float64=1e-5 ,
            safeguard_max::Int64=200 ,
            minmult::Int64=0 , 
            mine::Float64=-1.0 ,
            verbose::Bool=true ,
            impurity_operator::Dict{ IntTripleG , Array{ComplexF64,4} }=Dict{ IntTripleG , Array{ComplexF64,4} }()
    ) where {T<:Number}
    # input 
    # - irrEU : G => (E,U)
    # - type : "multiplet" --> keep a limited number of lowest-energy multiplets
    #          "energy"    --> upper bound to energy
    # - cutoff : upper bound in multiplets/energy
    # output: 
    # - [irrEU : altered, not returned]
    # - { kept } : kept block multiplets
    # - discarded : amount of multiplets discarded

    # order multiplets in terms of their energy 
    mm::Vector{Tuple{IntMultiplet,Float64}} = collect( ((G...,r),E[r]) for (G,(E,U)) in irrEU for r=1:length(E) )
    sort!( mm , by=x->x[2] )

    # discard, kept and safeguarded
    discarded::Vector{IntMultiplet} = IntMultiplet[]
    kept::Vector{IntMultiplet} = IntMultiplet[]
    sg::Vector{IntMultiplet} = IntMultiplet[]

    # multiplet cutoff
    if type=="multiplet"

        # excessive cutoff, no discard
        if cutoff>=length(mm)

            kept = map( x->x[1] , mm )
            discarded = []

        # discard needed
        else

            mine_idx::Int64 = mm[end][2]<mine ? cutoff : findfirst( x->x[2]>=mine , mm )
            cutoff = maximum([ mine_idx , cutoff ])
            cutoff_energy = mm[cutoff][2]
            #kept = map( x->x[1] , mm[1:cutoff] )

            mm_kept = filter( x->(x[2]<=cutoff_energy || (safeguard && isapprox(x[2],cutoff_energy;atol=safeguard_tol))) , mm ) 
            length(mm_kept)>(cutoff+safeguard_max) && (mm_kept=mm_kept[1:(cutoff+safeguard_max)])
            mm_discarded = mm[(length(mm_kept)+1):end]

            kept      = map( x->x[1] , mm_kept )
            discarded = map( x->x[1] , mm_discarded )
        end

    # energy cutoff
    elseif type=="energy" 
        kept      = collect( pair[1] for pair in mm if pair[2]<=cutoff  )
        cutoff    = length(kept)
        (cutoff<minmult && minmult<length(mm)) && (cutoff=minmult)
        kept      = collect( pair[1] for pair in mm[1:cutoff] )
        if cutoff<length(mm) 
            if safeguard
                safeguard = map( x->x[1] ,
                                [ m for m in mm[(cutoff+1):end] 
                                  if isapprox(m[2],mm[cutoff][2];atol=1e-1)] )
                append!( kept , safeguard )
            end
            cutoff = length(kept)
            discarded = map( x->x[1] , mm[(cutoff+1):end] )
        else 
            discarded = []
        end
    else 
        error( "type must be 'multiplet' or 'energy'" )
    end

    # apply cutoff to irrEU
    for (G,(E,U)) in irrEU 
        N = length(collect( k for k in kept if k[1:3]==G ))
        if N==0 
            pop!( irrEU , G )
            continue
        end
        # U does not need cutoff (multiplets run over)
        irrEU[G] = ( E[1:N] , U[:,1:N] )
    end

    # cut off impurity operator
    for ((G_u,G_a,G_v),uavmat) in impurity_operator
        if !(haskey(irrEU,G_u) && haskey(irrEU,G_v)) 
            pop!( impurity_operator , (G_u,G_a,G_v) )
            continue
        end
        R_u = length(irrEU[G_u][1]) 
        R_v = length(irrEU[G_v][1]) 
        impurity_operator[G_u,G_a,G_v] = uavmat[:,1:R_u,:,1:R_v]
    end

    return ( Set(kept) , discarded )::Tuple{Set{NTuple{4,Int64}},Vector{NTuple{4,Int64}}}
end

function NRG_doublegroups_nonsimple(
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
              impurity_operators::Dict{String,Dict{IntTripleG,Array{ComplexF64,4}}}=Dict{String,Dict{IntTripleG,Array{ComplexF64,3}}}() ,
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
        (multiplets_block, discarded) = cut_off_nonsimple!( 
            irrEU ; 
            type=cutoff_type , 
            cutoff=cutoff_magnitude , 
            safeguard=true ,
            minmult=minmult ,
            mine=mine ,
            verbose=false ,
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
            I_o => diagm(ComplexF64.(channels_codiagonals[n-1][I_o])) 
            for I_o in keys(channels_codiagonals[1])
        )
        println( "CONDUCTION COUPLING TERMS" )
        println( "Codiagonals (hoppings)")
        @printf "  %-8s  %-10s  %-s\n" "orbital" "multiplet" "amplitude"
        for (I,hop_matrix) in hop_symparams,
            cartesian_index in CartesianIndices(hop_matrix)

            i,j = Tuple(cartesian_index)
            isapprox(hop_matrix[i,j],0.0) && continue

            @printf "  %-8s  %-3i => %-3i  %-.3f\n" I i j hop_matrix[i,j]

        end
        println( "Diagonals (occupation energies)")
        @printf "  %-8s  %-10s  %-s\n" "irrep" "multiplet" "amplitude"
        for (I_o,I_o_multiplets_diagonals) in channels_diagonals[n],
            (r_o,r_o_multiplet_diagonal) in enumerate(I_o_multiplets_diagonals)

            @printf "  %-8s  %-3i => %-.3f\n" I_o r_o r_o_multiplet_diagonal

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
        diagonalization_performance = @timed (irrEU,combinations_Gu_muiualpha) = construct_and_diagonalize_uHv(
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

function nrg_full_doublegroup( 
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
     cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}}) = get_cg_o_info_nonsimple( cg_o_dir , atom_orbital_irreps )


    # for dmnrg
    shell_dimension = reduce( * , [4^(oirreps2dimensions[I]*R) for (I,R) in shell_config] )

    # spin symmetry
    cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} = Dict( (0,0,0)=>[1.0;;;] )

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
            get_symstates_basis_multiplets_doublegroups_nonsimple( 
                    impurity_config,
                    oirreps2dimensions,
                    identityrep,
                    multiplets_dir,
                    cg_o_fullmatint ,
                    cg_s_fullmatint ,
                    oirreps2indices ;
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
    lehmann_iaj::IntIrrepPCGNS = 
        calculation=="CLEAN" ? 
        IntIrrepPCG() :
        get_pcgred_nonsimple( basis_atom ,
                    symstates_atom_noint::ClearSymstateDict ,
                    multiplets_atom::IntMultipletSet ,
                    hiztegia ,
                    oirreps2indices::Dict{String,Int64} ,
                    cg_o_fullmatint::IntCGNS ,
                    cg_s_fullmatint::IntCG ;
                    verbose=false )::IntIrrepPCGNS

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
        get_symstates_basis_multiplets_doublegroups_nonsimple(
                shell_config,
                oirreps2dimensions,
                identityrep,
                multiplets_dir,
                cg_o_fullmatint ,
                cg_s_fullmatint ,
                oirreps2indices ;
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
    lehmann_muanu::Dict{NTuple{3,NTuple{3,Int64}},Array{ComplexF64,4}} = get_pcgred_nonsimple( 
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

    #   ------------------------------ #
    #%% precompute clebsch-gordan sums #
    #   ------------------------------ #
    #
    dsum,ksum,fsum = ClebschGordanSums(
        length(oindices2irreps),
        collect(union(Set(x[2] for x in multiplets_atom),Set(x[2] for x in multiplets_shell))) ,
        collect(Set(x[2] for x in multiplets_a_shell if x[1]==1)) ,
        cg_o_fullmatint,
        max_spin2 ,
        maximum(vcat([ x[3] for x in multiplets_atom ],[x[3] for x in multiplets_shell])) ,
        cg_s_fullmatint ;
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
    (irrEU,combinations_Gu_muiualpha) = construct_and_diagonalize_uHv( 
                    multiplets_atom , 
                    multiplets_shell ,
                    irrEU , 
                    hop_symparams_int , 
                    keys_as_dict_o ,
                    keys_as_dict_s ,
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

        nrg = NRG_doublegroups_nonsimple( 
            label,
            calculation,
            iterations, 
            cutoff_type, 
            cutoff_magnitude,
            L,
            irrEU,
            multiplets_shell, 
            keys_as_dict_o,
            keys_as_dict_s,
            cg_o_comb2A,
            dsum,
            ksum,
            lehmann_muanu,
            collect(multiplets_a_shell) ,
            combinations_Gu_muiualpha,
            betabar,
            oindex2dimensions,
            channels_codiagonals,
            max_spin2;
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

function nrg_full_doublegroup_nonsimple( 
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
     cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}}) = get_cg_o_info_nonsimple( cg_o_dir , atom_orbital_irreps )


    # for dmnrg
    shell_dimension = reduce( * , [4^(oirreps2dimensions[I]*R) for (I,R) in shell_config] )

    # spin symmetry
    cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} = Dict( (0,0,0)=>[1.0;;;] )

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
            get_symstates_basis_multiplets_doublegroups_nonsimple( 
                    impurity_config,
                    oirreps2dimensions,
                    identityrep,
                    multiplets_dir,
                    cg_o_fullmatint ,
                    cg_s_fullmatint ,
                    oirreps2indices ;
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
    lehmann_iaj::IntIrrepPCGNS = 
        calculation=="CLEAN" ? 
        IntIrrepPCG() :
        get_pcgred_nonsimple( basis_atom ,
                    symstates_atom_noint::ClearSymstateDict ,
                    multiplets_atom::IntMultipletSet ,
                    hiztegia ,
                    oirreps2indices::Dict{String,Int64} ,
                    cg_o_fullmatint::IntCGNS ,
                    cg_s_fullmatint::IntCG ;
                    verbose=false )::IntIrrepPCGNS

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
        get_symstates_basis_multiplets_doublegroups_nonsimple(
                shell_config,
                oirreps2dimensions,
                identityrep,
                multiplets_dir,
                cg_o_fullmatint ,
                cg_s_fullmatint ,
                oirreps2indices ;
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
    lehmann_muanu::Dict{NTuple{3,NTuple{3,Int64}},Array{ComplexF64,4}} = get_pcgred_nonsimple( 
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

    #   ------------------------------ #
    #%% precompute clebsch-gordan sums #
    #   ------------------------------ #
    #
    dsum,ksum,fsum = ClebschGordanSums(
        length(oindices2irreps),
        collect(union(Set(x[2] for x in multiplets_atom),Set(x[2] for x in multiplets_shell))) ,
        collect(Set(x[2] for x in multiplets_a_shell if x[1]==1)) ,
        cg_o_fullmatint,
        max_spin2 ,
        maximum(vcat([ x[3] for x in multiplets_atom ],[x[3] for x in multiplets_shell])) ,
        cg_s_fullmatint ;
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
    (irrEU,combinations_Gu_muiualpha) = construct_and_diagonalize_uHv( 
                    multiplets_atom , 
                    multiplets_shell ,
                    irrEU , 
                    hop_symparams_int , 
                    keys_as_dict_o ,
                    keys_as_dict_s ,
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

        nrg = NRG_doublegroups_nonsimple( 
            label,
            calculation,
            iterations, 
            cutoff_type, 
            cutoff_magnitude,
            L,
            irrEU,
            multiplets_shell, 
            keys_as_dict_o,
            keys_as_dict_s,
            cg_o_comb2A,
            dsum,
            ksum,
            lehmann_muanu,
            collect(multiplets_a_shell) ,
            combinations_Gu_muiualpha,
            betabar,
            oindex2dimensions,
            channels_codiagonals,
            max_spin2;
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
