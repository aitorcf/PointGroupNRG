function compute_multiplets( symmetry::String ; # "PS": Pointgroup-Spin, "D": Doublegroup, "J": Total angular momentum
                             multiplets_path::String="multiplets" , # where to store the output
                             irrep::SF="" ,
                             clebschgordan_path::String="" ,
                             verbose::Bool=false ) where {SF<:Union{String,Float64}}

    # add convention name to multiplet folder
    if multiplets_path[end]=="/"
        multiplets_path = multiplets_path[1:(end-1)] 
    end
    multiplets_path *= "/$(irrep)"

    # create asym dir if it does not exist
    isdir(multiplets_path) || mkpath( multiplets_path )

    nn_irrep = istotalangularmomentum(symmetry) ? (2:Int64(2*irrep+1)) : (2:get_N_nonsimple(irrep,clebschgordan_path,symmetry))

    # compute multiplet states
    for n in nn_irrep

        compute_asymstates_N_allsymmetries( 
            symmetry::String ,
            irrep::SF , 
            n::Int64 ,
            multiplets_path::String ;
            cg_path=clebschgordan_path ,
            verbose=verbose ,
        )

        #if simple
        #    if !doublegroup
        #        compute_asymstates_N( orbital , n , cg_o_dir , multiplets_path ; verbose=verbose )
        #    else 
        #        compute_asymstates_N_doublegroup( orbital , n , cg_o_dir , multiplets_path ; verbose=verbose )
        #    end
        #else 
        # if isdoublegroup(symmetry)
        #     compute_asymstates_N_doublegroup(
        #         irrep,
        #         n,
        #         clebschgordan_path,
        #         multiplets_path;
        #         verbose=verbose
        #     )
        # elseif ispointspin(symmetry)
        #     compute_asymstates_N_pointspin(
        #         irrep,
        #         n,
        #         clebschgordan_path,
        #         multiplets_path;
        #         verbose=verbose
        #     )
        # elseif istotalangularmomentum(symmetry)
        #     compute_asymstates_N_totalangularmomentum(
        #         irrep,
        #         n,
        #         multiplets_path;
        #         verbose=verbose
        #     )
        #
        # else
        #     error("Chosen symmetry $symmetry is none of the available.")
        # end
    end

end
