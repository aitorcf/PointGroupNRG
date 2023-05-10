include( "thermo.jl" )

# generate Nz values of Z
function generate_Z( Nz::Int64 )
    return Float64[ i/Nz for i in 0:(Nz-1) ]
end

# average thermodynamic calculations over values of z
function zavg_thermo( label::String , Z::Vector{Float64} )

    # data from all z
    thermo_data_clean_all_z = Matrix{Float64}[]
    thermo_data_imp_all_z   = Matrix{Float64}[]
    thermo_data_diff_all_z  = Matrix{Float64}[]

    # read and store data from all z
    for z in Z

        # clean
        thermo_clean_z_filename = thermo_filename_one_z( label , "clean" , z )
        thermo_clean_z = readdlm( thermo_clean_z_filename , skipstart=1 ) 
        push!( thermo_data_clean_all_z , thermo_clean_z )

        # imp
        thermo_imp_z_filename = thermo_filename_one_z( label , "imp" , z )
        thermo_imp_z   = readdlm( thermo_imp_z_filename , skipstart=1 ) 
        push!( thermo_data_imp_all_z , thermo_imp_z )

        # diff
        thermo_diff_z_filename = thermo_filename_one_z( label , "diff" , z )
        thermo_diff_z = readdlm( thermo_diff_z_filename , skipstart=1 )
        push!( thermo_data_diff_all_z , thermo_diff_z )

    end

    # interpolate all z data
    #
    # define new data variables
    thermo_data_clean_all_z_interpolated = Matrix{Float64}[]
    thermo_data_imp_all_z_interpolated = Matrix{Float64}[]
    thermo_data_diff_all_z_interpolated = Matrix{Float64}[]
    # collect temperatures from all z
    temperatures_zavg = sort(vcat(collect( data[:,1] for data in thermo_data_diff_all_z )...))
    # bound temperatures so that every z participates in average
    minmax_temperature = minimum([ maximum( data[:,1] ) for data in thermo_data_diff_all_z ])
    maxmin_temperature = maximum([ minimum( data[:,1] ) for data in thermo_data_diff_all_z ])
    filter!( x->(x<=minmax_temperature) , temperatures_zavg )
    filter!( x->(x>=maxmin_temperature) , temperatures_zavg )
    # interpolate for every z
    for i in 1:length(Z)

        # data and temperatures
        clean_old = thermo_data_clean_all_z[i]
        imp_old = thermo_data_imp_all_z[i]
        diff_old = thermo_data_diff_all_z[i]
        temperatures_old = diff_old[:,1]

        # interpolate 
        clean_interpolated = interpolate_thermo_matrix( clean_old , temperatures_zavg )
        imp_interpolated   = interpolate_thermo_matrix( imp_old   , temperatures_zavg )
        diff_interpolated  = interpolate_thermo_matrix( diff_old  , temperatures_zavg )

        # store
        push!( thermo_data_clean_all_z_interpolated , clean_interpolated )
        push!( thermo_data_imp_all_z_interpolated   , imp_interpolated   )
        push!( thermo_data_diff_all_z_interpolated  , diff_interpolated  )

    end

    # average over interpolated data
    #
    # create new matrices
    thermo_clean_zavg = zeros( Float64 , size(thermo_data_clean_all_z_interpolated[1])... )
    thermo_imp_zavg   = zeros( Float64 , size(thermo_data_clean_all_z_interpolated[1])... )
    thermo_diff_zavg  = zeros( Float64 , size(thermo_data_clean_all_z_interpolated[1])... )
    # store temperatures
    thermo_clean_zavg[:,1] = temperatures_zavg
    thermo_imp_zavg[:,1]   = temperatures_zavg
    thermo_diff_zavg[:,1]  = temperatures_zavg
    # average thermodynamic quantities
    thermo_clean_zavg[:,2:end] = length(Z)^-1 * sum( data_z[:,2:end] for data_z in thermo_data_clean_all_z_interpolated )
    thermo_imp_zavg[:,2:end]   = length(Z)^-1 * sum( data_z[:,2:end] for data_z in thermo_data_imp_all_z_interpolated )
    thermo_diff_zavg[:,2:end]  = length(Z)^-1 * sum( data_z[:,2:end] for data_z in thermo_data_diff_all_z_interpolated )

    # store z-averaged data
    #
    # z-averaged clean data
    zavg_clean_filename = thermo_filename_zavg( label , "clean" )
    write_thermo_data( zavg_clean_filename , thermo_clean_zavg )
    # z-averaged imp data
    zavg_imp_filename = thermo_filename_zavg( label , "imp" )
    write_thermo_data( zavg_imp_filename , thermo_imp_zavg )
    # z-averaged impurity contribution
    zavg_diff_filename = thermo_filename_zavg( label , "diff" )
    write_thermo_data( zavg_diff_filename , thermo_diff_zavg )

end


function zavg_spectral( label::String , 
                        Z::Vector{Float64} ; 
                        orbitalresolved::Bool=false ,
                        No::Int64=0 )

    println( "Averaging over values of z..." )

    Nz = length(Z)

    if orbitalresolved

        for i in 1:No

            data = Dict()
            for z in Z 
                filename = "spectral/spectral_$(label)_o$(i)_z$(z).dat"
                data[z] = readdlm( filename )
            end

            omegas = data[Z[1]][:,1]
            data_zavg = [0.0 for _=1:length(omegas)]
            for z in Z
                data_zavg .+= data[z][:,2]./Nz
            end

            zavgfile = "spectral/spectral_$(label)_o$(i)_zavg.dat"
            open( zavgfile , write=true ) do f
                writedlm( f , [omegas data_zavg] )
            end

        end

    else 

        data = Dict()
        for z in Z 
            data[z] = readdlm( "spectral/spectral_$(label)_z$(z).dat" )
        end

        omegas = data[Z[1]][:,1]

        data_zavg = [0.0 for _=1:length(omegas)]
        for z in Z
            data_zavg .+= data[z][:,2]./Nz
        end

        zavgfile = "spectral/spectral_$(label)_zavg.dat"
        open( zavgfile , write=true ) do f
            writedlm( f , [omegas data_zavg] )
        end

    end

end

using Interpolations


#function interpolate_spectral_function( 
#        xx::Vector{Float64} , 
#        yy::Vector{Float64} ;
#        step_reductor::Int64=100 )
#
#    # integer range for interpolation 
#    range_input = 1:length(yy)
#    range_output = 1:(1.0/(step_reductor-1)):length(yy)
#
#    # interpolate y
#    interpolator_y = cubic_spline_interpolation( range_input , yy )
#    interpolator_x = cubic_spline_interpolation( range_input , xx )
#    yy_dense = map( interpolator_y , range_output )
#    xx_dense = map( interpolator_x , range_output )
#
#    return xx_dense,yy_dense
#
#end
function zavg_spectral_new( label::String , 
                        Z::Vector{Float64} ; 
                        orbitalresolved::Bool=false ,
                        No::Int64=0 )

    println( "Averaging over values of z..." )

    Nz = length(Z)

    if orbitalresolved

        for i in 1:No

            data = Dict()
            for z in Z 
                filename = "spectral/spectral_$(label)_o$(i)_z$(z).dat"
                data[z] = readdlm( filename )
            end

            omegas = data[Z[1]][:,1]
            data_zavg = [0.0 for _=1:length(omegas)]
            for z in Z
                data_zavg .+= data[z][:,2]./Nz
            end

            zavgfile = "spectral/spectral_$(label)_o$(i)_zavg.dat"
            open( zavgfile , write=true ) do f
                writedlm( f , [omegas data_zavg] )
            end

        end

    else 

        data = Dict()
        omegas_lin_neg = collect(-1:0.001:-0.001)
        omegas_lin_pos = sort(-omegas_lin_neg)
        min_omega = maximum([minimum(abs.(readdlm( "spectral/spectral_$(label)_z$(z).dat" )[:,1])) for z in Z])
        omegas_log = collect( (sign*1.5^-i) for sign in [-1,1] for i in 1:100 )
        omegas_log = filter( x->abs(x)>min_omega , omegas_log )
        omegas = sort(vcat(omegas_lin_neg,omegas_lin_pos,omegas_log))
        @show omegas
        for z in Z 
            zmatrix_sparse = readdlm( "spectral/spectral_$(label)_z$(z).dat" )
            zmatrix_dense = hcat(interpolate_spectral_function( zmatrix_sparse[:,1] , zmatrix_sparse[:,2] )...)
            println( "***********************" )
            println()
            @show zmatrix_sparse[:,2]
            open( "spectral/test_z$(z).dat" , write=true ) do f
                writedlm( f , zmatrix_dense )
            end
            println()
            println( "***********************" )
            println()
            @show zmatrix_dense[:,2]
            println()
            println( "***********************" )
            println()
            data_interpolator = linear_interpolation( zmatrix_dense[:,1] , zmatrix_dense[:,2] , extrapolation_bc=Line() )
            data_interpolated = vcat(collect([o data_interpolator(o)] for o in omegas)...)
            @show data_interpolated
            println()
            println( "***********************" )
            data[z] = data_interpolated
            #data[z] = readdlm( "spectral/spectral_$(label)_z$(z).dat" )
        end

        #omegas = data[Z[1]][:,1]

        data_zavg = [0.0 for _=1:length(omegas)]
        for z in Z
            data_zavg .+= data[z][:,2]./Nz
        end
        
        #data_zavg = vcat(collect(values(data))...)
        #data_zavg = [ [data_zavg[i,1], data_zavg[i,2]] for i in 1:size(data_zavg,1) ]
        #@show data_zavg
        #println(); println()
        #sort!( data_zavg , by=x->x[1] )
        #data_zavg = hcat( data_zavg... )'

        zavgfile = "spectral/spectral_$(label)_zavg.dat"
        open( zavgfile , write=true ) do f
            #writedlm( f , [omegas data_zavg] )
            writedlm( f , [omegas data_zavg] )
        end

    end

end
