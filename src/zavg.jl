# generate Nz values of Z
function generate_Z( Nz::Int64 )
    #return Float64[ ( i/Nz - 0.5 + 1.0/(2*Nz) ) for i in 0:(Nz-1) ]
    return Float64[ (i/Nz)  for i in 0:(Nz-1) ]
end

# generate laps of Z containing number_of_cores z values
function generate_Zlaps( Z , number_of_cores )
    Nz = length(Z)
    cleandiv = Nz ÷ number_of_cores 
    rest = Nz % number_of_cores
    laps = rest==0 ? cleandiv : (cleandiv + 1)
    Zlaps = [ i=>[Z[j] for j in ((i-1)*number_of_cores+1):(i*number_of_cores) if j<=Nz] for i in 1:laps]
    return Zlaps
end

# average thermodynamic calculations over values of z
function zavg_thermo( label::String , 
                      Z::Vector{Float64} ;
                      average_from::String="" )

    source_label = average_from=="evenodd" ? label*"_evenodd" : label

    # data from all z
    thermo_data_clean_all_z = Matrix{Float64}[]
    thermo_data_imp_all_z   = Matrix{Float64}[]
    thermo_data_diff_all_z  = Matrix{Float64}[]

    # read and store data from all z
    for z in Z

        # clean
        thermo_clean_z_filename = thermo_filename_one_z( source_label , "clean" , z )
        thermo_clean_z = readdlm( thermo_clean_z_filename , skipstart=1 ) 
        push!( thermo_data_clean_all_z , thermo_clean_z )

        # imp
        thermo_imp_z_filename = thermo_filename_one_z( source_label , "imp" , z )
        thermo_imp_z   = readdlm( thermo_imp_z_filename , skipstart=1 ) 
        push!( thermo_data_imp_all_z , thermo_imp_z )

        # diff
        thermo_diff_z_filename = thermo_filename_one_z( source_label , "diff" , z )
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

function zavg_spectral( 
            label::String ,
            Z::Vector{Float64} ;
            orbitalresolved_number::Int64=0 )

    orbitalresolved = orbitalresolved_number>0

    # standard
    if !orbitalresolved

        data_all_z = Vector{Matrix{Float64}}()
        
        # gather spectral data
        for z in Z

            # read data
            filename = spectral_filename( label , z=z )
            data = readdlm( filename , skipstart=1 )

            # store data
            push!( data_all_z , data )

        end

        # linearly interpolate data for all z
        #
        # storage variable
        data_all_z_interpolated = Matrix{Float64}[]
        # collect energy values
        omegas_average = sort(vcat(collect( data[:,1] for data in data_all_z )...))
        filter!( x->(abs(x)<=1.0) , omegas_average )
        min_omega = maximum([ minimum(abs.(data[:,1])) for data in data_all_z ])
        filter!( x->(abs(x)>=min_omega) , omegas_average )
        sort!( omegas_average )
        # interpolate
        for data in data_all_z 

            # interpolating function
            interpolator = linear_interpolation( data[:,1] , data[:,2] , extrapolation_bc=Line() )

            # interpolation
            data_interpolated = [ omegas_average Float64[ interpolator(omega) for omega in omegas_average ]]

            # store interpolated data
            push!( data_all_z_interpolated , data_interpolated )

        end

        # average results 
        spectral_zavg = copy(data_all_z_interpolated[1])
        spectral_zavg[:,2] = sum( data[:,2] for data in data_all_z_interpolated )/length(Z)
            
        # remove duplicates
        spectral_zavg_noduplicates = [spectral_zavg[1,:]]
        for i in 1:size(spectral_zavg,1)
            row = spectral_zavg[i,:]
            duplicate = false
            for j in 1:length(spectral_zavg_noduplicates)
                row_noduplicates = spectral_zavg_noduplicates[j]
                if isapprox(row_noduplicates[1],row[1])
                    duplicate = true
                end
            end
            if !duplicate
                push!( spectral_zavg_noduplicates , row )
            end
        end
        spectral_zavg = Matrix(hcat(spectral_zavg_noduplicates...)')
            
        # smooth results 
        omegas_splined,spectral_zavg_splined = spline_interpolate_spectral_function( 
            spectral_zavg[:,1] ,
            spectral_zavg[:,2] ,
            orbitalresolved=orbitalresolved
       )

        # write results
        write_spectral_function(
            spectral_filename(label,zavg=true),
            spectral_zavg 
        )
        write_spectral_function(
            spectral_filename(label,zavg=true,tail="_splined") ,
            [omegas_splined spectral_zavg_splined]
        )

    elseif orbitalresolved 

        # same procedure for every orbital
        for orbital in 1:orbitalresolved_number

            # all z data for this orbital
            data_all_z = Vector{Matrix{Float64}}()

            # read header for orbital
            orbital_headervec = [""]
            open( spectral_filename(label,z=0.0,orbital=orbital) ) do f
                orbital_headervec[1] = readline(f)*"\n"
            end
            orbital_header = orbital_headervec[1]
            
            # gather spectral data
            for z in Z

                # read data
                filename = spectral_filename( label , z=z , orbital=orbital )
                data = read_spectral_data(filename)

                # store data
                push!( data_all_z , data )

            end

            # linearly interpolate data for all z
            #
            # storage variable
            data_all_z_interpolated = Matrix{Float64}[]
            # collect energy values
            omegas_average = sort(vcat(collect( data[:,1] for data in data_all_z )...))
            filter!( x->(abs(x)<=1.0) , omegas_average )
            min_omega = maximum([ minimum(abs.(data[:,1])) for data in data_all_z ])
            filter!( x->(abs(x)>=min_omega) , omegas_average )
            # interpolate
            for data in data_all_z 

                # interpolating function
                interpolator = linear_interpolation( data[:,1] , data[:,2] , extrapolation_bc=Line() )

                # interpolation
                data_interpolated = [ omegas_average Float64[ interpolator(omega) for omega in omegas_average ]]

                # store interpolated data
                push!( data_all_z_interpolated , data_interpolated )

            end


            # average results 
            spectral_zavg = copy(data_all_z_interpolated[1])
            spectral_zavg[:,2] = sum( data[:,2] for data in data_all_z_interpolated )/length(Z)

            # remove duplicates
            spectral_zavg_noduplicates = [spectral_zavg[1,:]]
            for i in 1:size(spectral_zavg,1)
                row = spectral_zavg[i,:]
                duplicate = false
                for j in 1:length(spectral_zavg_noduplicates)
                    row_noduplicates = spectral_zavg_noduplicates[j]
                    if isapprox(row_noduplicates[1],row[1])
                        duplicate = true
                    end
                end
                if !duplicate
                    push!( spectral_zavg_noduplicates , row )
                end
            end
            spectral_zavg = Matrix(hcat(spectral_zavg_noduplicates...)')

            # smooth results 
            #omegas_splined,spectral_zavg_splined = spline_interpolate_spectral_function( 
            #    spectral_zavg[:,1] ,
            #    spectral_zavg[:,2] ,
            #    orbitalresolved=orbitalresolved
            #)

            # write results
            write_spectral_function(
                spectral_filename(label,zavg=true,orbital=orbital),
                spectral_zavg,
                orbitalresolved_header=orbital_header
            )
            #write_spectral_function(
            #    spectral_filename(label,zavg=true,orbital=orbital,tail="_splined") ,
            #    [omegas_splined spectral_zavg_splined]
            #)
        end


    end

end
