function generate_Z( Nz::Int64 )
    return Float64[ i/Nz for i in 0:(Nz-1) ]
end

function zavg_thermo( label::String , Z::Vector{Float64} )

    Nz = length( Z )

    println( "Averaging thermodynamic functions over values of z..." )
    th_tot  = Dict()
    thclean_tot = Dict()
    thimp_tot = Dict()

    for z in Z
        th_z   = readdlm( "thermodata/th_diff_$(label)_z$z.dat" )
        thclean_z = readdlm( "thermodata/thermo_clean_$(label)_z$z.dat" ) 
        thimp_z = readdlm( "thermodata/thermo_imp_$(label)_z$z.dat" ) 
        t = th_z[:,1] 
        th_tot[z] = Dict( round(t[i],sigdigits=2)=>th_z[i,:] for i in 1:length(t) )
        thclean_tot[z] = Dict( round(t[i],sigdigits=2)=>thclean_z[i,:] for i in 1:length(t) ) 
        thimp_tot[z] = Dict( round(t[i],sigdigits=2)=>thimp_z[i,:] for i in 1:length(t) ) 
    end

    th_zavg = Dict()
    thclean_zavg = Dict()
    thimp_zavg = Dict()
    T = sort(collect(keys(th_tot[Z[end]])))
    T = T[Nz:(end-Nz)]
    for t in T
        th_zavg[t]  = sum( th_tot[z][t] for z in Z )/Nz
        thclean_zavg[t] = sum( thclean_tot[z][t] for z in Z )/Nz 
        thimp_zavg[t] = sum( thimp_tot[z][t] for z in Z )/Nz 
    end

    println()
    
    th_zavg_vec = [th_zavg[t] for t in T]
    thclean_zavg_vec = [thclean_zavg[t] for t in T]
    thimp_zavg_vec = [thimp_zavg[t] for t in T]
    writedlm( "thermodata/th_zavg_$(label).dat" , th_zavg_vec )
    writedlm( "thermodata/th_cleanavg_$(label).dat" , thclean_zavg_vec )
    writedlm( "thermodata/th_impavg_$(label).dat" , thimp_zavg_vec )
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
