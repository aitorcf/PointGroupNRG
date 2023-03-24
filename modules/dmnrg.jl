using LinearAlgebra 

function NRG( iterations::Int64, 
              cutoff_type::String, 
              cutoff_magnitude::Number,
              L::Float64,
              hop_symparams::Dict{ Int64 , Matrix{ComplexF64} },
              irrEU::Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}} },
              multiplets_shell::Set{NTuple{4,Int64}}, 
              cg_o_fullmatint::Dict{Tuple{Int64, Int64, Int64}, Array{ComplexF64, 3}},
              cg_s_fullmatint::Dict{Tuple{Int64, Int64, Int64}, Array{ComplexF64, 3}},
              Csum_o_array::Array{ComplexF64,6} ,
              Csum_s_array::Array{ComplexF64,6} ,
              Bsum_o_array::Array{ComplexF64,6} ,
              Bsum_s_array::Array{ComplexF64,6} ,
              pcgred_shell::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} },
              multiplets_a::Vector{NTuple{4,Int64}} , 
              combinations_uprima::Dict{NTuple{3,Int64}, Vector{NTuple{3,NTuple{4,Int64}}}},
              betabar::Float64 ,
              oindex2dimensions::Vector{Int64} ,
              xi_symparams::Dict{ Int64 , Vector{Vector{ComplexF64}} } ,
              mm_i ;
              verbose::Bool=false ,
              distributed::Bool=false ,
              method::String="distfor" ,
              minmult::Int64=0 , 
              mine::Float64=0.0 ,
              z::Float64=0.0 ,
              discretization::String="standard" ,
              spectral::Bool=false ,
              M::Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} }=Dict{ NTuple{3,NTuple{3,Int64}} , Array{ComplexF64,3} }() ,
              AA::Vector{Dict{ NTuple{4,Int64} , Tuple{Float64,Vector{Float64}} }}=Dict{ NTuple{4,Int64} , Tuple{Float64,Vector{Float64}} }[] ,
              Karray_orbital::Array{ComplexF64,6}=Array{ComplexF64,6}(undef,0,0,0,0,0,0) , 
              Karray_spin::Array{ComplexF64,6}=Array{ComplexF64,6}(undef,0,0,0,0,0,0) ,
              multiplets_atomhop::Vector{NTuple{4,Int64}}=NTuple{4,Int64}[] ,
              etafac::Float64=1.0 ,
              alpha::Float64=1.0 ,
              eta::Function=x->1.0 ,
              Nz=1 ,
              dmnrg=false )

    println( "=============" )
    println( "NRG PROCEDURE" )
    println( "=============" )
    println()

    temperatures   = []
    magnetizations = []
    energies       = []
    numbers        = []
    partitions     = []
    entropies      = []
    heatcaps       = []
    free_energies  = []
    impspins       = []
    impnums        = []
    impmults       = []

    performance = []
    maxes = []
    maxe_tot = 0.0
    maxn_tot = 0    
    maxs_tot = 0
    eigenvals5::Matrix{Float64} = zeros(Float64,iterations-1,5)

    if discretization=="lanczos" 
        e,ebar,ξ = get_hoppings( iterations , L , z , eta )
    else
        ξ = compute_xi_vector( iterations , z , L ; discretization="standard" )
    end

    # DMNRG 
    irrEUs = Dict{ NTuple{3,Int64} , Tuple{Vector{Float64},Matrix{ComplexF64}}}[]

    for n=2:iterations

        println( "ITERATION n=$n" )

        # cut off and block block multiplets
        print( "Applying cutoff to obtain block multiplets... " )
        @time (multiplets_block, discarded) = cut_off!( irrEU ; 
                                                        type=cutoff_type , 
                                                        cutoff=cutoff_magnitude , 
                                                        safeguard=true ,
                                                        minmult=minmult ,
                                                        mine=mine ,
                                                        verbose=false ,
                                                        M=M )

        # post-cutoff diagonalization info for DMNRG.
        # the first element corresponds to N=1 !!!
        push!( irrEUs , irrEU )

        # cutoff info
        equivstates = Int64( sum( (m[3]+1) for m in multiplets_block ) )
        println( "$(length(multiplets_block)) multiplets kept ($equivstates states), $(length(discarded)) multiplets discarded" )
        proportion = length(multiplets_block)/Float64( length(multiplets_block) + length(discarded) ) * 100
        # maximum energy after cutoff
        # energy
        maxe = maximum(collect( e for (G,(E,U)) in irrEU for e in E ))
        maxe_tot = maximum([ maxe , maxe_tot ])
        push!( maxes , maxe )

        println( "proportion: $(proportion)%. max energy: $maxe." )

        # renormalize by √Λ
        println( "Renormalizing eigenvalues...")
        for (G,(E,U)) in irrEU 
            irrEU[G] = ( E.*sqrt(L) , U )
        end

        # hopping parameter
        if discretization=="lanczos"
            hop_symparams = Dict( k=>diagm(ComplexF64.([xi_symparams[k][i][n-1]
                                                        for i in 1:length(xi_symparams[k])]))
                                  for (k,v) in xi_symparams )
        else
            hop_symparams = Dict( k=>ComplexF64.(ξ[n-1]*Matrix(LinearAlgebra.I,size(v)...)) # at n=2, we want ξ[1]=ξ_0
                                  for (k,v) in hop_symparams )
        end

        # construct and diagonalize ( m_u | H_1 | m_v )
        println( "Diagonalizing Hamiltonian..." )
        ppp = @timed (irrEU,combinations_uprima) = matdiag_redmat( 
                multiplets_block , 
                multiplets_shell ,
                irrEU , 
                hop_symparams , 
                cg_o_fullmatint , 
                cg_s_fullmatint ,
                Csum_o_array ,
                Csum_s_array ,
                Bsum_o_array ,
                Bsum_s_array ,
                pcgred_shell ,
                pcgred_shell ,
                multiplets_a , 
                multiplets_a ,
                combinations_uprima , 
                oindex2dimensions ;
                verbose=verbose ,
                distributed=distributed );
        @show ppp.time, ppp.bytes*10^-6, ppp.gctime
        push!( performance , ppp )
        
        # maximum values for symmetry quantum numbers
        # spin
        maxs = maximum(collect( G[3] for (G,(E,U)) in irrEU ))
        maxs_tot = maximum([ maxs , maxs_tot ])
        # occupation 
        maxn = maximum(collect( G[1] for (G,(E,U)) in irrEU ))
        maxn_tot = maximum([ maxn , maxn_tot ])

        # impurity info
        mm_i = imp_mults( irrEU ,
                          oindex2dimensions ,
                          combinations_uprima ,
                          mm_i )
        m_imp = mult_thermo( irrEU ,
                             betabar ,
                             oindex2dimensions ,
                             mm_i )
        @show m_imp
        push!( impmults , m_imp )

        # thermodynamics 
        println( "THERMODYNAMICS" )
        fac = L
        if discretization!=="lanczos" 
            t = temperature( n , L , betabar*fac ; z=z , discretization=discretization )
        elseif discretization=="lanczos" 
            t = temperature( n , L , betabar*fac ; z=z , discretization=discretization , ebar0=ebar[1] )
        end
        ρ    = partition(  irrEU , betabar*fac , oindex2dimensions )
        entr = entropy(    irrEU , betabar*fac , oindex2dimensions )
        mag  = magsusc(    irrEU , betabar*fac , oindex2dimensions )
        N    = number(     irrEU , betabar*fac , oindex2dimensions )
        en   = energy(     irrEU , betabar*fac , oindex2dimensions )
        c    = heatcap(    irrEU , betabar*fac , oindex2dimensions )
        f    = free_energy(irrEU , betabar*fac , oindex2dimensions )
        @show t 
        @show ρ 
        @show entr
        @show mag
        @show N
        @show en
        @show c 
        @show f
        println()
        push!( magnetizations , mag )
        push!( temperatures , t )
        push!( energies , en )
        push!( partitions , ρ )
        push!( numbers , N )
        push!( entropies , entr )
        push!( heatcaps , c )
        push!( free_energies , f )

        eigs = sort([ e for (E,U) in values(irrEU) for e in E ])[1:5]
        eigenvals5[(n-1),:] = [(e-eigs[1]) for e in eigs[1:5]]

        if spectral 

                print( "Updating M and AA... " )
                @time M, AA = update_redmat_AA_CGsummethod(
                            M ,
                            irrEU ,
                            combinations_uprima ,
                            collect(multiplets_atomhop) ,
                            cg_o_fullmatint ,
                            cg_s_fullmatint ,
                            Karray_orbital ,
                            Karray_spin ,
                            AA ,
                            oindex2dimensions ;
                            verbose=false )
                #M,AA = update_redmat_AA(
                #            M,
                #            irrEU ,
                #            combinations_uprima ,
                #            multiplets_atomhop ,
                #            cg_o_fullmatint ,
                #            cg_s_fullmatint ,
                #            AA ,
                #            oindex2dimensions ;
                #            verbose=false )
                println()
                 
        end

    end

    # add last irrEU
    push!( irrEUs , irrEU )

    maxe_avg = sum(maxes)/length(maxes)
    @show maxe_avg
    println( "Maximum irrep qnums: N=$maxn_tot , 2S=$maxs_tot\n" )

    spec = zeros(Float64,1,1)
    if spectral 

        number_of_hoppers = sum( oindex2dimensions[m_a[2]]*(m_a[3]+1)
                                 for m_a in multiplets_atomhop )
        spec = compute_spectral_function(
            AA ,
            L ,
            iterations ,
            alpha ,
            etafac ;
            widthfac=maximum((L,maxe_avg)) ,
            number_of_hoppers=number_of_hoppers )
    end

    return ( t=temperatures , 
             m=magnetizations , 
             e=energies , 
             p=partitions ,
             n=numbers , 
             c=heatcaps ,
             f=free_energies ,
             entr=entropies,
             perf=performance ,
             impmults=impmults ,
             specfunc=spec ,
             irrEUs=irrEUs )
end

function compute_reddens( 
            irrEUs ,
            multiplets_shell ) 

    firstrho = Dict( G=>Matrix{Float64}(I,size(U)...)
                     for (G,(E,U)) in irrEUs[end] )
    N = length(irrEUs)

    for n=(N-1):-1:1,
        (G,(E,U)) in irrEUs[n] 

        # kept and total sizes
        K = length(E)
        T = size(U,1)

        # kept part 

        
        # discarded part 


    end

end

