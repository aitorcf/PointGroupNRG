using PointGroupNRG
using PointGroupNRG.NRGCalculator
using DelimitedFiles
using Interpolations
using QuadGK
using Distributed

import Base.length
import Base.real
import Base.imag
import Base.:-
import Base.iszero

const integral_steps = 100000
const delta = 1e-4
const integral_tolerance = 1e-10
const convergence_tolerance = 1e-1


# ------------- #
# Hybridization #
# ------------- #

struct Hybridization
    y::Vector{Float64}
    x::Vector{Float64}
end
function Hybridization( hybridization::Function ,
                        energies::Vector{Float64} )
    x = copy(energies)
    y = map( hybridization , energies ) 
    return Hybridization(y,x)
end
function Hybridization( hybridization::Function ;
                        x::Vector{Float64}=Float64[],
                        nsteps::Int64=integral_steps )
    if length(x)==0
        return Hybridization( hybridization , collect(-1.0:(2.0/nsteps):1.0) )
    else
        return Hybridization( map(hybridization,x) , x )
    end
end
function Hybridization( hybridization::Vector{Float64} )
    nsteps = length(hybridization)-1
    energies = collect(-1.0:(2.0/nsteps):1.0)
    return Hybridization( hybridization , energies )
end

# array utilities
length( h::Hybridization ) = length(h.x)
getindex( h::Hybridization , idx::Int64 ) = (h.y[idx],h.x[idx])
function setindex!( h::Hybridization , v::Float64 , idx::Int64 )
    h.y[idx] = v
end
function setindex!( h::Hybridization , v::NTuple{2,Float64} , idx::Int64 )
    h.y[idx],h.x[idx] = v
end

function filesave( filename::String ,
                   h::Hybridization )
    writedlm( filename , hcat(h.x,h.y) )
end

function functionalize( h::Hybridization )::Function
    interpolator = linear_interpolation( h.x , h.y ; extrapolation_bc=Line() )
    return x->interpolator(x)
end


# ----------- #
# Self-energy #
# ----------- #

struct SelfEnergy
    y_real::Vector{Float64}
    y_imag::Vector{Float64}
    x::Vector{Float64}
end
# Σ(ω) = 0
function SelfEnergy( x::Vector{Float64} )
    y_real = zeros(ComplexF64,size(x))
    y_imag = zeros(ComplexF64,size(x))
    return SelfEnergy(y_real,y_imag,x)
end
# Σ_h(ω) = 1/π⋅∫dϵΔ(ϵ)/(ω-ϵ+iδ)
# coleman p. 593
function SelfEnergy( h::Hybridization ;
                     integral_tolerance::Float64=integral_tolerance )

    # imaginary part
    selfenergy_imag = - copy(h.y)

    # real part
    hybridization_function = linear_interpolation(h.x,h.y;extrapolation_bc=Line())
    integrand(o,e) = hybridization_function(e)/(o-e+im*delta)
    #convolution = [ integrate_dmft(e->integrand(o,e),(-1,1)) for o in h.x ] 
    convolution = [ quadgk(e->integrand(o,e),-1,1,rtol=integral_tolerance)[1] for o in h.x ] 
    selfenergy_real = real.(convolution)./pi

    ## real part with integral
    #d = [(h.x[i+1]-h.x[i]) for i in 1:(length(h.x)-1)]
    #push!( d , 0.0 )
    #kernel = reduce(hcat,[
    #    [1.0/(omega-epsilon+im*d[i]) for omega in h.x]
    #    for (i,epsilon) in enumerate(h.x)
    #]) 
    #selfenergy_real = real.(convolute(h.y,kernel,h.x))./pi

    return SelfEnergy( selfenergy_real , selfenergy_imag , copy(h.x) )
end

function filesave( filename::String ,
                   se::SelfEnergy )
    writedlm( filename , hcat(se.x,se.y_real,se.y_imag) )
end

-(se::SelfEnergy) = SelfEnergy( -se.y_real , -se.y_imag , se.x )
iszero(se::SelfEnergy) = iszero(sum(abs2.(se.y_real+se.y_imag)))
imag(se::SelfEnergy) = copy(se.y_imag)

# extract hybridization from hybridization self-energy
function Hybridization( se::SelfEnergy )
    Hybridization( -imag(se) , copy(se.x) )
end

# --------------- #
# Green functions #
# --------------- #

struct Green
    y_denominator_real::Vector{Float64} 
    y_denominator_imag::Vector{Float64}
    x::Vector{Float64}
end
function Green( mu::Float64 ,
                x::Vector{Float64} ;
                delta::Float64=delta )
    y_denominator_real = x .- mu
    y_denominator_imag = map( a->delta , x )
    return Green( y_denominator_real , y_denominator_imag , copy(x) )
end
# G(ω) = 1/(ω-μ-Σ(ω)+iδ)
function Green( mu::Float64 ,
                se::SelfEnergy ;
                delta::Float64=delta )
    y_denominator_real = @. se.x-mu-se.y_real
    y_denominator_imag = -se.y_imag .+ delta
    x = copy(se.x)
    return Green( y_denominator_real , y_denominator_imag , x )
end
# G = 1/(ω-μ-1/π∫dϵΔ(ϵ)/(ω-μ+iδ))
function Green( mu::Float64 , 
                h::Hybridization )
    return Green(mu,SelfEnergy(h))
end
function Green( mu::Float64 ,
                h::Function ,
                x::Vector{Float64} ;
                integral_tolerance::Float64=integral_tolerance )

    # hybridization self-energy convolution
    denominator(o,e) = (o-mu-e)
    integrand(o,e) = h(e)/denominator(o,e)
    convolution = zeros(ComplexF64,size(x))
    for (i,o) in enumerate(x)

        # normal integration
        margin = 0.01
        blowup = o-mu
        h_blowup = h(blowup)
        subtracted = h_blowup*log(abs(denominator(o,x[1]-margin)/denominator(o,x[end]+margin)))
        preblow  = blowup<=x[1] ? 0.0 : quadgk( 
            e->(h(e)-h(blowup))/denominator(o,e) +
            subtracted/(blowup-x[1]) ,
            x[1] , 
            blowup )[1] 
        postblow = blowup>=x[end] ? 0.0 : quadgk(
            e->(h(e)-h_blowup)/denominator(o,e) +
            subtracted/(x[end]-blowup) ,
            blowup,
            x[end] )[1]
        convolution[i] = preblow + postblow - im*pi*h_blowup

    end

    # denominators
    y_denominator_real = x .- mu .- real.(convolution)./pi
    y_denominator_imag = map( h , x )

    return Green( y_denominator_real , y_denominator_imag , x )
end
# G(ω)^{-1} = W(ω)^{-1} - Σ(ω)
function Green( g::Green , 
                se::SelfEnergy ; 
                delta::Float64=delta )
    denominator_real = g.y_denominator_real .- se.y_real
    denominator_imag = g.y_denominator_imag .- se.y_imag .+ delta
    return Green( denominator_real , denominator_imag , copy(g.x) )
end

# Σ = W^{-1} - G^{-1}
function SelfEnergy( w::Green , g::Green )

    # energy values
    x = copy(w.x)

    # interpolate if needed
    interpolate = false
    if (length(w.x)!==length(g.x) || !isapprox(w.x,g.x))
        interpolate = true 
        x = sort(collect(Set(vcat(g.x,w.x)))) # combined energy values 
        weiss_interpolator_real = linear_interpolation(w.x,w.y_denominator_real;extrapolation_bc=Line())
        weiss_interpolator_imag = linear_interpolation(w.x,w.y_denominator_imag;extrapolation_bc=Line())
        weiss_real = map( weiss_interpolator_real , x )
        weiss_imag = map( weiss_interpolator_imag , x )
        green_interpolator_real = linear_interpolation(g.x,g.y_denominator_real;extrapolation_bc=Line())
        green_interpolator_imag = linear_interpolation(g.x,g.y_denominator_imag;extrapolation_bc=Line())
        green_real = map( green_interpolator_real , x )
        green_imag = map( green_interpolator_imag , x )
        w = Green( weiss_real , weiss_imag , x )
        g = Green( green_real , green_imag , x )
    end

    y_real = w.y_denominator_real - g.y_denominator_real
    y_imag = w.y_denominator_imag - g.y_denominator_imag

    return SelfEnergy(y_real,y_imag,x)
end



imag(w::Green) = @. -w.y_denominator_imag/
                    (w.y_denominator_real^2+w.y_denominator_imag^2)
real(w::Green) = @. w.y_denominator_real/
                   (w.y_denominator_real^2+w.y_denominator_imag^2)

function diff_absolute( g1::Green ,
                        g2::Green )

    if ( length(g1.x)!==length(g2.x) || !isapprox(sum(abs2.(g1.x-g2.x)),zero(Float64)) )
        g1_real_interpolator = linear_interpolation( g1.x , g1.y_denominator_real ;extrapolation_bc=Line())
        g1_imag_interpolator = linear_interpolation( g1.x , g1.y_denominator_imag ;extrapolation_bc=Line())
        g1_real = map( g1_real_interpolator , g2.x )
        g1_imag = map( g1_imag_interpolator , g2.x )
        g1 = Green( g1_real , g1_imag , g2.x )
    end

    diff_real = sum(abs2.(g1.y_denominator_real-g2.y_denominator_real))
    diff_imag = sum(abs2.(g1.y_denominator_imag-g2.y_denominator_imag))

    return sqrt(diff_real + diff_imag)
end
function diff_relative( g1::Green ,
                        g2::Green )

    g1_interpolator = linear_interpolation( g1.x , imag(g1) )
    g2_interpolator = linear_interpolation( g2.x , imag(g2) )
    diff_interpolated(x) = abs(g1_interpolator(x)-g2_interpolator(x))

    integral  = quadgk( diff_interpolated      , g1.x[1],g1.x[end] )[1]
    reference = quadgk( x->-g1_interpolator(x) , g1.x[1],g1.x[end] )[1]

    return integral/reference
end

# Δ(ω) = Im G1^{-1}
function Hybridization( w::Green )
    return Hybridization( w.y_denominator_imag , w.x )
end

function filesave( filename::String ,
                   g::Green )
    writedlm( filename , hcat(g.x,real(g),imag(g),g.y_denominator_real,g.y_denominator_imag) )
end


# ----------------------------- #
# Local DOS / Spectral function #
# ----------------------------- #

struct DOS
    y::Vector{Float64}
    x::Vector{Float64}
end
# ρ(ϵ) = -1/π Im G^{-1}
DOS( g::Green ) = DOS( -imag(g)./pi , g.x )
# from file
function DOS( filename::String )
    data = readdlm( filename , comments=true , comment_char='#' )
    data_x = data[:,1]
    data_y = data[:,2] # 2 = spin degeneracy of one-electron excitations 
    return DOS(data_y,data_x)
end
# from function
function DOS( f::Function ;
              x::Vector{Float64}=Float64[] ,
              nsteps::Int64=integral_steps )
    if length(x)==0
        xx = collect(-1.0:(2.0/nsteps):1.0)
    else 
        xx = copy(x)
    end
    y::Vector{Float64} = map( f , x )
    return DOS( y , xx )
end

integrate_dmft( dos::DOS ) = integrate_dmft(dos.x,dos.y)
function normalize( dos::DOS ; impurity_occupation::Float64=1.0 ) 
    xx = filter( x->x<0 , dos.x )
    yy = dos.y[eachindex(xx)]
    return DOS( dos.y ./ integrate_dmft(xx,yy) .* impurity_occupation , dos.x )
end

function filesave( filename::String ,
                   dos::DOS )
    writedlm( filename , hcat(dos.x,dos.y) )
end

# G(ω) = ∫dϵ ρ(ϵ)/(ω-μ-ϵ-Σ(ω)+iδ)
function Green( mu::Float64 ,
                dos::DOS ;
                se::SelfEnergy=SelfEnergy(dos.x) ,
                delta::Float64=delta ,
                integral_tolerance::Float64=integral_tolerance )

    # set energy values to self-energy basis (NRG)
    x = copy(se.x)

    # functions
    dos_function     = linear_interpolation(dos.x,dos.y;extrapolation_bc=Line())
    se_real_function = linear_interpolation(se.x,se.y_real;extrapolation_bc=Line())
    se_imag_function = linear_interpolation(se.x,se.y_imag;extrapolation_bc=Line())

    # convolution intergral
    #integrand(o,e) = dos_function(e)/(o-mu-e-se_real_function(o)-im*se_imag_function(o)+im*delta)
    #convolution = [ quadgk(e->integrand(o,e),x[1],x[end],rtol=integral_tolerance)[1] for o in x ] 
    denominator(o,e) = (o-mu-e-se_real_function(o)-im*(se_imag_function(o)-delta))
    integrand(o,e) = dos_function(e)/denominator(o,e)
    convolution = zeros(ComplexF64,size(x))
    for (i,o) in enumerate(x)

        # normal integration
        if !iszero(se_imag_function(o))

            convolution[i] = quadgk( e->integrand(o,e) , x[1] , x[end] )[1]
            #peak = o-mu-se.y_real[i]
            #subtracted = dos_function(peak)*
            #             ( 0.5*log((o-mu-x[end]-se_real_function(o))^2+se_imag_function(o)^2) -
            #               im*atan((o-mu-x[end]-se_real_function(o))/se_imag_function(o)) -
            #               0.5*log((o-mu-x[1]  -se_real_function(o))^2+se_imag_function(o)^2) +
            #               im*atan((o-mu-x[1]  -se_real_function(o))/se_imag_function(o)) )
            #prepeak_subtracted  = subtracted/(peak-x[1])
            #postpeak_subtracted = subtracted/(x[end]-peak)
            #prepeak = peak==x[1] ? 0.0 : quadgk( 
            #    e->(dos_function(e)-dos_function(peak))/denominator(o,e) +
            #    prepeak_subtracted ,
            #    x[1] , 
            #    peak 
            #)[1]
            #postpeak = peak==x[end] ? 0.0 : quadgk( 
            #    e->(dos_function(e)-dos_function(peak))/denominator(o,e) +
            #    postpeak_subtracted ,
            #    peak , 
            #    x[end] 
            #)[1]
            #convolution[i] = prepeak + postpeak

        # improper integration
        else

            margin::Float64 = 1e-6

            blowup = o-mu-se.y_real[i]
            subtracted = dos_function(blowup)*log(abs(denominator(o,x[1]-margin)/denominator(o,x[end]+margin)))
            preblow  = blowup==x[1] ? 0.0 : quadgk( 
                e->(dos_function(e)-dos_function(blowup))/denominator(o,e) +
                    subtracted/(blowup-x[1]) ,
                x[1] , 
                blowup )[1] 
            postblow = blowup==x[end] ? 0.0 : quadgk(
                e->(dos_function(e)-dos_function(blowup))/denominator(o,e) +
                    subtracted/(x[end]-blowup) ,
                blowup,
                x[end])[1]
            convolution[i] = preblow + postblow - im*pi*dos_function(blowup)

        end
    end

    # denominators
    y_denominator_real = real.(convolution.^-1)
    y_denominator_imag = imag.(convolution.^-1)

    return Green( y_denominator_real , y_denominator_imag , x )
end

function functionalize( dos::DOS )::Function
    interpolator = linear_interpolation( dos.x , dos.y )
    return x->interpolator(x)
end

function bound_energies( dos::DOS ;
                         boundary::Float64=1.0 )
    x_and_y = collect(zip(dos.x,dos.y))
    filter!( x->abs(x[1])<=boundary , x_and_y )
    x = first.(x_and_y)
    y = last.(x_and_y)

    y_interpolator = linear_interpolation( dos.x , dos.y )
    positive_boundary_value = y_interpolator(boundary)
    negative_boundary_value = y_interpolator(-boundary)
    push!( x , boundary )
    push!( y , positive_boundary_value )
    insert!( x , 1 , -boundary )
    insert!( y , 1 , negative_boundary_value)

    return DOS(y,x)
end

# self-energy trick
function SelfEnergy( T::Green ,
                     G::Green ,
                     u::Float64 )
    x = copy(T.x)

    t = real(T) .+ im*imag(T)
    g = real(G) .+ im*imag(G)

    se = u .* t ./ g
    y_real = real(se)
    y_imag = imag(se)

    return SelfEnergy( y_real , y_imag , x )
end

# ----------- #
# Integration #
# ----------- #

function integrate_dmft(data::Matrix{<:Number}) 
    sum( (data[i+1,1]-data[i,1])*data[i,2] for i in 1:(size(data,1)-1) )
end
function integrate_dmft( xx::Vector{<:Number} , yy::Vector{<:Number} )
    return integrate_dmft(hcat(xx,yy))
end
function integrate_dmft( f::Function , 
                         range::NTuple{2,Float64} ; 
                         nsteps::Int64=integral_steps )

    # integral boundaries
    bound_lower = range[1]
    bound_upper = range[2]

    # interval
    step = (bound_upper-bound_lower)/nsteps

    # main loop
    x = bound_lower
    integral::ComplexF64 = 0.0
    for _ in 1:nsteps
        integral += step*f(x)
        x += step
    end

    return integral
end

# C(ω) = ∫dϵ K(ω,ϵ) f(ϵ)
function convolute( f::Vector{<:Number} , 
                    k::Matrix{<:Number} ,
                    epsilons::Vector{<:Number} ) 

    # C(ω)
    result = zeros(ComplexF64,size(f))

    # integral loop
    for i in eachindex(result)
        result[i] = integrate_dmft(epsilons,k[i,:].*f[:])
    end


    return result
end

# Σ_h(ω) = 1/π⋅∫dϵΔ(ϵ)/(ω-ϵ+iδ)
function hybridization2selfenergy( hybridization::Vector{<:Number} ; 
                                   nsteps::Int64=integral_steps , 
                                   delta::Float64=delta )

    # hybridization array
    omegas = collect(-1.0:(2.0/nsteps):1.0)

    # self-energy
    selfenergy = zeros(ComplexF64,size(hybridization))

    # integrand
    integrand = zeros(ComplexF64,size(hybridization))

    # self energy loop
    for (io,omega) in enumerate(omegas)
        @. integrand = hybridization/(omega-omegas+im*delta)
        selfenergy[io] = integrate_dmft(omegas,integrand)/pi
    end

    return selfenergy
end
function hybridization2selfenergy( hybridization::Function ; 
                                   nsteps::Int64=integral_steps , 
                                   delta::Float64=delta )

    # hybridization array
    omegas = collect(-1.0:(2.0/nsteps):1.0)
    hybridization_vector = map( hybridization , omegas )

    return hybridization2selfenergy(hybridization_vector)
end

# W(ω) = 1/(ω-μ-Σ_h(ω))
function hybridization2weiss( hybridization::Vector{<:Number} ,
                              mu::Float64 )

    # integral intervals
    nsteps = length(selfenergy)-1

    # self-energy for the denominator
    selfenergy = hybridization2selfenergy( hybridization )

    # weiss field and energies
    weiss = zeros(ComplexF64,size(selfenergy))
    omegas = collect(-1.0:(2.0/nsteps):1.0)
    @. weiss = 1/( omegas - mu - selfenergy + delta*im )

    return weiss 
end
function hybridization2weiss( hybridization::Function ,
                              mu::Float64 ;
                              nsteps::Float64=integral_steps )
    omegas = collect(-1.0:(2.0/nsteps):1.0)
    hybridization_vector = map( hybridization , omegas )
    return hybridization2weiss( hybridization_vector , mu )
end

# Δ(ω) = Im 1/W(ω)
function weiss2hybridization( weiss::Vector{<:Number} )
    return imag.(weiss.^-1)
end


# --------------- #
# DMFT Calculator #
# --------------- #

function dmft(
            label::String ,
            L::Float64 ,
            iterations::Int64 , 
            cutoff_type::String ,
            cutoff_magnitude::R ,
            cg_o_dir::String ,
            multiplet_dir::String ,
            input_file::String ,
            chemical_potential::Float64 , # fermi level
            hop_symparams::Dict{ String , Matrix{ComplexF64} } ;
            z::Float64=0.0 ,
            max_spin2::Int64=10 ,
            channels_dos::Dict{ String , Vector{Function} }=Dict{ String , Vector{Function} }() ,
            enforce_particle_hole_symmetry::Bool=true,
            minmult::Int64=0 ,
            mine::Float64=0.0 ,
            K_factor::Float64=2.0 ,
            etafac::Float64=0.5 ,
            broadening_distribution::String="gaussian",
            band_width::Float64=1.0 ,
            Nz::Int64=4 ,
            convergence_tolerance::Float64=convergence_tolerance ,
            scale_asymptotic::Bool=true ,
            distributed::Bool=false ,
            input_noninteracting::String="" ,
            compute_selfenergy::Bool=false ,
            u::Float64=0.0 ) where {R<:Real}

    println( "---------------" )
    println( "*    DMFT     *")
    println( "---------------" )

    Z = NRGCalculator.generate_Z(Nz)
    etafac /= Nz
    if distributed
        @everywhere label=$label
        @everywhere L=$L
        @everywhere iterations=$iterations
        @everywhere cutoff_type=$cutoff_type
        @everywhere cutoff_magnitude=$cutoff_magnitude
        @everywhere cg_o_dir=$cg_o_dir
        @everywhere multiplet_dir=$multiplet_dir
        @everywhere input_file=$input_file
        @everywhere K_factor=$K_factor
        @everywhere etafac=$etafac
        @everywhere broadening_distribution=$broadening_distribution
    end

    convergence = false
    iteration = 0
    while !convergence

        iteration += 1
        println( "====================" )
        println( "DMFT iteration $iteration" )
        println( "====================" )

        # weiss nrg
        if !compute_selfenergy
            for z in Z
                PointGroupNRG.NRGCalculator.nrg_molecule( 
                        label,
                        "IMP",
                        L,
                        iterations,
                        cutoff_type,
                        cutoff_magnitude,
                        cg_o_dir,
                        multiplet_dir,
                        input_noninteracting,
                        0.0,
                        hop_symparams;
                        z=z,
                        channels_dos=channels_dos,
                        spectral=true,
                        K_factor=K_factor,
                        spectral_broadening=etafac,
                        broadening_distribution=broadening_distribution,
                        scale_asymptotic=true
                )
            end
            NRGCalculator.zavg_spectral(label,Z;orbitalresolved_number=1)

            # initial LDOS and Weiss field
            ldos_initial = DOS("spectral/spectral_$(label)_zavg_o1.dat")
            ldos_initial = bound_energies(ldos_initial)
            W_initial = Green( chemical_potential , ldos_initial )
            filesave( "green/$(iteration)_ldos_initial.dat" , ldos_initial )
            filesave( "green/$(iteration)_weiss_initial.dat" , W_initial )

        end

        # interacting nrg
        if distributed
            @everywhere hop_symparams=$hop_symparams
            @everywhere channels_dos=$channels_dos
            @sync @distributed for z in Z
                NRGCalculator.nrg_molecule( 
                    label,
                    "IMP",
                    L,
                    iterations,
                    cutoff_type,
                    cutoff_magnitude,
                    cg_o_dir,
                    multiplet_dir,
                    input_file,
                    0.0,
                    hop_symparams;
                    z=z,
                    channels_dos=channels_dos,
                    spectral=true,
                    K_factor=K_factor,
                    spectral_broadening=etafac,
                    broadening_distribution=broadening_distribution,
                    scale_asymptotic=true
                )
            end
        else
            for z in Z
                PointGroupNRG.NRGCalculator.nrg_molecule( 
                    label,
                    "IMP",
                    L,
                    iterations,
                    cutoff_type,
                    cutoff_magnitude,
                    cg_o_dir,
                    multiplet_dir,
                    input_file,
                    0.0,
                    hop_symparams;
                    z=z,
                    channels_dos=channels_dos,
                    spectral=true,
                    K_factor=K_factor,
                    spectral_broadening=etafac,
                    broadening_distribution=broadening_distribution,
                    scale_asymptotic=true,
                    compute_selfenergy=compute_selfenergy
                )
            end
        end
        NRGCalculator.zavg_spectral(label,Z;orbitalresolved_number=1)

        # spectral function from nrg
        A = DOS( "spectral/spectral_$(label)_zavg_o1.dat" )
        A = bound_energies(A)
        #A = bound_energies(A)
        filesave( "green/$(iteration)_spectral.dat" , A )

        # compute initial hybridization and DOS for NRG energy values (for the record)
        band_dos_function = channels_dos["A1g"][1]
        band_dos = DOS( band_dos_function ; x=A.x )
        filesave( "green/$(iteration)_band_dos_initial.dat" , band_dos )
        hybridization_function_initial(x) =  pi * hop_symparams["A1g"][1,1]^2 * channels_dos["A1g"][1](x)
        hybridization = Hybridization(hybridization_function_initial;x=A.x)
        filesave( "green/$(iteration)_hybridization_initial.dat" , hybridization )

        # NRG Green function
        G_NRG = Green( 0.0 , A ) 
        ldos_nrg = DOS( G_NRG )
        filesave( "green/$(iteration)_green_nrg.dat" , G_NRG )
        filesave( "green/$(iteration)_ldos_nrg.dat" , ldos_nrg )

        # NRG Self-energy
        if compute_selfenergy

            # analytic initial weiss field
            W_initial = Green( chemical_potential , hybridization_function_initial , A.x )
            ldos_initial = DOS(W_initial)
            filesave( "green/$(iteration)_weiss_initial.dat" , W_initial )
            filesave( "green/$(iteration)_ldos_initial.dat" , ldos_initial )

            # triple excitation green function
            ldos_triple = DOS("spectral/spectral_$(label)_z0.0_triple.dat")
            #ldos_triple = bound_energies(ldos_triple)
            T_NRG = Green( 0.0 , ldos_triple )
            filesave( "green/$(iteration)_ldos_triple.dat" , ldos_triple )

            # self-energy trick
            SE_NRG = SelfEnergy( T_NRG , G_NRG , u )


        else # initial weiss field already computed

            SE_NRG = SelfEnergy( W_initial , G_NRG )

        end
        filesave( "green/$(iteration)_selfenergy.dat" , SE_NRG )

        # extended interacting Green function
        G_extended = Green( 0.0 , ldos_initial ; se=SE_NRG )
        ldos_extended = DOS(G_extended)
        filesave( "green/$(iteration)_green_extended.dat" , G_extended )
        filesave( "green/$(iteration)_ldos_extended.dat" , ldos_extended )

        # final Weiss field
        W_final = Green( G_extended , -SE_NRG )
        filesave( "green/$(iteration)_weiss_final.dat" , W_final )
        ldos_final = DOS(W_final)
        filesave( "green/$(iteration)_ldos_final.dat" , ldos_final )

        # check convergence
        difference = diff_relative( W_initial , W_final )
        println( "Difference: $difference , Tolerance: $convergence_tolerance" )
        convergence = difference<convergence_tolerance
        if convergence
            filesave( "green/solution_green_extended.dat" , G_extended )
            filesave( "green/solution_weiss.dat" , W_final )
            filesave( "green/solution_selfenergy.dat" , SE_NRG )
            filesave( "green/solution_ldos.dat" , ldos_extended )
            filesave( "green/solution_hybridization.dat" , hybridization )
            println( "CONVERGENCE ACHIEVED IN ITERATION $iteration" )
            return
        end
        sleep(2)

        # extract hybridization function
        hybridization = Hybridization(W_final)
        filesave( "green/$(iteration)_hybridization_final.dat" , hybridization )
        if any(hybridization.y .< 0)
            error( "Negative hybridization at iteration $(iteration). Aborting...")
        end
        hybridization_function_weiss = functionalize(hybridization)
        hybridization_function_final(x) = abs(x)<=1 ? hybridization_function_weiss(x) : 0.0

        # isolate hopping from dos
        piV2 = real(integrate_dmft( hybridization_function_final , (W_final.x[1],W_final.x[end]) ))
        dos_function(x) = hybridization_function_final(x)/piV2
        band_dos = DOS(dos_function;x=W_final.x)
        filesave( "green/$(iteration)_band_dos_final.dat" , band_dos )
        channels_dos["A1g"][1] = dos_function
        # V
        hop_symparams["A1g"][1,1] = sqrt(piV2/pi)

    end


end
