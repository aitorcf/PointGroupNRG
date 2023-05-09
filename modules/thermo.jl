# ##########################################################
# This module contains functions that compute thermodynamic
# averages and other related quantities with the information 
# of a given NRG step. 
# ##########################################################

include("discretization.jl")
using DelimitedFiles

function compute_temperature(N::Int64,
    L::Float64,
    betabar::Float64;
    z::Float64=0.0,
    discretization::String="standard",
    ebar0::Float64=1.0)

    if discretization == "standard"

        return 0.5 * (1 + L^(-1)) * L^(z - (N - 2) / 2) / (betabar)


    elseif discretization == "co2005"

        eco0_z = compute_epsilon_co2005(8, z, L)[9] * L^4
        return eco0_z * L^(-(N - 2) / 2) / betabar

    elseif discretization == "lanczos"

        return ebar0 * L^(-(N - 2) / 2) / betabar

    end
end

function compute_partition_function(
    irrEU::Dict{NTuple{3,Int64},Tuple{Vector{Float64},Matrix{ComplexF64}}},
    betabar::Float64,
    oirreps2dimensions::Dict)::Float64

    #   Z = Tr exp( - betabar * H_N )
    #     = Sum_G( D * ( Sum_i exp(-betabar*E_i) ))
    #

    part::Float64 = 0.0
    for (G, (E, U)) in irrEU
        # orbital and spin
        (I, S) = G[2:3]

        # degeneracy/dimensionality
        Do = oirreps2dimensions[I]
        Ds = 2 * S + 1
        D = Do * Ds


        # weighted sum over G 
        part += D * sum(exp(-betabar * e) for e in E)
    end
    return part
end

# int method
function compute_partition_function(
    irrEU::Dict{NTuple{3,Int64},Tuple{Vector{Float64},Matrix{ComplexF64}}},
    betabar::Float64,
    oindex2dimensions::Vector{Int64})::Float64
    #   Z = Tr exp( - betabar * H_N )
    #     = Sum_G( D * ( Sum_i exp(-betabar*E_i) ))

    part::Float64 = 0.0
    for (G, (E, U)) in irrEU
        # orbital and spin
        (I, S) = G[2:3]


        # degeneracy/dimensionality
        Do = oindex2dimensions[I]::Int64
        Ds = S + 1
        D = Do * Ds

        # weighted sum over G 
        part += D * sum(exp(-betabar * e) for e in E)
    end
    return part
end

function compute_magnetic_susceptibility(
            irrEU, 
            betabar, 
            oirreps2dimensions::Dict; 
            verbose=false)

    verbose && println("MAGNETIC SUSCEPTIBILITY CALCULATION")

    part = compute_partition_function(irrEU, betabar, oirreps2dimensions)

    mag = 0
    for (G, (E, U)) in irrEU
        # orbital and spin
        (I, S) = G[2:3]

        # degeneracy/dimensionality
        Do = oirreps2dimensions[I]
        Ds = 2 * S + 1
        D = Do * Ds

        # weighted sum over G 
        # S^2 = S(S+1)
        # Sz^2 = S^2/3
        contrib = D * (S * (S + 1) / 3.0) * sum(exp(-betabar * e) for e in E)
        mag += contrib

        if verbose
            println("G = $G, D = $D, E = $E")
            println("contribution = $(contrib/part)")
            println()
        end
    end
    return mag / part
end

# int method
function compute_magnetic_susceptibility(
    irrEU::Dict{NTuple{3,Int64},Tuple{Vector{Float64},Matrix{ComplexF64}}},
    betabar::Float64,
    oindex2dimensions::Vector{Int64};
    verbose=false)

    verbose && println("MAGNETIC SUSCEPTIBILITY CALCULATION")

    part::Float64 = compute_partition_function(irrEU, betabar, oindex2dimensions)

    mag::Float64 = 0
    for (G, (E, U)) in irrEU
        # orbital and spin
        (I, S) = G[2:3]

        # degeneracy/dimensionality
        Do = oindex2dimensions[I]::Int64
        Ds = S + 1
        D = Do * Ds

        # weighted sum over G 
        # S^2 = S(S+1)
        # Sz^2 = S^2/3
        contrib = D * (S / 2.0 * (S / 2.0 + 1) / 3.0) * sum(exp(-betabar * e) for e in E)
        mag += contrib

        if verbose
            println("G = $G, D = $D, E = $E")
            println("contribution = $(contrib/part)")
            println()
        end
    end
    return mag / part
end

# string method
function compute_average_number_of_particles(
        irrEU, 
        betabar, 
        oirreps2dimensions::Dict; 
        verbose=false)

    part = compute_partition_function(irrEU, betabar, oirreps2dimensions)

    N = 0
    for (G, (E, U)) in irrEU
        # orbital and spin
        (N, I, S) = G

        # degeneracy/dimensionality
        Do = oirreps2dimensions[I]
        Ds = 2 * S + 1
        D = Do * Ds

        # weighted sum over G 
        contrib = D * N * sum(exp(-betabar * e) for e in E)
        N += contrib
    end
    return N / part
end

# int method 
function compute_average_number_of_particles(
        irrEU::Dict{NTuple{3,Int64},Tuple{Vector{Float64},Matrix{ComplexF64}}},
        betabar::Float64,
        oindex2dimensions::Vector{Int64};
        verbose=false)

    part::Float64 = compute_partition_function(irrEU, betabar, oindex2dimensions)

    N::Float64 = 0.0
    for (G, (E, U)) in irrEU
        # orbital and spin
        (N, I, S) = G

        # degeneracy/dimensionality
        Do = oindex2dimensions[I]
        Ds = S + 1
        D = Do * Ds

        # weighted sum over G 
        contrib = D * N * sum(exp(-betabar * e) for e in E)
        N += contrib
    end
    return N / part
end

# string method
function compute_energy(irrEU, betabar, oirreps2dimensions::Dict; verbose=false)

    part = compute_partition_function(irrEU, betabar, oirreps2dimensions)

    energy = 0
    for (G, (E, U)) in irrEU
        # orbital and spin
        (N, I, S) = G

        # degeneracy/dimensionality
        Do = oirreps2dimensions[I]
        Ds = 2 * S + 1
        D = Do * Ds

        # weighted sum over G 
        contrib = D * sum((e * exp(-betabar * e)) for e in E)
        energy += contrib
    end
    return energy / part
end

# int method
function compute_energy(
    irrEU::Dict{NTuple{3,Int64},Tuple{Vector{Float64},Matrix{ComplexF64}}},
    betabar::Float64,
    oindex2dimensions::Vector{Int64};
    verbose=false)

    part::Float64 = compute_partition_function(irrEU, betabar, oindex2dimensions)

    energy::Float64 = 0
    for (G, (E, U)) in irrEU
        # orbital and spin
        (N, I, S) = G

        # degeneracy/dimensionality
        Do = oindex2dimensions[I]
        Ds = S + 1
        D = Do * Ds

        # weighted sum over G 
        contrib = D * sum((e * exp(-betabar * e)) for e in E)
        energy += contrib
    end
    return energy / part
end
# int method
function energy2(
    irrEU::Dict{NTuple{3,Int64},Tuple{Vector{Float64},Matrix{ComplexF64}}},
    betabar::Float64,
    oindex2dimensions::Vector{Int64};
    verbose=false)

    part::Float64 = compute_partition_function(irrEU, betabar, oindex2dimensions)

    energy::Float64 = 0
    for (G, (E, U)) in irrEU
        # orbital and spin
        (N, I, S) = G

        # degeneracy/dimensionality
        Do = oindex2dimensions[I]
        Ds = S + 1
        D = Do * Ds

        # weighted sum over G 
        contrib = D * sum((e * e * exp(-betabar * e)) for e in E)
        energy += contrib
    end
    return energy / part
end

# int method 
function compute_entropy(
    irrEU::Dict{NTuple{3,Int64},Tuple{Vector{Float64},Matrix{ComplexF64}}},
    betabar::Float64,
    oindex2dimensions::Vector{Int64};
    verbose=false)

    part::Float64 = compute_partition_function(irrEU, betabar, oindex2dimensions)
    e::Float64 = compute_energy(irrEU, betabar, oindex2dimensions)
    return log(part) + betabar * e

end

# int method 
function compute_heat_capacity(
    irrEU::Dict{NTuple{3,Int64},Tuple{Vector{Float64},Matrix{ComplexF64}}},
    betabar::Float64,
    oindex2dimensions::Vector{Int64};
    verbose=false)

    e = compute_energy(irrEU, betabar, oindex2dimensions)
    e2 = energy2(irrEU, betabar, oindex2dimensions)
    return betabar * betabar * (e2 - e * e)

end

# int method 
function compute_free_energy(
    irrEU::Dict{NTuple{3,Int64},Tuple{Vector{Float64},Matrix{ComplexF64}}},
    betabar::Float64,
    oindex2dimensions::Vector{Int64};
    verbose=false)

    return log(compute_partition_function(irrEU, betabar, oindex2dimensions))

end

# ----------------------- #
# IMPURITY THERMODYNAMICS #
# ----------------------- #
function imp_qnums(
    irrEU::Dict{NTuple{3,Int64},Tuple{Vector{Float64},Matrix{ComplexF64}}},
    oindex2dimensions::Vector{Int64},
    combinations_uprima::Dict{NTuple{3,Int64},Vector{NTuple{3,NTuple{4,Int64}}}},
    ss_ip,
    nn_ip)
    # compute averages of impurity quantum numbers
    # in order to keep track of the atomic low-energy
    # states

    # average impurity spin^2 and N per multiplet
    ss_u = Dict()
    nn_u = Dict()

    # iterate over irrep blocks
    for (G, (E, U)) in irrEU

        # iterate over multiplets in irrep block
        # (store the block-shell combination)
        for (m_u, m_mu, m_i) in combinations_uprima[G]

            # qnum setup
            ss_u[m_u] = 0.0
            nn_u[m_u] = 0.0

            # outer multiplicity
            r_u = m_u[4]

            ## iterate over primed multiplets
            for c_up in combinations_uprima[G]

                # U(Γ) term
                m_up = c_up[1]
                r_up = m_up[4]
                uterm = abs2(U[r_up, r_u])

                # primed term 
                m_ip = c_up[3]

                # quantum number averages
                ss_u[m_u] += uterm * ss_ip[m_ip]
                nn_u[m_u] += uterm * nn_ip[m_ip]

            end
        end
    end

    return (ss_u, nn_u)
end

function impspin(
    irrEU::Dict{NTuple{3,Int64},Tuple{Vector{Float64},Matrix{ComplexF64}}},
    betabar::Float64,
    oindex2dimensions::Vector{Int64},
    ss_u;
    verbose=false)

    s_imp = 0.0

    part::Float64 = compute_partition_function(irrEU, betabar, oindex2dimensions)

    energy = 0
    for (G, (E, U)) in irrEU

        # number of multiplets in Γ
        R_G = length(E)

        # orbital and spin
        (N, I, S) = G

        # degeneracy/dimensionality
        Do = oindex2dimensions[I]
        Ds = S + 1
        D = Do * Ds

        # weighted sum over Γ
        for r_g in 1:R_G
            m_u = (G..., r_g)
            s_imp += D * ss_u[m_u] * exp(-betabar * E[r_g])
        end

    end
    return s_imp / part
end

function impnum(
    irrEU::Dict{NTuple{3,Int64},Tuple{Vector{Float64},Matrix{ComplexF64}}},
    betabar::Float64,
    oindex2dimensions::Vector{Int64},
    nn_u;
    verbose=false)

    n_imp = 0.0

    part::Float64 = compute_partition_function(irrEU, betabar, oindex2dimensions)

    energy = 0
    for (G, (E, U)) in irrEU

        # number of multiplets in Γ
        R_G = length(E)

        # orbital and spin
        (N, I, S) = G

        # degeneracy/dimensionality
        Do = oindex2dimensions[I]
        Ds = S + 1
        D = Do * Ds

        # weighted sum over Γ
        for r_g in 1:R_G
            m_u = (G..., r_g)
            n_imp += D * nn_u[m_u] * exp(-betabar * E[r_g])
        end

    end
    return n_imp / part
end

# ------------------- #
# IMPURITY MULTIPLETS #
# ------------------- #
function imp_mults(
    irrEU::Dict{NTuple{3,Int64},Tuple{Vector{Float64},Matrix{ComplexF64}}},
    oindex2dimensions::Vector{Int64},
    combinations_uprima::Dict{NTuple{3,Int64},Vector{NTuple{3,NTuple{4,Int64}}}},
    mm_ip::Dict{NTuple{4,Int64},Vector{Float64}})

    # average impurity spin^2 and N per multiplet
    M::Int64 = length(collect(values(mm_ip))[1])
    mm_u::Dict{NTuple{4,Int64},Vector{Float64}} = Dict()

    # iterate over irrep blocks
    for (G, (E, U)) in irrEU

        # iterate over multiplets in irrep block
        # (store the block-shell combination)
        for (m_u, m_mu, m_i) in combinations_uprima[G]

            # qnum temp setup
            tmp_m_u = zeros(Float64,M)

            # outer multiplicity
            r_u = m_u[4]


            ## iterate over primed multiplets
            for c_up in combinations_uprima[G]

                # U(Γ) term
                m_up = c_up[1]
                r_up = m_up[4]
                uterm = abs2(U[r_up, r_u])

                # primed term 
                m_ip::NTuple{4,Int64} = c_up[3]

                # quantum number averages
                tmp_m_u .+= uterm .* mm_ip[m_ip]

            end
            mm_u[m_u] = tmp_m_u
        end
    end

    return mm_u
end
function mult_thermo(
    irrEU::Dict{NTuple{3,Int64},Tuple{Vector{Float64},Matrix{ComplexF64}}},
    betabar::Float64,
    oindex2dimensions::Vector{Int64},
    mm_u::Dict{NTuple{4,Int64},Vector{Float64}};
    verbose=false)::Vector{Float64}

    M::Int64 = length(collect(values(mm_u))[1])
    mult_imp::Vector{Float64} = [0.0 for i in 1:M]

    part::Float64 = compute_partition_function(irrEU, betabar, oindex2dimensions)

    for (G, (E, U)) in irrEU

        # number of multiplets in Γ
        R_G = length(E)

        # orbital and spin
        (N, I, S) = G

        # degeneracy/dimensionality
        Do = oindex2dimensions[I]
        Ds = S + 1
        D = Do * Ds

        # weighted sum over Γ
        for r_g in 1:R_G
            m_u = (G..., r_g)
            mult_imp .+= D * mm_u[m_u] * exp(-betabar * E[r_g])
        end

    end
    return mult_imp ./ part
end

# ~~~~~~~~~~~~~~~ #
# WRITING TO FILE #
# ~~~~~~~~~~~~~~~ #

function write_impurity_info(
        nrg,
        omults,
        mult2index,
        orbital,
        z)
    multfile = "thermodata/impmult_$(orbital)_z$z.dat"
    open(multfile, write=true) do f
        for m in omults
            write(f, "$(m[1]),$(m[2]),$(m[3]) ")
            for n in 1:length(nrg.impmults)
                write(f, "$(nrg.impmults[end+1-n][mult2index[m]]) ")
            end
            write(f, "\n")
        end
    end
end

function write_thermodata_onez(
        nrg,
        calculation,
        orbital,
        z)

    filename = "thermodata/thermo_$(lowercase(calculation))_$(orbital)_z$z.dat"

    println("Saving thermodynamic data to $filename...\n\n" )

    open( filename , write=true ) do f
        write( f , "# 1 temperature | 2 magnetic susceptibility | 3 entropy | 4 heat capacity | 5 free energy | 6 number of particles | 7 energy | 8 partition function\n" )
        writedlm(
            f ,
            [nrg.t nrg.m nrg.entr nrg.c nrg.f nrg.n nrg.e nrg.p]
        )
    end
end

function write_thermodiff(orbital, z)

    th_clean = readdlm( "thermodata/thermo_clean_$(orbital)_z$z.dat" , skipstart=1 )
    t = th_clean[:, 1]

    th_imp = readdlm( "thermodata/thermo_imp_$(orbital)_z$z.dat" , skipstart=1 )

    if size(th_clean) !== size(th_imp)
        println( "Different number of iterations with respect to clean calculation. Thermodiff will not be computed." ) 
        return nothing
    end

    th_diff = th_imp .- th_clean
    th_diff[:, 4] .= th_imp[:, 4] ./ th_clean[:, 4]
    th_diff[:, 1] .= th_imp[:, 1]

    mfile = "thermodata/th_diff_$(orbital)_z$z.dat"
    writedlm(mfile, th_diff)

end
