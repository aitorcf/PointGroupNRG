import Base.getindex
import Base.show

abstract type ClebschGordanSum end

function println( cgs::CGS ) where {CGS<:ClebschGordanSum}
    println("***")
    println(CGS)
    println("***")
    println()
    println("Orbital part")
    println()
    for (irrep_combination,values) in cgs.orbitalsum
        @show irrep_combination
        @show values
        println()
    end
    println("Spin part")
    println()
    for (irrep_combination,value) in cgs.spinsum
        @show irrep_combination
        @show value
        println()
    end
end

# general empty constructor
function EmptyClebschGordanSum(CGS)
    return CGS(Dict{NTuple{3,Int64},Array{ComplexF64,4}}(),Dict{NTuple{3,Int64},Array{ComplexF64,3}}())
end

function ClebschGordanSums( number_of_orbital_irreps::Int64 ,
                            orbital_irreps_impurity_shell::Vector{Int64} ,
                            orbital_irreps_one_electron::Vector{Int64} ,
                            cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}} ,
                            maximum_spin2::Int64 ,
                            maximum_spin2_impurity_shell::Int64 ,
                            cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} ;
                            verbose=false )
    dsum = DSum( number_of_orbital_irreps,
                 orbital_irreps_impurity_shell,
                 orbital_irreps_one_electron,
                 cg_o_fullmatint,
                 maximum_spin2,
                 maximum_spin2_impurity_shell,
                 cg_s_fullmatint;
                 verbose=false )
    ksum = KSum( number_of_orbital_irreps,
                 orbital_irreps_impurity_shell,
                 orbital_irreps_one_electron,
                 cg_o_fullmatint,
                 maximum_spin2,
                 maximum_spin2_impurity_shell,
                 cg_s_fullmatint;
                 verbose=false )
    fsum = FSum( number_of_orbital_irreps,
                 orbital_irreps_impurity_shell,
                 orbital_irreps_one_electron,
                 cg_o_fullmatint,
                 maximum_spin2,
                 maximum_spin2_impurity_shell,
                 cg_s_fullmatint;
                 verbose=false )

    return dsum,ksum,fsum
end

# ----- #
# D SUM #
# ----- #
#
# The D sum appears in the calculation of matrix elements between
# new multiplets:
#
#       ⟨u,αu||H||v,αv⟩ ~ ∑_αβ ∑_a ⟨i||f†_a||j⟩_α × ⟨ν||f†a||μ⟩_β^* × D + h.c.
#
# D[Γuv,Γa,Γi,Γj,Γμ,Γν] = [ αu , αv , α , β ]
struct DSum <: ClebschGordanSum
    orbitalsum::Dict{ NTuple{6,Int64} , Array{ComplexF64,4} }
    spinsum::Dict{ NTuple{6,Int64} , ComplexF64 }
end

function getindex( cgs::DSum , k::NTuple{6,IntIrrep} )::Tuple{ComplexF64,Array{ComplexF64,4}}

    # orbital and spin irreps
    k_orbital::NTuple{6,Int64} = map( x->x[2] , k )
    k_spin::NTuple{6,Int64}    = map( x->x[3] , k )

    # zero elements
    if (!haskey(cgs.orbitalsum,k_orbital) || !haskey(cgs.spinsum,k_spin))
        return (zero(ComplexF64),zeros(ComplexF64,0,0,0,0))
    end

    # non-zero elements
    array_orbital::Array{ComplexF64,4} = cgs.orbitalsum[k_orbital]
    sign_and_spin::ComplexF64 = (-1)^(k[5][1])*cgs.spinsum[k_spin]

    return (sign_and_spin,array_orbital)
end

function DSum( number_of_orbital_irreps::Int64 ,
               orbital_irreps_impurity_shell::Vector{Int64} ,
               orbital_irreps_one_electron::Vector{Int64} ,
               cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}} ,
               maximum_spin2::Int64 ,
               maximum_spin2_impurity_shell::Int64 ,
               cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} ;
               verbose=false )

    # ORBITAL SECTOR
    #
    orbitalsum = Dict{ NTuple{6,Int64} , Array{ComplexF64,4} }()
    i_uv::Int64 = 1
    # iterate through irreps
    for I_uv in 1:number_of_orbital_irreps,
        I_a  in orbital_irreps_one_electron,
        I_i  in 1:number_of_orbital_irreps,
        I_j  in 1:number_of_orbital_irreps,
        I_mu in orbital_irreps_impurity_shell,
        I_nu in orbital_irreps_impurity_shell

        # early discard
        haskey( cg_o_fullmatint , (I_mu,I_i, I_uv) ) || continue
        haskey( cg_o_fullmatint , (I_nu,I_j, I_uv) ) || continue
        haskey( cg_o_fullmatint , (I_a, I_j, I_i)  ) || continue
        haskey( cg_o_fullmatint , (I_a, I_mu,I_nu) ) || continue

        # clebsch-gordan array views
        @views begin
            cg_o_muiuv = cg_o_fullmatint[I_mu,I_i, I_uv][:,:,:,:]
            cg_o_nujuv = cg_o_fullmatint[I_nu,I_j, I_uv][:,:,:,:]
            cg_o_aji   = cg_o_fullmatint[I_a, I_j, I_i][:,:,:,:]
            cg_o_amunu = cg_o_fullmatint[I_a, I_mu,I_nu][:,:,:,:]
        end

        # construct dictionary entry
        orbitalarray = zeros(ComplexF64,size(cg_o_muiuv,1),size(cg_o_nujuv,1),size(cg_o_aji,1),size(cg_o_amunu,1))

        # iterate through outer multiplicities
        for αu in axes(orbitalarray,1),
            αv in axes(orbitalarray,2),
            α  in axes(orbitalarray,3),
            β  in axes(orbitalarray,4)

            # array entry
            element::ComplexF64 = 0.0

            # iterate through partners
            for i_a  in axes(cg_o_aji,2),
                i_i  in axes(cg_o_aji,4),
                i_j  in axes(cg_o_aji,3),
                i_mu in axes(cg_o_amunu,3),
                i_nu in axes(cg_o_amunu,4)

                element += conj(cg_o_muiuv[αu,i_mu,i_i,i_uv])*cg_o_nujuv[αv,i_nu,i_j,i_uv]*
                           conj(cg_o_aji[α,i_a,i_j,i_i])*     cg_o_amunu[β,i_a,i_mu,i_nu]

            end

            orbitalarray[αu,αv,α,β] = element

        end

        all(iszero.(orbitalarray)) && continue
        orbitalsum[I_uv,I_a,I_i,I_j,I_mu,I_nu] = orbitalarray

    end

    # SPIN SECTOR
    #
    maximum_spin2_onelectron = iszero(maximum_spin2) ? 0 : 1
    spinsum = Dict{ NTuple{6,Int64} , ComplexF64 }()
    si_uv::Int64 = 1
    # iterate through irreps
    for S_uv in 0:maximum_spin2,
        S_a  in 0:maximum_spin2_onelectron,
        S_i  in 0:maximum_spin2,
        S_j  in 0:maximum_spin2,
        S_mu in 0:maximum_spin2_impurity_shell,
        S_nu in 0:maximum_spin2_impurity_shell

        # early discard
        haskey( cg_s_fullmatint , (S_mu,S_i, S_uv) ) || continue
        haskey( cg_s_fullmatint , (S_nu,S_j, S_uv) ) || continue
        haskey( cg_s_fullmatint , (S_a, S_j, S_i)  ) || continue
        haskey( cg_s_fullmatint , (S_a, S_mu,S_nu) ) || continue

        # clebsch-gordan array views
        @views begin
            cg_s_muiuv = cg_s_fullmatint[S_mu,S_i, S_uv][:,:,:]
            cg_s_nujuv = cg_s_fullmatint[S_nu,S_j, S_uv][:,:,:]
            cg_s_aji   = cg_s_fullmatint[S_a, S_j, S_i][:,:,:]
            cg_s_amunu = cg_s_fullmatint[S_a, S_mu,S_nu][:,:,:]
        end

        # dictionary entry
        element::ComplexF64 = 0.0

        # iterate through partners
        for si_a  in axes(cg_s_aji,1),
            si_i  in axes(cg_s_aji,3),
            si_j  in axes(cg_s_aji,2),
            si_mu in axes(cg_s_amunu,2),
            si_nu in axes(cg_s_amunu,3)

            element += conj(cg_s_muiuv[si_mu,si_i,si_uv])*cg_s_nujuv[si_nu,si_j,si_uv]*
                       conj(cg_s_aji[si_a,si_j,si_i])*    cg_s_amunu[si_a,si_mu,si_nu]

        end

        iszero(element) && continue
        spinsum[S_uv,S_a,S_i,S_j,S_mu,S_nu] = element

    end

    return DSum(orbitalsum,spinsum)
end

# ----- #
# K SUM #
# ----- #
#
# Used for constructing new matrix elements of 
# local creation operators:
#
#   ⟨ u , αu || f†_a || v , αv ⟩_β = ∑_α ⟨ μ || f†_a || ν ⟩_α × K(Γu,Γv,Γa,Γij,Γμ,Γν,αu,αv,α,β)
#
# 
struct KSum <: ClebschGordanSum
    orbitalsum::Dict{ NTuple{6,Int64} , Array{ComplexF64,4} }
    spinsum::Dict{ NTuple{6,Int64} , ComplexF64 }
end

function getindex( cgs::KSum , k::NTuple{6,IntIrrep} )::Tuple{ComplexF64,Array{ComplexF64,4}}

    # orbital and spin irreps
    k_orbital::NTuple{6,Int64} = map( x->x[2] , k )
    k_spin::NTuple{6,Int64}    = map( x->x[3] , k )

    # zero elements
    if (!haskey(cgs.orbitalsum,k_orbital) || !haskey(cgs.spinsum,k_spin))
        return (zero(ComplexF64),zeros(ComplexF64,0,0,0,0))
    end

    # non-zero elements
    array_orbital::Array{ComplexF64,4} = cgs.orbitalsum[k_orbital]
    value_spin::ComplexF64 = cgs.spinsum[k_spin]

    return (value_spin,array_orbital)
end

function KSum( number_of_orbital_irreps::Int64 ,
               orbital_irreps_impurity_shell::Vector{Int64} ,
               orbital_irreps_one_electron::Vector{Int64} ,
               cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}} ,
               maximum_spin2::Int64 ,
               maximum_spin2_impurity_shell::Int64 ,
               cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} ;
               verbose=false )

    # ORBITAL SECTOR
    #
    orbitalsum = Dict{ NTuple{6,Int64} , Array{ComplexF64,4} }()
    i_u::Int64 = 1
    # iterate through irreps
    for I_u  in 1:number_of_orbital_irreps,
        I_v  in 1:number_of_orbital_irreps,
        I_a  in orbital_irreps_one_electron,
        I_ij in 1:number_of_orbital_irreps,
        I_mu in orbital_irreps_impurity_shell,
        I_nu in orbital_irreps_impurity_shell

        # early discard
        haskey( cg_o_fullmatint , (I_mu,I_ij,I_u)  ) || continue
        haskey( cg_o_fullmatint , (I_nu,I_ij,I_v)  ) || continue
        haskey( cg_o_fullmatint , (I_a, I_nu,I_mu) ) || continue
        haskey( cg_o_fullmatint , (I_a, I_v, I_u)  ) || continue

        # clebsch-gordan array views
        @views begin
            cg_o_muiju = cg_o_fullmatint[I_mu,I_ij,I_u][:,:,:,:]
            cg_o_nuijv = cg_o_fullmatint[I_nu,I_ij,I_v][:,:,:,:]
            cg_o_anumu = cg_o_fullmatint[I_a, I_nu,I_mu][:,:,:,:]
            cg_o_avu   = cg_o_fullmatint[I_a, I_v, I_u][:,:,:,:]
        end

        # construct dictionary entry
        orbitalarray = zeros(ComplexF64,size(cg_o_muiju,1),size(cg_o_nuijv,1),size(cg_o_anumu,1),size(cg_o_avu,1))

        # iterate through outer multiplicities
        for αu in axes(orbitalarray,1),
            αv in axes(orbitalarray,2),
            α  in axes(orbitalarray,3),
            β  in axes(orbitalarray,4)

            # array entry
            element::ComplexF64 = 0.0

            # iterate through partners
            for i_v  in axes(cg_o_nuijv,4),
                i_a  in axes(cg_o_avu,2),
                i_ij in axes(cg_o_nuijv,3),
                i_mu in axes(cg_o_anumu,4),
                i_nu in axes(cg_o_anumu,3)

                element += conj(cg_o_muiju[αu,i_mu,i_ij,i_u])*cg_o_nuijv[αv,i_nu,i_ij,i_v]*
                           conj(cg_o_anumu[α,i_a,i_nu,i_mu])* cg_o_avu[β,i_a,i_v,i_u]

            end

            orbitalarray[αu,αv,α,β] = element

        end

        all(iszero.(orbitalarray)) && continue
        orbitalsum[I_u,I_v,I_a,I_ij,I_mu,I_nu] = orbitalarray

    end

    # SPIN SECTOR
    #
    maximum_spin2_onelectron = iszero(maximum_spin2) ? 0 : 1
    spinsum = Dict{ NTuple{6,Int64} , ComplexF64 }()
    si_u::Int64 = 1
    # iterate through irreps
    for S_u  in 0:maximum_spin2,
        S_v  in 0:maximum_spin2,
        S_a  in 0:maximum_spin2_onelectron,
        S_ij in 0:maximum_spin2,
        S_mu in 0:maximum_spin2_impurity_shell,
        S_nu in 0:maximum_spin2_impurity_shell

        # early discard
        haskey( cg_s_fullmatint , (S_mu,S_ij,S_u)  ) || continue
        haskey( cg_s_fullmatint , (S_nu,S_ij,S_v)  ) || continue
        haskey( cg_s_fullmatint , (S_a, S_nu,S_mu) ) || continue
        haskey( cg_s_fullmatint , (S_a, S_v, S_u)  ) || continue

        # clebsch-gordan array views
        @views begin
            cg_s_muiju = cg_s_fullmatint[S_mu,S_ij,S_u][:,:,:]
            cg_s_nuijv = cg_s_fullmatint[S_nu,S_ij,S_v][:,:,:]
            cg_s_anumu = cg_s_fullmatint[S_a, S_nu,S_mu][:,:,:]
            cg_s_avu   = cg_s_fullmatint[S_a, S_v, S_u][:,:,:]
        end

        # dictionary entry
        element::ComplexF64 = 0.0

        # iterate through partners
        for si_v  in axes(cg_s_nuijv,3),
            si_a  in axes(cg_s_anumu,1),
            si_ij in axes(cg_s_muiju,2),
            si_mu in axes(cg_s_muiju,1),
            si_nu in axes(cg_s_nuijv,1)

            element += conj(cg_s_muiju[si_mu,si_ij,si_u])*cg_s_nuijv[si_nu,si_ij,si_v]*
                       conj(cg_s_anumu[si_a,si_nu,si_mu])*cg_s_avu[si_a,si_v,si_u]

        end

        iszero(element) && continue
        spinsum[S_u,S_v,S_a,S_ij,S_mu,S_nu] = element

    end

    return KSum(orbitalsum,spinsum)
end

# ----- #
# F SUM #
# ----- #
#
# Used to construct matrix elements of block operators
#
#   ⟨u||f†_a||v⟩_β = ∑_α ⟨i||f†_a||j⟩_α × δ_μν × F(Γu,Γv,Γa,Γi,Γj,Γμν)
#
struct FSum <: ClebschGordanSum
    orbitalsum::Dict{ NTuple{6,Int64} , Array{ComplexF64,4} }
    spinsum::Dict{ NTuple{6,Int64} , ComplexF64 }
end

function getindex( cgs::FSum , k::NTuple{6,IntIrrep} )::Tuple{ComplexF64,Array{ComplexF64,4}}

    # orbital and spin irreps
    k_orbital::NTuple{6,Int64} = map( x->x[2] , k )
    k_spin::NTuple{6,Int64}    = map( x->x[3] , k )

    # zero elements
    if (!haskey(cgs.orbitalsum,k_orbital) || !haskey(cgs.spinsum,k_spin))
        return (zero(ComplexF64),zeros(ComplexF64,0,0,0,0))
    end

    # non-zero elements
    sign::ComplexF64 = (-1)^k[6][1]
    array_orbital::Array{ComplexF64,4} = cgs.orbitalsum[k_orbital]
    value_spin::ComplexF64 = cgs.spinsum[k_spin]

    return (sign*value_spin,array_orbital)
end

function FSum( number_of_orbital_irreps::Int64 ,
               orbital_irreps_impurity_shell::Vector{Int64} ,
               orbital_irreps_one_electron::Vector{Int64} ,
               cg_o_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,4}} ,
               maximum_spin2::Int64 ,
               maximum_spin2_impurity_shell::Int64 ,
               cg_s_fullmatint::Dict{NTuple{3,Int64},Array{ComplexF64,3}} ;
               verbose=false )

    # ORBITAL SECTOR
    #
    orbitalsum = Dict{ NTuple{6,Int64} , Array{ComplexF64,4} }()
    i_u::Int64 = 1
    # iterate through irreps
    for I_u    in 1:number_of_orbital_irreps,
        I_v    in 1:number_of_orbital_irreps,
        I_a    in orbital_irreps_one_electron,
        I_i    in 1:number_of_orbital_irreps,
        I_j    in 1:number_of_orbital_irreps,
        I_munu in orbital_irreps_impurity_shell

        # early discard
        haskey( cg_o_fullmatint , (I_munu,I_i,I_u) ) || continue
        haskey( cg_o_fullmatint , (I_munu,I_j,I_v) ) || continue
        haskey( cg_o_fullmatint , (I_a,   I_j,I_i) ) || continue
        haskey( cg_o_fullmatint , (I_a,   I_v,I_u) ) || continue

        # clebsch-gordan array views
        @views begin
            cg_o_munuiu = cg_o_fullmatint[I_munu,I_i,I_u][:,:,:,:]
            cg_o_munujv = cg_o_fullmatint[I_munu,I_j,I_v][:,:,:,:]
            cg_o_aji    = cg_o_fullmatint[I_a,   I_j,I_i][:,:,:,:]
            cg_o_avu    = cg_o_fullmatint[I_a,   I_v,I_u][:,:,:,:]
        end

        # construct dictionary entry
        orbitalarray = zeros(ComplexF64,size(cg_o_munuiu,1),size(cg_o_munujv,1),size(cg_o_aji,1),size(cg_o_avu,1))

        # iterate through outer multiplicities
        for αu in axes(orbitalarray,1),
            αv in axes(orbitalarray,2),
            α  in axes(orbitalarray,3),
            β  in axes(orbitalarray,4)

            # array entry
            element::ComplexF64 = 0.0

            # iterate through partners
            for i_v    in axes(cg_o_munujv,4),
                i_a    in axes(cg_o_aji,2),
                i_i    in axes(cg_o_aji,4),
                i_j    in axes(cg_o_aji,3),
                i_munu in axes(cg_o_munuiu,2)

                element += conj(cg_o_munuiu[αu,i_munu,i_i,i_u])*cg_o_munujv[αv,i_munu,i_j,i_v]*
                           conj(cg_o_aji[α,i_a,i_j,i_i])*       cg_o_avu[β,i_a,i_v,i_u]

            end

            orbitalarray[αu,αv,α,β] = element

        end

        all(iszero.(orbitalarray)) && continue
        orbitalsum[I_u,I_v,I_a,I_i,I_j,I_munu] = orbitalarray

    end

    # SPIN SECTOR
    #
    maximum_spin2_onelectron = iszero(maximum_spin2) ? 0 : 1
    spinsum = Dict{ NTuple{6,Int64} , ComplexF64 }()
    si_u::Int64 = 1
    # iterate through irreps
    for S_u    in 0:maximum_spin2,
        S_v    in 0:maximum_spin2,
        S_a    in 0:maximum_spin2_onelectron,
        S_i    in 0:maximum_spin2,
        S_j    in 0:maximum_spin2,
        S_munu in 0:maximum_spin2_impurity_shell

        # early discard
        haskey( cg_s_fullmatint , (S_munu,S_i,S_u) ) || continue
        haskey( cg_s_fullmatint , (S_munu,S_j,S_v) ) || continue
        haskey( cg_s_fullmatint , (S_a,   S_j,S_i) ) || continue
        haskey( cg_s_fullmatint , (S_a,   S_v,S_u) ) || continue

        # clebsch-gordan array views
        @views begin
            cg_s_munuiu = cg_s_fullmatint[S_munu,S_i,S_u][:,:,:]
            cg_s_munujv = cg_s_fullmatint[S_munu,S_j,S_v][:,:,:]
            cg_s_aji    = cg_s_fullmatint[S_a,   S_j,S_i][:,:,:]
            cg_s_avu    = cg_s_fullmatint[S_a,   S_v,S_u][:,:,:]
        end

        # dictionary entry
        element::ComplexF64 = 0.0

        # iterate through partners
        for si_v    in axes(cg_s_munujv,3),
            si_a    in axes(cg_s_aji,1),
            si_i    in axes(cg_s_aji,3),
            si_j    in axes(cg_s_aji,2),
            si_munu in axes(cg_s_munujv,1)

            element += conj(cg_s_munuiu[si_munu,si_i,si_u])*cg_s_munujv[si_munu,si_j,si_v]*
                       conj(cg_s_aji[si_a,si_j,si_i])*      cg_s_avu[si_a,si_v,si_u]

        end

        iszero(element) && continue
        spinsum[S_u,S_v,S_a,S_i,S_j,S_munu] = element

    end

    return FSum(orbitalsum,spinsum)
end
