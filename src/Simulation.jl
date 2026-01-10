module Simulation

using LinearAlgebra
using Base.Threads
using ..LinearModel: ModelParams, coeff_matrix

export integration!, run_full_model!

"""
    integration!(state, t, k, init; params, mode=:full)

Integrate the linear system in spectral space for all wavenumbers.

- `state` : ComplexF64 array (Nt, Nv, Nk)
- `t`     : vector of time
- `k`     : vector of wavenumbers
- `init`  : (Nv, Nk) initial state at t[1]
"""
function integration!(state::Array{ComplexF64,3},
                      t::AbstractVector{<:Real},
                      k::AbstractVector{<:Real},
                      init::AbstractMatrix{ComplexF64};
                      params::ModelParams,
                      mode::Symbol = :full)

    Nt, Nv, Nk = size(state)
    Δt = t[2] - t[1]

    @assert Nt == length(t)
    @assert Nv == size(init, 1)
    @assert Nk == length(k)

    # Precompute exponentials for each k
    Φs = Vector{Matrix{ComplexF64}}(undef, Nk)
    for j in eachindex(k)
        L = coeff_matrix(Float64(k[j]); param=params, mode=mode)
        Φs[j] = exp(Δt * L)
    end

    @threads for j in eachindex(k)
        Φ = Φs[j]
        @views state[1, :, j] .= init[:, j]
        tmp = similar(state, ComplexF64, Nv)
        for n in 2:Nt
            @views mul!(tmp, Φ, state[n-1, :, j])
            @views state[n, :, j] .= tmp
        end
    end

    return state
end
end