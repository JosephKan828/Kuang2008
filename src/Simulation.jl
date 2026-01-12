module Simulation

using LinearAlgebra
using Base.Threads
using ..LinearModel: ModelParams, coeff_matrix

export integration!

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
    # Φs = Vector{Matrix{ComplexF64}}(undef, Nk)
    Φs = [Matrix{ComplexF64}(undef, Nv, Nv) for _ in 1:Nk]

    @threads for j in 1:Nk
        L = coeff_matrix(Float64(k[j]); param=params, mode=mode)
        Φs[j] .= exp(Δt * L)
    end

    # Integrate
    workspaces = [Vector{ComplexF64}(undef, Nv) for _ in 1:nthreads()]

    @threads for j in 1:Nk
        tid = threadid()
        tmp = workspaces[tid]

        Φ = Φs[j]

        @views state[1, :, j] .= init[:, j]
        
        for n in 2:Nt
            prev_slice = view(state, n-1, :, j)

            mul!(tmp, Φ, prev_slice)

            @views state[n, :, j] .= tmp
        end
    end

    return state
end
end