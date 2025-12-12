module Diagnostics

using LinearAlgebra
using HDF5
using Base.Threads

# This assumes Diagnostics.jl is included from within module Kuang2008
# and LinearModel.jl defines module LinearModel.
using ..LinearModel: ModelParams, default_params, coeff_matrix

export diagnose_growth_phase, save_growth_phase

"""
    diagnose_growth_phase(params; λ = 540.0:540.0:43200.0, Nv=6)

Compute growth rate and phase speed for a set of wavelengths `λ`
given model parameters `params`.

Returns `(λ_vec, k_vec, growth, speed)`:

- `λ_vec`   : Vector of wavelengths (km)
- `k_vec`   : Vector of nondimensional wavenumbers
- `growth`  : (Nv, Nk) growth rates (1/day)
- `speed`   : (Nv, Nk) phase speeds (m/s)
"""
function diagnose_growth_phase(params::ModelParams;
                               λ = collect(540.0:540.0:43200.0),
                               Nv::Int = 6)

    # k = 2π * Lref / λ (same as your original formula)
    k = 2π .* 4320.0 ./ λ
    Nk = length(k)

    growth = zeros(Float64, Nv, Nk)
    speed  = zeros(Float64, Nv, Nk)

    @threads for i in eachindex(k)
        mat = coeff_matrix(k[i]; param = params)
        eigval = eigvals(mat)

        σ = real.(eigval)
        c = -imag.(eigval) ./ k[i] * (4320000.0 / 86400.0)

        growth[:, i] .= σ[1:Nv]
        speed[:, i]  .= c[1:Nv]
    end

    return λ, k, growth, speed
end


"""
    save_growth_phase(outfile, λ, k, growth, speed)

Save growth rate and phase speed diagnostics to an HDF5 file
with basic metadata.
"""
function save_growth_phase(outfile::AbstractString,
                           λ::AbstractVector,
                           k::AbstractVector,
                           growth::AbstractArray,
                           speed::AbstractArray)

    h5open(outfile, "w") do f
        write(f, "λ", λ)
        write(f, "k", k)
        write(f, "growth_rate", growth)
        write(f, "phase_speed", speed)

        attributes(f["λ"])["standard_name"] = "wavelength"
        attributes(f["λ"])["units"] = "km"

        attributes(f["k"])["standard_name"] = "non-dimensional wavenumber"
        attributes(f["k"])["units"] = "km/km"

        attributes(f["growth_rate"])["standard_name"] = "modal growth rate"
        attributes(f["growth_rate"])["units"] = "1/day"

        attributes(f["phase_speed"])["standard_name"] = "phase speed"
        attributes(f["phase_speed"])["units"] = "m/s"
    end

    return nothing
end

end # module Diagnostics
