module Kuang2008

using LinearAlgebra

include("LinearModel.jl")
include("Simulation.jl")
include("IO.jl")
include("Diagnostics.jl")

using .LinearModel
using .Simulation
using .IO
using .Diagnostics

export ModelParams, default_params, coeff_matrix
export integration!, run_full_model!
export load_background, load_vertical_modes, load_domain, load_wavenumbers, save_state

# New exports for diagnostics
export diagnose_growth_phase, save_growth_phase, run_diagnostics

"""
    run_diagnostics(; mode_set="conv_only",
                     λ = collect(540.0:540.0:43200.0),
                     Nv::Int = 6,
                     outfile = "output/Full/Origin/diagnose.h5")

Convenience wrapper that:

1. Builds `params = default_params(mode_set)`
2. Runs `diagnose_growth_phase`
3. Saves to `outfile` via `save_growth_phase`
"""
function run_diagnostics(; mode_set::AbstractString = "conv_only",
                         λ = collect(540.0:540.0:43200.0),
                         Nv::Int = 6,
                         outfile::AbstractString = "output/Full/Origin/diagnose.h5")

    params = default_params(mode_set)
    λ_vec, k_vec, growth, speed = Diagnostics.diagnose_growth_phase(params;
                                                                    λ = λ,
                                                                    Nv = Nv)
    Diagnostics.save_growth_phase(outfile, λ_vec, k_vec, growth, speed)
    return nothing
end

end # module Kuang2008
