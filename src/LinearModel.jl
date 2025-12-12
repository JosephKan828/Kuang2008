module LinearModel

using LinearAlgebra

include("params.jl")
using .Params: ModelParams, default_params

export ModelParams, default_params, coeff_matrix

# ============================================
# Coefficient matrix generator
# ============================================

"""
    coeff_matrix(kn; param = default_params(), mode = :full)

Compute the complex coefficient matrix for a given wavenumber `kn`.

Keyword arguments
-----------------
- `param` : `ModelParams` struct (default = `default_params()`)
- `mode`  : `:full` (6×6 system) or `:oneway` (4×4 reduced system)

Returns
-------
- `Matrix{ComplexF64}` of size (6,6) or (4,4)
"""
function coeff_matrix(kn::Float64;
                      param::ModelParams = default_params(),
                      mode::Symbol = :full)

    if mode === :full
        # ---------- FULL 6×6 SYSTEM ----------

        # Precompute auxiliary constants
        A = 1.0 - 2.0*param.f + (param.b2 - param.b1)/param.F
        B = 1.0 + (param.b2 + param.b1)/param.F - A*param.r0

        α = -param.d1*(-1.5*param.rq + param.RT11) - param.d2*( 1.5*param.rq + param.RT12)
        β = -param.d1*param.RT21                      - param.d2*param.RT22
        γ = -param.d1*( param.rq + param.Rq1)        - param.d2*(-param.rq + param.Rq2)
        δ = -param.d1*( 1.0 + param.r0)              - param.d2*( 1.0 - param.r0)

        mat = ComplexF64[
            -param.ϵ                  0.0                     (kn*param.c1)^2              0.0                          0.0                      0.0;
             0.0                   -param.ϵ                  0.0                     (kn*param.c2)^2                  0.0                      0.0;
            -1.0                     0.0                   -1.5*param.rq + param.RT11   param.RT21                    param.rq + param.Rq1      1.0 + param.r0;
             0.0                    -1.0                    1.5*param.rq + param.RT12   param.RT22                   -param.rq + param.Rq2      1.0 - param.r0;
             param.a1               param.a2                α                           β                            γ                          δ;
             param.f /B/param.τL   (1-param.f)/B/param.τL  -1.5*A*param.rq/B/param.τL   0.0                          A*param.rq/B/param.τL    -1.0/param.τL
        ]

        return mat

    elseif mode === :oneway
        # ---------- 4×4 REDUCED “ONE-WAY” SYSTEM ----------

        γ0 = (1.0 - param.r0) / (1.0 + param.r0)
        γq = (2.0 * param.rq) / (1.0 + param.r0)

        # denominator D = b1 + b2 * γ0
        D = param.b1 + param.b2 * γ0

        # helper coefficients from d/dt (f T1 + (1-f) T2)
        α0 = param.f * (-1im * kn * param.c1 - param.ϵ) +
             1.5 * γq * (1.0 - param.f)

        β0 = (1.0 - param.f) * (-1im * kn * param.c2 - param.ϵ)

        δ0 = -(1.0 - param.f) * γq

        # coefficients in J1_eq = coeff_T1*T1 + coeff_T2*T2 + coeff_q*q + coeff_J1*J1
        coeff_T1 = (-1.5 * param.b2 * γq - param.F * α0) / D
        coeff_T2 = (-param.F * β0) / D
        coeff_q  = ( param.b2 * γq - param.F * δ0) / D
        coeff_J1 = (-param.F * param.f) / D

        # bottom row of the matrix: dJ1/dt = α*T1 + β*T2 + γ*q + δ*J1
        α = coeff_T1 / param.τL
        β = coeff_T2 / param.τL
        γ = coeff_q  / param.τL
        δ = (coeff_J1 - 1.0) / param.τL

        mat = ComplexF64[
            -1im * kn * param.c1 - param.ϵ      0.0                             0.0            1.0;
             1.5 * γq                           -1im * kn * param.c2 - param.ϵ  -γq           0.0;
             1.5 * param.m2 * γq                0.0                             -param.m2*γq   param.m1;
             α                                  β                               γ              δ
        ]

        return mat

    else
        error("Unknown mode = $mode. Use :full or :oneway.")
    end
end

end # module LinearModel
