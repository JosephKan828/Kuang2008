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
        # ==================================================
        # FULL 6×6 LINEAR SYSTEM
        # ==================================================
        
        # --------------------------------------------------
        # 1. Derived nondimensional constants
        # --------------------------------------------------
        ϵ  = param.ϵ
        rq = param.rq
        r0 = param.r0
        τL = param.τL
        f  = param.f
        
        c1 = param.c1
        c2 = param.c2

        A = 1.0 - 2.0*f + (param.b2 - param.b1)/param.F
        B = 1.0 + (param.b2 + param.b1)/param.F - A*r0
        
        invBτ = 1.0 / (B * τL)
        
        # --------------------------------------------------
        # 2. Radiative tendency aggregation (LW + SW)
        # --------------------------------------------------
        Rq1  = param.Rq1_LW  + param.Rq1_SW
        Rq2  = param.Rq2_LW  + param.Rq2_SW
        
        RT11 = param.RT11_LW + param.RT11_SW
        RT12 = param.RT12_LW + param.RT12_SW
        RT21 = param.RT21_LW + param.RT21_SW
        RT22 = param.RT22_LW + param.RT22_SW
        
        Rw11 = param.Rw11_LW + param.Rw11_SW
        Rw12 = param.Rw12_LW + param.Rw12_SW
        Rw21 = param.Rw21_LW + param.Rw21_SW
        Rw22 = param.Rw22_LW + param.Rw22_SW
        
        # --------------------------------------------------
        # 3. Moisture–dynamics coupling coefficients
        # --------------------------------------------------
        d1, d2 = param.d1, param.d2
        
        α = param.a1 - d1*Rw11 - d2*Rw12
        β = param.a2 - d1*Rw21 - d2*Rw22
        
        γ = -d1*(RT11 - 1.5*rq) - d2*(RT12 + 1.5*rq)
        δ = -d1*RT21 - d2*RT22
        
        ζ = -d1*(Rq1 + rq) - d2*(Rq2 - rq)
        η = -d1*(1.0 + r0) - d2*(1.0 - r0)
        
        # --------------------------------------------------
        # 4. Linear operator assembly
        # --------------------------------------------------
        mat = ComplexF64[
            # (T₁, T₂)
            -ϵ        0.0       (c1*kn)^2                 0.0        0.0              0.0;
            0.0      -ϵ          0.0        (c2*kn)^2                 0.0              0.0;
        
            # (w₁, w₂)
            Rw11-1.0  Rw21   RT11 - 1.5*rq   RT21   rq + Rq1    1.0 + r0;
            Rw12      Rw22-1.0 RT12 + 1.5*rq RT22  -rq + Rq2   1.0 - r0;
        
            # (q)
            α         β         γ            δ         ζ                 η;
        
            # (L)
            f*invBτ   (1-f)*invBτ  -1.5*A*rq*invBτ   0.0   A*rq*invBτ   -1.0/τL
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
