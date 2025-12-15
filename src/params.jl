module Params

export ModelParams, default_params

# ============================================
# 1. Model parameter structure
# ============================================

"""
    ModelParams

Holds physical and empirical parameters for the linear dynamical model.

Fields (all `Float64`):

- `a1, a2`   : relations between vertical motion and moisture
- `b1, b2`   : relations between hb and temperature
- `c1, c2`   : phase speeds (e.g. for two vertical modes)
- `d1, d2`   : convection / moisture feedback coefficients
- `m1, m2`   : moisture–heating coupling parameters
- `r0`       : relation between L and U
- `rq`       : relation between U and q+
- `F`        : C.-C. scaling on MSE
- `f`        : weighting between mode-1 and mode-2 temperature
- `τL`       : relaxation timescale on convective heating
- `ϵ`        : linear damping rate
- `RT11,...` : radiative feedbacks on temperature
- `Rq1,Rq2`  : radiative feedbacks on moisture
"""
struct ModelParams
    a1      :: Float64
    a2      :: Float64
    b1      :: Float64
    b2      :: Float64
    c1      :: Float64
    c2      :: Float64
    d1      :: Float64
    d2      :: Float64
    m1      :: Float64
    m2      :: Float64
    r0      :: Float64
    rq      :: Float64
    F       :: Float64
    f       :: Float64
    τL      :: Float64
    ϵ       :: Float64
    RT11_LW :: Float64
    RT11_SW :: Float64
    RT12_LW :: Float64
    RT12_SW :: Float64
    RT21_LW :: Float64
    RT21_SW :: Float64
    RT22_LW :: Float64
    RT22_SW :: Float64
    Rq1_LW  :: Float64
    Rq1_SW  :: Float64
    Rq2_LW  :: Float64
    Rq2_SW  :: Float64
    Rw11_LW :: Float64
    Rw11_SW :: Float64
    Rw12_LW :: Float64
    Rw12_SW :: Float64
    Rw21_LW :: Float64
    Rw21_SW :: Float64
    Rw22_LW :: Float64
    Rw22_SW :: Float64
end


# ============================================
# 2. Default parameter sets
# ============================================

"""
    default_params(exptype="conv_only", rad_scaling=0.0) -> ModelParams

Return a `ModelParams` instance for a given experiment type.

Arguments
---------
- `exptype`:
    - `"conv_only"` or `:conv_only` — Convection only  
    - `"conv_radiation_full"` or `:conv_radiation_full`
        Convection + (moisture + temperature) radiation  
    - `"conv_radiation_moist"` or `:conv_radiation_moist`
        Convection + (moisture-only) radiation

- `rad_scaling`:
    Scalar factor applied to all radiative feedback coefficients.
"""
function default_params(exptype::Union{String,Symbol} = "conv_only",
                        rad_scaling::Float64 = 0.0)

    exptype_sym = exptype isa String ? Symbol(exptype) : exptype

    if exptype_sym === :conv_only
        # --- Convection only ---
        return ModelParams(
            1.4, 0.0, 1.0, 2.0, 1.0, 0.5, 1.1, -1.0, 0.3, 1.0,
            1.0, 0.7, 4.0, 0.5, 1/12, 0.1,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0
        )

    elseif exptype_sym === :conv_radiation_full
        # --- Convection + (moisture + temperature) radiation ---
        return ModelParams(
            1.4, 0.0, 1.0, 2.0, 1.0, 0.5, 1.1, -1.0, 0.3, 1.0,
            1.0, 0.7, 4.0, 0.5, 1/12, 0.1,
            -0.021   * rad_scaling, 6.04e-5   * rad_scaling,
            -0.0032  * rad_scaling, -0.00017  * rad_scaling,
            -0.0046  * rad_scaling, 8.78e-5   * rad_scaling,
            -0.038   * rad_scaling, 5.98e-5   * rad_scaling,
            4.18     * rad_scaling, 2.52      * rad_scaling,
            11.08    * rad_scaling, -5.32     * rad_scaling,
            3416.61  * rad_scaling, 2269.21   * rad_scaling,
            7331.88  * rad_scaling, -3748.63  * rad_scaling,
            -1395.51 * rad_scaling, -162.58   * rad_scaling,
            -2577.28 * rad_scaling, 1627.12   * rad_scaling
        )

    # elseif exptype_sym === :conv_radiation_moist
    #     # --- Convection + (moisture-only) radiation ---
    #     return ModelParams(
    #         1.4, 0.0, 1.0, 2.0, 1.0, 0.5, 1.1, -1.0, 0.3, 1.0,
    #         1.0, 0.7, 4.0, 0.5, 1/12, 0.1,
    #         0.0, 0.0, 0.0, 0.0,
    #         5.61 * rad_scaling, 3.36 * rad_scaling
    #     )

    else
        error("Invalid experiment type: $exptype. " *
                "Choose one of :conv_only, :conv_radiation_full, :conv_radiation_moist.")
    end
end

end # module Params
