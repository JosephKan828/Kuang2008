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
    a1   :: Float64
    a2   :: Float64
    b1   :: Float64
    b2   :: Float64
    c1   :: Float64
    c2   :: Float64
    d1   :: Float64
    d2   :: Float64
    m1   :: Float64
    m2   :: Float64
    r0   :: Float64
    rq   :: Float64
    F    :: Float64
    f    :: Float64
    τL   :: Float64
    ϵ    :: Float64
    RT11 :: Float64
    RT12 :: Float64
    RT21 :: Float64
    RT22 :: Float64
    Rq1  :: Float64
    Rq2  :: Float64
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
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        )

    elseif exptype_sym === :conv_radiation_full
        # --- Convection + (moisture + temperature) radiation ---
        return ModelParams(
            1.4, 0.0, 1.0, 2.0, 1.0, 0.5, 1.1, -1.0, 0.3, 1.0,
            1.0, 0.7, 4.0, 0.5, 1/12, 0.1,
            -0.042 * rad_scaling, -0.0087 * rad_scaling, -0.011 * rad_scaling,
            -0.069 * rad_scaling,  5.61  * rad_scaling,  3.36  * rad_scaling
        )

    elseif exptype_sym === :conv_radiation_moist
        # --- Convection + (moisture-only) radiation ---
        return ModelParams(
            1.4, 0.0, 1.0, 2.0, 1.0, 0.5, 1.1, -1.0, 0.3, 1.0,
            1.0, 0.7, 4.0, 0.5, 1/12, 0.1,
            0.0, 0.0, 0.0, 0.0,
            5.61 * rad_scaling, 3.36 * rad_scaling
        )

    else
        error("Invalid experiment type: $exptype. " *
                "Choose one of :conv_only, :conv_radiation_full, :conv_radiation_moist.")
    end
end

end # module Params
