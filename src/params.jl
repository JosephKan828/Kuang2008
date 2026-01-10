module Params

using HDF5

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
function default_params(exptype::Union{String,Symbol} = "no_rad",
                        rad_scaling::Float64 = 0.0)

    qt_file :: String = "/work/b11209013/2025_Research/MSI/Rad_Stuff/qt_coeff.h5"
    w_file  :: String = "/work/b11209013/2025_Research/MSI/Rad_Stuff/w_coeff.h5"

    RT1_LW, RT1_SW, RT2_LW, RT2_SW, Rq_LW, Rq_SW = h5open(qt_file, "r") do qt_f
        return (
        read(qt_f, "RT1_lw"),
        read(qt_f, "RT1_sw"),
        read(qt_f, "RT2_lw"),
        read(qt_f, "RT2_sw"),
        read(qt_f, "Rq_lw"),
        read(qt_f, "Rq_sw")
       )
    end
    
    Rw1_LW, Rw1_SW, Rw2_LW, Rw2_SW = h5open(w_file, "r") do w_f
        return (
        read(w_f, "Rw1_lw"),
        read(w_f, "Rw1_sw"),
        read(w_f, "Rw2_lw"),
        read(w_f, "Rw2_sw")
       )
    end
 

    exptype_sym = exptype isa String ? Symbol(exptype) : exptype

    if exptype_sym === :no_rad
        # --- Convection only ---
        return ModelParams(
            1.4, 0.0, 1.0, 2.0, 1.0, 0.5, 1.1, -1.0, 0.3, 1.0,
            1.0, 0.7, 4.0, 0.5, 1/12, 0.1,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0
        )

    elseif exptype_sym === :qt_cld_rad
        # --- Convection + (moisture + temperature) radiation ---
        return ModelParams(
            1.4, 0.0, 1.0, 2.0, 1.0, 0.5, 1.1, -1.0, 0.3, 1.0,
            1.0, 0.7, 4.0, 0.5, 1/12, 0.1,
            RT1_LW[1] * rad_scaling, RT1_SW[1] * rad_scaling,
            RT1_LW[2] * rad_scaling, RT1_SW[2] * rad_scaling,
            RT2_LW[1] * rad_scaling, RT2_SW[1] * rad_scaling,
            RT2_LW[2] * rad_scaling, RT2_SW[2] * rad_scaling,
            Rq_LW[1]  * rad_scaling, Rq_SW[1]  * rad_scaling,
            Rq_LW[2]  * rad_scaling, Rq_SW[2]  * rad_scaling,
            Rw1_LW[1] * rad_scaling, Rw1_SW[1] * rad_scaling,
            Rw1_LW[2] * rad_scaling, Rw1_SW[2] * rad_scaling,
            Rw2_LW[1] * rad_scaling, Rw2_SW[1] * rad_scaling,
            Rw2_LW[2] * rad_scaling, Rw2_SW[2] * rad_scaling
        )

    elseif exptype_sym === :qt_rad
        # --- Convection + (moisture + temperature) radiation ---
        return ModelParams(
            1.4, 0.0, 1.0, 2.0, 1.0, 0.5, 1.1, -1.0, 0.3, 1.0,
            1.0, 0.7, 4.0, 0.5, 1/12, 0.1,
            RT1_LW[1] * rad_scaling, RT1_SW[1] * rad_scaling,
            RT1_LW[2] * rad_scaling, RT1_SW[2] * rad_scaling,
            RT2_LW[1] * rad_scaling, RT2_SW[1] * rad_scaling,
            RT2_LW[2] * rad_scaling, RT2_SW[2] * rad_scaling,
            Rq_LW[1]  * rad_scaling, Rq_SW[1]  * rad_scaling,
            Rq_LW[2]  * rad_scaling, Rq_SW[2]  * rad_scaling,
            0.0 * rad_scaling      , 0.0 * rad_scaling,
            0.0 * rad_scaling      , 0.0 * rad_scaling,
            0.0 * rad_scaling      , 0.0 * rad_scaling,
            0.0 * rad_scaling      , 0.0 * rad_scaling
        )
    
    elseif exptype_sym === :cld_rad
        # --- Convection + (moisture + temperature) radiation ---
        return ModelParams(
            1.4, 0.0, 1.0, 2.0, 1.0, 0.5, 1.1, -1.0, 0.3, 1.0,
            1.0, 0.7, 4.0, 0.5, 1/12, 0.1,
            0.0 * rad_scaling      , 0.0 * rad_scaling,
            0.0 * rad_scaling      , 0.0 * rad_scaling,
            0.0 * rad_scaling      , 0.0 * rad_scaling,
            0.0 * rad_scaling      , 0.0 * rad_scaling,
            0.0 * rad_scaling      , 0.0 * rad_scaling,
            0.0 * rad_scaling      , 0.0 * rad_scaling,
            Rw1_LW[1] * rad_scaling, Rw1_SW[1] * rad_scaling,
            Rw1_LW[2] * rad_scaling, Rw1_SW[2] * rad_scaling,
            Rw2_LW[1] * rad_scaling, Rw2_SW[1] * rad_scaling,
            Rw2_LW[2] * rad_scaling, Rw2_SW[2] * rad_scaling
        )


    else
        error("Invalid experiment type: $exptype. " *
                "Choose one of :conv_only, :conv_radiation_full, :conv_radiation_moist.")
    end
end

end # module Params
