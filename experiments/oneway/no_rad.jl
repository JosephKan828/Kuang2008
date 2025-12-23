using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using LinearAlgebra
using Base.Threads
using HDF5
using Kuang2008
using Measures
using Plots

# ============================================
# Setup
# ============================================
BLAS.set_num_threads(12)

const ROOT     = "/home/b11209013/2025_Research/Kuang2008/"
const WORKPATH = "/work/b11209013/2025_Research/Kuang2008/"
const DATAPATH = ROOT * "data/"
const OUTPATH  = WORKPATH * "output/Oneway/"
const FIGPATH  = ROOT * "figures/Oneway/"

mkpath(OUTPATH)
mkpath(FIGPATH)

# ============================================
# Load data
# ============================================
ρ0, p0, T0, z_bg = load_background(DATAPATH * "background.h5")
G1, G2           = load_vertical_modes(DATAPATH * "vertical_mode.h5")
x, z, t          = load_domain(DATAPATH * "domain.h5")
k                = load_wavenumbers(DATAPATH * "inv_mat.h5")

Nt = length(t)
Nk = length(k)
Nv = 4  # T1, T2, q, J_1

params = default_params("conv_only")

# ============================================
# Initial condition & state container
# ============================================
init = randn(ComplexF64, Nv, Nk)*0.1 .+
       1im .* randn(ComplexF64, Nv, Nk)*0.1

state_vec = zeros(ComplexF64, Nt, Nv, Nk)

# ============================================
# 1. Run full model
# ============================================
run_full_model!(state_vec, t, k, init; params=params, mode=:oneway)

save_state(OUTPATH * "state.h5", state_vec, t, k,
           ["T1","T2","q","J1"])

# ============================================
# 2. Diagnostics: growth rate & phase speed
#    (eigen-analysis of linear operator)
# ============================================
function compute_growth_phase(k::AbstractVector{<:Real},
                              params; Nv::Int = 6)

    Nk = length(k)
    growth = zeros(Float64, Nv, Nk)
    speed  = zeros(Float64, Nv, Nk)

    @threads for j in eachindex(k)
        eigval = eigvals(coeff_matrix(k[j]; param=params))

        σ = real.(eigval)
        c = -imag.(eigval) ./ k[j] * (4320000.0 / 86400.0)  # m/s, same as before

        growth[:, j] .= σ[1:Nv]
        speed[:,  j] .= c[1:Nv]
    end

    return growth, speed
end

growth, speed = compute_growth_phase(k, params; Nv=Nv)

# reconstruct λ consistent with your old script: λ = 2π * 4320 / k
λ = @. 2.0 * π * 4320.0 / k

# Save diagnostics to HDF5
h5open(OUTPATH * "diagnose.h5", "w") do f
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

# ============================================
# 3. Plot diagnostics: most unstable mode
# ============================================

# mimic Python: x = 40000 / λ
x = 40000.0 ./ λ     # "non-dimensional wavelength" style

# find most unstable mode (max growth over modes, per k)
idx_max = [argmax(view(growth, :, j)) for j in 1:Nk]

growth_max = [growth[idx_max[j], j] for j in 1:Nk]
speed_max  = [speed[idx_max[j],  j] for j in 1:Nk]

# --- plotting style (rough rcParams analogue) ---
default(
    dpi            = 500,
    guidefontsize  = 14,
    titlefontsize  = 16,
    tickfontsize   = 12,
    legendfontsize = 11,
    framestyle     = :box,
    grid           = :both,
    gridalpha      = 0.25,
    gridlinewidth  = 0.5,
    gridstyle      = :dash,
)

# 3.1 Growth rate of most unstable mode
p1 = plot(
    x, growth_max;
    lw     = 2.5,
    color  = :black,
    label  = "Max growth over modes",
    xlabel = "",
    ylabel = "Growth Rate (1/day)",
    title  = "Growth Rate of the Most Unstable Mode",
    xlim   = (0, 30),
    # ylim   = (0, 0.13),
    xticks = 0:5:30,
    yticks = 0:0.02:0.12,
)

# 3.2 Phase speed of all modes + highlight most unstable
p2 = plot(
    xlabel = "Non-dimensional Wavelength (40000 km / λ)",
    ylabel = "Phase Speed (m s^{-1})",
    title  = "Phase Speed of All Modes",
    xlim   = (0, 30),
    ylim   = (-5, 60),
    xticks = 0:5:30,
    yticks = 0:10:60,
)

# scatter all modes
for i in 1:Nv
    scatter!(
        p2,
        x, speed[i, :];
        markersize       = 3,
        markercolor       = :blue,
        markerstrokecolor = :black,
        alpha            = 0.3,
        markerstrokewidth = 0,
        label            = (i == 1 ? "All modes" : ""),
    )
end

# highlight most unstable
scatter!(
    p2,
    x, speed_max;
    markersize        = 6,
    markercolor       = :white,
    markerstrokecolor = :black,
    markerstrokewidth = 1.2,
    label             = "Most-unstable (max growth)",
)

plt_diag = plot(p1, p2; layout=(2, 1), size=(900, 1100), margin=10mm)
savefig(plt_diag, FIGPATH * "diagnostic.png")