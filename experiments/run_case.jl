# ============================================
# Executing simulations
# ============================================

# ============================================
# Load packages
# ============================================
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using LinearAlgebra
using Base.Threads
using HDF5
using Kuang2008
using Measures
using Plots

import TOML

# ============================================
# Helper functions
# ============================================

# functions for loading configuration
function load_case_config(root::AbstractString, case::AbstractString)

    """
    load case configuration
    """

    cfg_path = joinpath(root, "configs", case * ".toml")
    if isfile(cfg_path)
        return TOML.parsefile(cfg_path)
    end
    return nothing
end

# ============================================
# Main
# ============================================

function main()
    # -------------------------------------------
    # Get case information
    # -------------------------------------------
    @assert length(ARGS) <= 3 "Usage: julia run_case.jl [case] [rad_scaling]"

    case            = ARGS[1] # refer to specific case
    rad_scaling_str = ARGS[2]
    rad_scaling     = parse(Float64, rad_scaling_str)  # radiation scaling factor

    root = abspath(joinpath(@__DIR__, ".."))      # root directory

    # -------------------------------------------
    # Configuration
    # -------------------------------------------
    cfg = load_case_config(root, case)

    @assert cfg !== nothing "No configuration for case = ", case

    data_cfg = cfg["data_files"]
    path_cfg = cfg["paths"]
    comp_cfg = cfg["compute"]
    mod_cfg  = cfg["model"]    

    DATAPATH = joinpath(path_cfg["root"], path_cfg["data_dir"])

    # -------------------------------------------
    # Compute setting
    # -------------------------------------------
    blas_threads = Int(comp_cfg["blas_threads"])
    BLAS.set_num_threads(blas_threads)
    
    # -------------------------------------------
    # Load data
    # -------------------------------------------

    ρ0, p0, T0          = load_background(joinpath(DATAPATH, data_cfg["background"]))
    G1, G2              = load_vertical_modes(joinpath(DATAPATH, data_cfg["vertical_mode"]))
    x, z, t             = load_domain(joinpath(DATAPATH, data_cfg["domain"]))
    λ, kcal, kdis, Finv = load_inv_mat(joinpath(DATAPATH, data_cfg["inv_mat"]))

    Nt, Nv, Nk = length(t), Int(mod_cfg["nv"]), length(kcal)

    state_names = mod_cfg["state_names"]

    println("Finish loading data.")

    # -------------------------------------------
    # Initial conditions
    # -------------------------------------------
    init = zeros(ComplexF64, Nv, Nk)

    scales = [10.0, 10.0, 10.0, 10.0, 0.1, 100.0] # scale of initial field

    for i in 1:Nv
        init[i, :] .= rand(ComplexF64, Nk) .* scales[i]
    end

    state_vec = zeros(ComplexF64, Nt, Nv, Nk)
    params = default_params(case, rad_scaling) # setup parameter set

    println("Finish initial conditions.")

    # -------------------------------------------
    # Run full model
    # -------------------------------------------
    println("Running model for $case...")

    @time integration!(state_vec, t, kcal, init; params=params, mode=:full)

    println("Finish running full model.")

    subfolder = case == "no_rad" ? case : joinpath(case, "rad_scaling=$rad_scaling_str")
    outdir = joinpath(path_cfg["work"], "output", subfolder)
    mkpath(outdir)

    println("Saving data to ", outdir)

    save_state(joinpath(outdir, "state.h5"), state_vec, t, kcal, state_names)

    # -------------------------------------------
    # Save linear operator
    # -------------------------------------------
    
    optrs = Array{ComplexF64}(undef, Nk, Nv, Nv)

    @threads for j in eachindex( kcal )
        optrs[ j, :, : ] = coeff_matrix( kcal[ j ]; param=params )
    end
    save_optrs( joinpath( outdir, "optrs.h5" ), optrs, kcal )

end

main()