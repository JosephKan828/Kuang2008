using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using LinearAlgebra
using Base.Threads
using HDF5
using Kuang2008
using Measures
using Plots

import TOML

const PROJECT_ROOT = "/home/b11209013/2025_Research/Kuang2008/"
const WORK_ROOT    = "/work/b11209013/Kuang2008/"

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
    cfg = load_case_config(PROJECT_ROOT, case)

    @assert cfg !== nothing "No configuration for case = ", case

    data_cfg = cfg["data_files"]
    path_cfg = cfg["paths"]
    comp_cfg = cfg["compute"]
    mod_cfg  = cfg["model"]    

    DATAPATH = joinpath(PROJECT_ROOT, path_cfg["data_dir"])

    # -------------------------------------------
    # Compute setting
    # -------------------------------------------
    blas_threads = Int(comp_cfg["blas_threads"])
    BLAS.set_num_threads(1)
    
    # -------------------------------------------
    # Load data
    # -------------------------------------------

    ρ0, p0, T0, z_bg = load_background(joinpath(DATAPATH, data_cfg["background"]))
    G1, G2           = load_vertical_modes(joinpath(DATAPATH, data_cfg["vertical_mode"]))
    x, z, t          = load_domain(joinpath(DATAPATH, data_cfg["domain"]))
    k                = load_wavenumbers(joinpath(DATAPATH, data_cfg["inv_mat"]))

    Nt, Nv, Nk = length(t), Int(mod_cfg["nv"]), length(k)

    state_names = mod_cfg["state_names"]

    println("Finish loading data.")

    # -------------------------------------------
    # Initial conditions
    # -------------------------------------------
    init = zeros(ComplexF64, Nv, Nk)

    scales = [10.0, 10.0, 10.0, 10.0, 0.1, 10] # scale of initial field

    for i in 1:Nv
        init[i, :] .= (2.0.*rand(ComplexF64, Nk) .- 1.0) .* (1.0 + 1.0im) .* scales[i]
    end

    state_vec = zeros(ComplexF64, Nt, Nv, Nk)
    params = default_params(case, rad_scaling) # setup parameter set

    # dump( params ) # show parameters in this structure field

    println("Finish initial conditions.")

    # -------------------------------------------
    # Run full model
    # -------------------------------------------
    println("Running model for $case...")

    @time integration!(state_vec, t, k, init; params=params, mode=:full)

    println("Finish running full model.")

    subfolder = case == "no_rad" ? case : joinpath(case, "rad_scaling=$rad_scaling_str")
    outdir = joinpath(WORK_ROOT, "output", subfolder)
    mkpath(outdir)

    println("Saving data to ", outdir)

    save_state(joinpath(outdir, "state.h5"), state_vec, t, k, state_names)

    # -------------------------------------------
    # Save linear operator
    # -------------------------------------------
    
    optrs = Array{ComplexF64}(undef, Nk, Nv, Nv)

    @threads for j in eachindex( k )
        optrs[ j, :, : ] = coeff_matrix( k[ j ]; param=params )
    end
    
    save_optrs( joinpath( outdir, "optrs.h5" ), optrs, k )

end

main()
