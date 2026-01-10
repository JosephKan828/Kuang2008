using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using LinearAlgebra
using Base.Threads
using HDF5
using Kuang2008
using Measures
using Plots

import TOML

function params_from_case(case::AbstractString, rad_scaling::Float64)
    """
    Load case name and find corresponding parameters with params.jl
    """
    return default_params(case, rad_scaling)
end

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

"""
Get a value from TOML Dict with nested keys, else default.
Example: get_cfg(cfg, ("compute","blas_threads"), 12)
"""
function get_cfg(cfg, keys::Tuple, default)
    cfg === nothing && return default
    d = cfg
    for k in keys[1:end-1]
        if !(haskey(d, k))
            return default
        end
        d = d[k]
    end
    lastk = keys[end]
    return haskey(d, lastk) ? d[lastk] : default
end

# ============================================
# Main
# ============================================

function main()
    # -------------------------------------------
    # Get case information
    # -------------------------------------------
    case = length(ARGS) >= 1 ? ARGS[1] : "no_rad" # refer to specific case
    rad_scaling = length(ARGS) >= 2 ? ARGS[2] : "0.001"
    println( ARGS )

    println( "rad_scaling = ", rad_scaling )
    root = abspath(joinpath(@__DIR__, ".."))      # root directory

    # -------------------------------------------
    # Configuration
    # -------------------------------------------
    cfg = load_case_config(root, case)

    if cfg === nothing
        println("No configuration for case = ", case)
        return
    end

    # Specify directories
    ROOT     = "/home/b11209013/2025_Research/Kuang2008/"
    WORKPATH = "/work/b11209013/2025_Research/Kuang2008/"
    DATAPATH = joinpath(ROOT, get_cfg(cfg, ("paths","data_dir"), "data"))
    OUTPATH  = joinpath(WORKPATH, get_cfg(cfg, ("paths","output_dir"), "output"))
    FIGPATH  = joinpath(ROOT, get_cfg(cfg, ("paths","figures_dir"), "figures"))

    # -------------------------------------------
    # Compute setting
    # -------------------------------------------
    blas_threads = Int(get_cfg(cfg, ("compute","blas_threads"), 12))
    BLAS.set_num_threads(blas_threads)
    
    # -------------------------------------------
    # Load data
    # -------------------------------------------
    background_file     = joinpath(DATAPATH, get_cfg(cfg, ("data_files","background"), "background.h5"))
    vertical_mode_file  = joinpath(DATAPATH, get_cfg(cfg, ("data_files","vertical_mode"), "vertical_mode.h5"))
    domain_file         = joinpath(DATAPATH, get_cfg(cfg, ("data_files","domain"), "domain.h5"))
    invmat_file         = joinpath(DATAPATH, get_cfg(cfg, ("data_files","inv_mat"), "inv_mat.h5"))

    ρ0, p0, T0, z_bg = load_background(background_file)
    G1, G2           = load_vertical_modes(vertical_mode_file)
    x, z, t          = load_domain(domain_file)
    k                = load_wavenumbers(invmat_file)

    Nt = length(t)
    Nk = length(k)
    Nv = Int(get_cfg(cfg, ("model","nv"), 6))
    state_names = get_cfg(cfg, ("model","state_names"), ["w1","w2","T1","T2","q","L"])

    println("Finish loading data.")

    # -------------------------------------------
    # Parameters
    # -------------------------------------------
    rad_scaling_float = parse(Float64, rad_scaling)
    params = params_from_case(case, rad_scaling_float)

    # -------------------------------------------
    # Initial conditions
    # -------------------------------------------
    init = Array{ComplexF64}(undef, Nv, Nk)

    scales = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

    for i in 1:Nv
        init[i, :] = randn(ComplexF64, Nk).*(1 + 1im) * scales[i]
    end

    state_vec = zeros(ComplexF64, Nt, Nv, Nk)

    println("Finish initial conditions.")

    # -------------------------------------------
    # Run full model
    # -------------------------------------------
    integration!(state_vec, t, k, init; params=params, mode=:full)

    println("Finish running full model.")


    if case == "no_rad"
        outdir = joinpath(ROOT, "output", case); mkpath(outdir)
    else
        outdir = joinpath(ROOT, "output", case, "rad_scaling=$(rad_scaling)"); mkpath(outdir)
    end

    println("Saving data to ", outdir)

    save_state(joinpath(outdir, "state.h5"), state_vec, t, k,
            ["w1","w2","T1","T2","q","L"])

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
