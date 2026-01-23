# ==========================================================
# This program is to generate background field of simulation
# ==========================================================

# --------------------------------------------------
# Import necessary packages
# --------------------------------------------------

using Pkg
Pkg.activate(joinpath(@__DIR__))

using HDF5
using LinearAlgebra

# --------------------------------------------------
# Set constants
# --------------------------------------------------

const Γ  :: Float64 = 6.5/1000.0 # Background lapse rate of temperature (K/m)
const g  :: Float64 = 9.81       # Gravitational acceleration (m/s²)
const R  :: Float64 = 287.5      # Gas constant for dry air (J/(kg·K))
const cp :: Float64 = 1004.5     # Specific heat at constant pressure for dry air (J/(kg·K))

const DATA :: String = "/home/b11209013/Kuang2008/data/" # Directory to save data files

const x_scale :: Float64 = 4320e3  # Scaling factor for x-coordinate (km)
# --------------------------------------------------
# Customize domain settings
# --------------------------------------------------

# define grid parameters
## horizontal direction
Δx :: Float64 = 1e3 # Grid spacing in x-direction (m)
Lx :: Float64 = 4e7 # Half domain length in x-direction (m)

## vertical direction
Δz :: Float64 = 1e2  # Grid spacing in z-direction (m)
Lz :: Float64 = 14e3 # Domain height in z-direction (m)

# Time steps
Δt :: Float64 = 0.1 # Time step (day)
Lt :: Float64 = 60  # Total simulation time (day)

# Thermodynamics conditions
T0 :: Float64 = 300.0 # Reference temperature at surface (K)
p0 :: Float64 = 1.0e5 # Reference pressure at surface (Pa)

# wavenumber sampling
λ_min :: Float64 = 540.0e3   # Minimum wavelength (m)
λ_max :: Float64 = 43200.0e3 # Maximum wavelength (m)
dλ    :: Float64 = 540.0e3   # Wavelength interval (m)

# --------------------------------------------------
# Main function
# --------------------------------------------------

function main(
    Δx :: Float64, Lx :: Float64,
    Δz :: Float64, Lz :: Float64,
    Δt :: Float64, Lt :: Float64,
    T0 :: Float64, p0 :: Float64,
    λ_min :: Float64, λ_max :: Float64, dλ :: Float64,
)
    # Design domain
    x :: Vector{Float64} = -Lx:Δx:Lx # x-coordinates (m)
    z :: Vector{Float64} = 0.0:Δz:Lz # z-coordinates (m)
    t :: Vector{Float64} = 0.0:Δt:Lt # time steps (day)

    Nx :: Int64 = length(x) # Number of grid points in x-direction
    Nz :: Int64 = length(z) # Number of grid points in z-direction
    Nt :: Int64 = length(t) # Number of time steps
    
    # Design vertical modes
    G1 :: Matrix{Float64} = reshape(π/2.0*sin.(π*z./(Lz)), Nz, 1) # First vertical mode
    G2 :: Matrix{Float64} = reshape(π/2.0*sin.(2.0*π*z./(Lz)), Nz, 1)   # Second vertical mode

    # Design background thermodynamic profiles
    Tbar :: Vector{Float64} = @. T0 - Γ*z                  # Background temperature (K)
    pbar :: Vector{Float64} = @. p0*(1.0-Γ*z/T0)^(g/(R*Γ)) # Background pressure (Pa)
    ρbar :: Vector{Float64} = @. pbar/(R*Tbar)             # Background density (kg/m³)

    # Design sampling in wavenumber space
    λ :: Vector{Float64} = collect(λ_min:dλ:λ_max) # wavelength of waves

    wnum_cal :: Vector{Float64} = @. 2.0*π*x_scale / λ # wavenumber for calculations
    wnum_dis :: Vector{Float64} = @. Lx / λ            # wavenumber for display in figure

    Nk :: Int64 = length(wnum_cal) # Number of wavenumber samples

    # Design inverse matrix for Fourier transform
    F_inv :: Matrix{ComplexF64} = zeros(ComplexF64, Nk, Nx)

    @inbounds for j in 1:Nk
        for i in 1:Nx
            F_inv[j, i] = exp(-1im * wnum_cal[j] * x[i]/x_scale)
        end
    end

    # Save files
    ## Domain
    h5open(DATA*"domain.h5", "w") do file
        write(file, "x", x)
        write(file, "z", z)
        write(file, "t", t)

        attributes(file)["Description"] = "Domain settings for simulations"

        attributes(file["x"])["Units"] = "m"
        attributes(file["x"])["standard_name"] = "x-coordinate"

        attributes(file["z"])["Units"] = "m"
        attributes(file["z"])["standard_name"] = "z-coordinate"

        attributes(file["t"])["Units"] = "day"
        attributes(file["t"])["standard_name"] = "time"
    end

    ## Vertical modes
    h5open(DATA*"vertical_mode.h5", "w") do file
        write(file, "G1", G1)
        write(file, "G2", G2)

        attributes(file)["Description"] = "Vertical modes for simulations"
    end

    ## Background thermodynamic profiles
    h5open(DATA*"background.h5", "w") do file
        write(file, "Tbar", Tbar)
        write(file, "pbar", pbar)
        write(file, "rho_bar", ρbar)

        attributes(file)["Description"] = "Background thermodynamic profiles for simulations"

        attributes(file["Tbar"])["Units"] = "K"
        attributes(file["Tbar"])["standard_name"] = "temperature"

        attributes(file["pbar"])["Units"] = "Pa"
        attributes(file["pbar"])["standard_name"] = "pressure"

        attributes(file["rho_bar"])["Units"] = "kg/m^3"
        attributes(file["rho_bar"])["standard_name"] = "density"
    end

    # Wavenumber samples and inverse Fourier matrix
    h5open(DATA*"inv_mat.h5", "w") do file
        write(file, "lambda", λ)
        write(file, "wnum_cal", wnum_cal)
        write(file, "wnum_dis", wnum_dis)
        write(file, "F_inv", F_inv)

        attributes(file)["Description"] = "Wavenumber samples and inverse Fourier matrix for simulations"

        attributes(file["lambda"])["Units"] = "m"
        attributes(file["lambda"])["standard_name"] = "wavelength"

        attributes(file["wnum_cal"])["Units"] = "None"
        attributes(file["wnum_cal"])["standard_name"] = "wavenumber"

        attributes(file["wnum_dis"])["Units"] = "None"
        attributes(file["wnum_dis"])["standard_name"] = "wavenumber"
    end
end

# Execute main function
main(Δx, Lx, Δz, Lz, Δt, Lt, T0, p0, λ_min, λ_max, dλ)