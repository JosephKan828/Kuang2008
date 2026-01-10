# =============================================
# Post processing for Kuang2008 Model
# =============================================

import sys
import h5py
import numpy as np
import argparse
from pathlib import Path

def is_venv():
    return sys.prefix != sys.base_prefix

if is_venv():
    print( f"Running in venv: {sys.prefix}" )
else:
    print( "Running in GLOBAL environment." )

sys.path.append(
    "/home/b11209013/2025_Research/Kuang2008/post/lib"
)

import Diagnostics #type: ignore
import Reconstruct #type: ignore
import Plot        #type: ignore


def parse_args():
    parser = argparse.ArgumentParser(
        description="Post-processing pipeline for Kuang2008 simulations"
    )

    parser.add_argument("--case", required=True,
                        help="Case name (e.g., no_rad, qt_rad, cld_rad)")
    parser.add_argument("--in-dir", required=True,
                        help="Input directory containing simulation output")
    parser.add_argument("--fig-dir", required=True,
                        help="Directory to save figures")
    parser.add_argument("--post-dir", required=True,
                        help="Directory to save derived diagnostics")

    return parser.parse_args()

# =============================================
# Main
# =============================================

def main() -> None:

    # --------------------------------------------
    # Get case information
    # -------------------------------------------

    args = parse_args()

    case     = args.case
    in_dir   = Path(args.in_dir)
    fig_dir  = Path(args.fig_dir)
    post_dir = Path(args.post_dir)

    # --------------------------------------------
    # Create directories
    # --------------------------------------------

    # Create figure and post directories
    fig_dir.mkdir(parents=True, exist_ok=True)
    post_dir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------
    # Load data
    # --------------------------------------------

    # Load state matrix
    with h5py.File(
        in_dir / "state.h5", "r"
    ) as f:
        state: np.ndarray =  f.get( "state" )[ :, :, :300 ] # Shape: ( Nk, Nv, Nt )
        time : np.ndarray =  f.get( "time" )[ :300 ]
        var  : np.ndarray =  f.get( "variables" )[ ... ]
        wnum : np.ndarray =  f.get( "wavenumber" )[ ... ]

    state = state.astype( np.complex64, copy=False )
    state: np.ndarray = state.transpose( 1, 0, 2 )

    # Load linear operators
    with h5py.File(
        in_dir / "optrs.h5", "r"
    ) as f:
        optrs: np.ndarray = f.get( "operators" )[ ... ] # shape: ( Nv, Nv, Nk )

    optrs: np.ndarray = optrs.astype( np.complex64, copy=False )

    # Load inverse matrix
    with h5py.File(
    "/home/b11209013/2025_Research/Kuang2008/data/inv_mat.h5",
        "r"
    ) as f:
        x      : np.ndarray = f.get( "x" )[ ... ] # Shape: ( Nx, )
        inv_mat: np.ndarray = f.get( "inverse matrix" )[ ... ] # shape: ( Nx, Nk )

    inv_mat: np.ndarray = inv_mat.astype( np.complex64, copy=False )

    print( "Finish Loading Data" )

    # ---------------------------------------------
    # Calculate Modal Growth and Phase Speed
    # ---------------------------------------------

    # Calculate growth rate
    σ: np.ndarray = Diagnostics.modal_growth_rate(
        optrs, wnum
    ) # Shape: ( Nv, Nk )

    # Calculate phase speed
    c: np.ndarray = Diagnostics.phase_speed(
        optrs, wnum
    ) # Shape: ( Nv, Nk )

    print( "Finish Calculating Diagnostics" )

    # ------------------------------------------------
    # Plot Diagnostics
    # ------------------------------------------------

    # Plot.plot_growth_rate( σ   , wnum, fig_dir / "growth_rate.png" )
    Plot.plot_diagnostics( σ, c, wnum, fig_dir / "growth_rate.png", fig_dir / "phase_speed.png" )

    print( "Finish Plotting Diagnostics" )

    # -------------------------------------------------
    # Reconstruct Fields
    # -------------------------------------------------

    # Calculate J1 and J2
    L: np.ndarray = state[ -1, :, : ]
    U: np.ndarray = L + 0.7*( state[ -2, :, : ]-1.5*state[ 2, :, : ] )

    J1: np.ndarray = ( L + U )[ None, :, : ]
    J2: np.ndarray = ( L - U )[ None, :, : ]

    state_J: np.ndarray = np.concatenate( [ state, J1, J2 ], axis=0 ).astype( np.complex64, copy=False )

    # design basis
    levels   : np.ndarray = np.linspace( 0.0, 14000.0, 141, dtype=np.float32 )
    level_max: np.float32 = levels.max()
    G1    : np.ndarray = np.pi/2 * np.sin( np.pi*levels/level_max )
    G2    : np.ndarray = np.pi/2 * np.sin( 2*np.pi*levels/level_max )

    basis_G : np.ndarray = np.stack( [ G1, G2 ], axis=0 )
    basis_T : np.ndarray = basis_G * 0.0033

    basis_v : np.ndarray = np.concatenate( [ basis_G, basis_T, basis_T ], axis=0 )

    # Choose wavelength and wavenumber
    target_λ: np.float64 = np.float64( 8640.0 )
    target_k: np.float64 = 2 * np.pi * 4320.0/target_λ

    kidx    : np.int64   = np.argmin( np.abs( target_k - wnum ) )

    # constrain x domain of inverse matrix
    x_mask : np.ndarray = np.logical_and( x>=-4320000, x<=4320000 )

    inv_mat_lim : np.ndarray = inv_mat[ x_mask, : ].astype( np.complex64, copy=False )
    x_lim       : np.ndarray = x[ x_mask ]

    # Reconstruct
    FourierBasis : np.ndarray = inv_mat_lim[ :, kidx ][ None, :, None ]
    pcs          : np.ndarray = state[ :, kidx, : ][ :, None, : ]

    profile      : np.ndarray = np.einsum( "vz,vxt->vzxt", basis_v, FourierBasis*pcs, optimize=True )

    # w_prof : list[np.ndarray] = Reconstruct.reconstruct(
    #     kidx,
    #     state[ :2, :, :],
    #     np.array( [ G1, G2 ] ),
    #     inv_mat_lim
    # )
    # print( "Finish computing w" )

    # T_prof : list[np.ndarray] = Reconstruct.reconstruct(
    #     kidx,
    #     state[ 2:4, :, : ],
    #     np.array( [ G1, G2 ] )*0.0033,
    #     inv_mat_lim
    # )

    # print( "Finish computing T" )

    # J_prof : list[np.ndarray] = Reconstruct.reconstruct(
    #     kidx,
    #     np.array( [ J1, J2 ] ),
    #     np.array( [ G1, G2 ] ),
    #     inv_mat_lim
    # )
    # print( "Finish computing J" )

    print( "Finish reconstructing profiles" )

    # --------------------------------------------
    # Plot profile
    # --------------------------------------------

    # combining data 
    w: np.ndarray = profile[ 0 ] + profile[ 1 ]
    T: np.ndarray = profile[ 2 ] + profile[ 3 ]
    J: np.ndarray = profile[ 4 ] + profile[ 5 ]

    # Design levels

    Tmax: np.float64 = np.nanmax( np.abs( T ) ) // 0.005 * 0.005
    wmax: np.float64 = np.nanmax( np.abs( w ) ) // 0.01 * 0.01
    Jmax: np.float64 = np.nanmax( np.abs( J ) ) // 1

    print( "Maximum T", np.nanmax( np.abs( T ) ) )

    Tlevel: np.ndarray = np.linspace( -Tmax, Tmax, 21 )
    wlevel: np.ndarray = np.linspace( -wmax, wmax, 21 )
    Jlevel: np.ndarray = np.linspace( -Jmax, Jmax, 21 )

    Tlevel = Tlevel[ np.abs( Tlevel ) >= 1e-5 ]
    wlevel = wlevel[ np.abs( wlevel ) >= 1e-5 ]
    Jlevel = Jlevel[ np.abs( Jlevel ) >= 1e-5 ]

    # Plot animation

    Plot.plot_animation(
        np.linspace( -4320000, 4320000, w.shape[ 1 ] ),
        np.linspace( 0, 14000, w.shape[ 0 ] ),
        J, T, w, "RdBu_r", None, Tlevel, wlevel,
        "Convective-heating Only",
        fig_dir / "profile_evo.mp4")

    print( "Finish plotting animation" )

if __name__ == "__main__":
    main()
