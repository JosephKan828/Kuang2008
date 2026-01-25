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
    print(f"Running in venv: {sys.prefix}")
else:
    print("Running in GLOBAL environment.")

sys.path.append("/home/b11209013/Kuang2008/post/lib")

import Diagnostics  # type: ignore
import Reconstruct  # type: ignore
import Plot  # type: ignore


def parse_args():
    parser = argparse.ArgumentParser(
        description="Post-processing pipeline for Kuang2008 simulations"
    )

    parser.add_argument(
        "--case", required=True, help="Case name (e.g., no_rad, qt_rad, cld_rad)"
    )
    parser.add_argument(
        "--in-dir", required=True, help="Input directory containing simulation output"
    )
    parser.add_argument("--fig-dir", required=True, help="Directory to save figures")

    return parser.parse_args()


# =============================================
# Main
# =============================================


def main() -> None:

    # --------------------------------------------
    # Get case information
    # -------------------------------------------

    args = parse_args()

    case = args.case
    in_dir = Path(args.in_dir)
    fig_dir = Path(args.fig_dir)

    # --------------------------------------------
    # Create directories
    # --------------------------------------------

    # Create figure and post directories
    fig_dir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------
    # Load data
    # --------------------------------------------

    # Load state matrix
    with h5py.File(in_dir / "state.h5", "r") as f:
        state: np.ndarray = np.array(f.get("state"))[
            :, :, :300
        ]  # Shape: ( Nk, Nv, Nt )
        time: np.ndarray = np.array(f.get("time"))[:300]
        var: np.ndarray = np.array(f.get("variables"))

    state = state.astype(np.complex64, copy=False)
    state: np.ndarray = state.transpose(1, 0, 2)

    # Load linear operators
    with h5py.File(in_dir / "optrs.h5", "r") as f:
        optrs: np.ndarray = np.array(f.get("operators"))  # shape: ( Nv, Nv, Nk )

    optrs: np.ndarray = optrs.astype(np.complex64, copy=False)

    # Load inverse matrix
    with h5py.File("/home/b11209013/Kuang2008/data/inv_mat.h5", "r") as f:
        λ: np.ndarray = np.array(f.get("lambda"))
        kcal: np.ndarray = np.array(f.get("wnum_cal"))
        kdis: np.ndarray = np.array(f.get("wnum_dis"))
        inv_mat: np.ndarray = np.array(f.get("F_inv"))  # shape: ( Nx, Nk )

    inv_mat: np.ndarray = inv_mat.astype(np.complex64, copy=False)

    # load domain
    with h5py.File("/home/b11209013/Kuang2008/data/domain.h5", "r") as f:
        x: np.ndarray = np.array(f.get("x"))  # shape: ( Nx, )
        z: np.ndarray = np.array(f.get("z"))  # shape: ( Nz, )

    # Load vertical modes
    with h5py.File("/home/b11209013/Kuang2008/data/vertical_mode.h5", "r") as f:
        G1: np.ndarray = np.array(f.get("G1")).squeeze()
        G2: np.ndarray = np.array(f.get("G2")).squeeze()

    print("Finish Loading Data")

    # ---------------------------------------------
    # Calculate Modal Growth and Phase Speed
    # ---------------------------------------------

    # Calculate growth rate
    σ: np.ndarray = Diagnostics.modal_growth_rate(optrs, kcal)  # Shape: ( Nv, Nk )

    # Calculate phase speed
    c: np.ndarray = Diagnostics.phase_speed(optrs, kcal)  # Shape: ( Nv, Nk )

    print("Finish Calculating Diagnostics")

    # ------------------------------------------------
    # Identify maximum growth and phase speed
    # ---------------------------------------------

    max_idx: np.ndarray = np.argmax(σ, axis=0)

    σmax: np.ndarray = np.take_along_axis(σ, max_idx[None, :], axis=0)[0]
    cmax: np.ndarray = np.take_along_axis(c, max_idx[None, :], axis=0)[0]

    # ------------------------------------------------
    # Plot Diagnostics
    # ------------------------------------------------
    print("Figure is saved in: ", fig_dir)
    Plot.plot_diagnostics(
        σmax, cmax, c, kdis, fig_dir / "growth_rate.png", fig_dir / "phase_speed.png"
    )

    print("Finish Plotting Diagnostics")

    # -------------------------------------------------
    # Reconstruct Fields
    # -------------------------------------------------

    # drop moisture
    state = np.array([state[0], state[1], state[2], state[3], state[5], state[6]])

    # design basis
    basis_G: np.ndarray = np.stack([G1, G2], axis=0)
    basis_T: np.ndarray = basis_G * 0.0033

    basis_v: np.ndarray = np.concatenate([basis_G, basis_T, basis_T], axis=0)

    # Choose wavelength and wavenumber
    target_λ: np.float64 = np.float64(8640.0)
    target_k: np.float64 = 2 * np.pi * 4320.0 / target_λ

    kidx: np.int64 = np.argmin(np.abs(target_k - kcal))

    # constrain x domain of inverse matrix
    x_mask: np.ndarray = np.logical_and(x >= -4320000, x <= 4320000)

    inv_mat_lim: np.ndarray = inv_mat[x_mask, :].astype(np.complex64, copy=False)
    x_lim: np.ndarray = x[x_mask]

    # Reconstruct
    FourierBasis: np.ndarray = inv_mat_lim[:, kidx][None, :, None]
    pcs: np.ndarray = state[:, kidx, :][:, None, :]

    profile: np.ndarray = np.einsum(
        "vz,vxt->vzxt", basis_v, FourierBasis * pcs, optimize=True
    )

    print("Finish reconstructing profiles")

    # --------------------------------------------
    # Plot profile
    # --------------------------------------------

    # combining data
    w: np.ndarray = profile[0] + profile[1]
    T: np.ndarray = profile[2] + profile[3]
    J: np.ndarray = profile[4] + profile[5]

    # Design levels

    Tmax: np.float64 = np.nanmax(np.abs(T)) * 0.8 // 0.0005 * 0.0005
    wmax: np.float64 = np.nanmax(np.abs(w)) * 0.8 // 0.01 * 0.01
    Jmax: np.float64 = np.nanmax(np.abs(J)) * 0.8 // 1

    print("Maximum T", np.nanmax(np.abs(T)))

    Tlevel: np.ndarray = np.linspace(-Tmax, Tmax, 21)
    wlevel: np.ndarray = np.linspace(-wmax, wmax, 21)
    Jlevel: np.ndarray = np.linspace(-Jmax, Jmax, 21)

    Tlevel = Tlevel[np.abs(Tlevel) >= 1e-5]
    wlevel = wlevel[np.abs(wlevel) >= 1e-5]
    Jlevel = Jlevel[np.abs(Jlevel) >= 1e-5]

    # Plot animation

    Plot.plot_animation(
        np.linspace(-4320000, 4320000, w.shape[1]),
        np.linspace(0, 14000, w.shape[0]),
        J,
        T,
        w,
        "RdBu_r",
        None,
        Tlevel,
        wlevel,
        "Convective-heating Only",
        fig_dir / "profile_evo.mp4",
        skip=5,
    )

    print("Finish plotting animation")


if __name__ == "__main__":
    main()
