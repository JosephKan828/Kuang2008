from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")  # Use non-interactive backend
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter


def _apply_scientific_style(ax, xlabel, ylabel):
    """Internal helper to apply consistent styling to plots."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlabel(xlabel, fontsize=14, fontweight="medium")
    ax.set_ylabel(ylabel, fontsize=14, fontweight="medium")

    # Set tick parameters (inside pointing ticks are common in physics)
    ax.tick_params(axis="both", labelsize=12, direction="in", top=False, right=False)
    ax.grid(True, linestyle="--", alpha=0.4)  # Subtle grid for readability

    ax.set_xlim(0, 30)
    ax.set_ylim(bottom=0)


def plot_diagnostics(
    max_σ: np.ndarray,
    max_c: np.ndarray,
    phase_speed: np.ndarray,
    wnum: np.ndarray,
    path_growth,
    path_speed,
) -> None:

    # --- 1. Phase speed ---

    fig, ax = plt.subplots(figsize=(8, 5))

    for j in range(phase_speed.shape[0]):
        ax.scatter(wnum, phase_speed[j, :], color="gray", alpha=0.15, s=3)

    ax.plot(wnum, max_c, color="#1f77b4", linewidth=2, label="Most Unstable Mode")
    ax.scatter(wnum, max_c, s=15, color="#1f77b4", edgecolor="white")

    _apply_scientific_style(
        ax, r"Zonal Wavenumber ($k$)", r"Phase Speed $c$ [m s$^{-1}$]"
    )

    plt.savefig(path_speed, dpi=300, bbox_inches="tight")
    plt.close()

    # --- 2. Growth rate ---
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(wnum, max_σ, color="firebrick", linewidth=2.5)
    ax.fill_between(wnum, max_σ, color="firebrick", alpha=0.1)

    _apply_scientific_style(ax, "Zonal Wavenumber", "Growth Rate [ day$^{-1}$ ]")

    # Add annotation for max growth rate
    idx_max = np.argmax(max_σ)
    max_val = max_σ[idx_max]
    max_k = wnum[idx_max]

    annotation_text = (
        rf"$\sigma_{{max}} = {max_val:.3f}$ d$^{{-1}}$"
        + "\n"
        + rf"$k_{{max}} = {max_k:.2f}$"
    )
    ax.annotate(
        annotation_text,
        xy=(20, max_val),
        xytext=(20 + 2, max_val * 0.9),
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8),
        # arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"),
        fontsize=12,
    )

    plt.savefig(path_growth, dpi=600, bbox_inches="tight")
    plt.close()


def plot_animation(
    x,
    z,
    J,
    T,
    w,
    cmap,
    J_level,
    T_level,
    w_level,
    title,
    path,
    frames: int = 300,
    fps: int = 40,
    skip: int = 1,
) -> None:
    """
    Scientific animation with full redraw per frame.
    No use of contour collections or artist mutation.
    """

    x = x[::skip]
    z = z[::skip]

    # ---- Cast to real ----
    J = np.real(J)[::skip, ::skip, :]
    T = np.real(T)[::skip, ::skip, :]
    w = np.real(w)[::skip, ::skip, :]

    nz, nx, nt = J.shape

    # ---- Fixed color scaling (critical for science) ----
    if J_level is None:
        vmin, vmax = np.nanmin(J), np.nanmax(J)
    else:
        vmin, vmax = J_level

    # ---- Figure ----
    fig, ax = plt.subplots(figsize=(10.5, 6.2), constrained_layout=True)

    ax.set_xlim(-4_320_000, 4_320_000)
    ax.set_ylim(0, 14_000)

    ax.set_xticks(np.linspace(-4_000_000, 4_000_000, 9))
    ax.set_xticklabels(
        ["-40", "-30", "-20", "-10", "0", "10", "20", "30", "40"], fontsize=16
    )

    ax.set_yticks(np.linspace(0, 14_000, 8))
    ax.set_yticklabels(["0", "2", "4", "6", "8", "10", "12", "14"], fontsize=16)

    ax.set_xlabel("X [100 km]", fontsize=18)
    ax.set_ylabel("Z [km]", fontsize=18)
    ax.minorticks_on()
    ax.set_aspect("auto")

    # ---- Create pcolormesh ONCE ----
    pcm = ax.pcolormesh(
        x, z, J[:, :, 0], cmap=cmap, vmin=vmin, vmax=vmax, shading="auto"
    )

    cbar = fig.colorbar(pcm, ax=ax, orientation="horizontal", shrink=0.85, aspect=45)
    cbar.ax.tick_params(labelsize=14)

    # ---- Contours: keep references so we can remove them ----
    cont_T = None
    cont_w = None

    curr_ax = ax

    def update(i: int):
        nonlocal cont_T, cont_w

        pcm.set_array(J[:, :, i].ravel())

        # Remove old contours (NO .collections usage)
        if cont_T is not None:
            cont_T.remove()
        if cont_w is not None:
            cont_w.remove()

        # Draw new contours
        cont_T = curr_ax.contour(
            x, z, T[:, :, i], levels=T_level, colors="k", linewidths=2, zorder=2
        )
        cont_w = curr_ax.contour(
            x, z, w[:, :, i], levels=w_level, colors="seagreen", linewidths=2, zorder=2
        )

        return (pcm,)  # keep blit simple; contours are re-created

    # blit_ok = not draw_contours  # contours break the benefit of blitting
    ani = FuncAnimation(fig, update, frames=frames, interval=1000 // fps, blit=False)

    # Encoding settings also affect time a lot.
    writer = FFMpegWriter(
        fps=fps,
        bitrate=5_000,
        codec="libx264",
        extra_args=[
            "-pix_fmt",
            "yuv420p",
            "-preset",
            "ultrafast",
            "-vf",
            "scale=trunc(iw/2)*2:trunc(ih/2)*2",
        ],
    )
    ani.save(path, writer=writer, dpi=300)  # dpi=150 usually enough for videos

    plt.close(fig)
