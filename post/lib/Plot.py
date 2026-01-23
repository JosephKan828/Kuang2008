from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")  # Use non-interactive backend
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

def plot_diagnostics(
    max_σ : np.ndarray,
    max_c : np.ndarray,
    phase_speed : np.ndarray,
    wnum        : np.ndarray,
    path_growth ,
    path_speed
) -> None:

    # Plot

    fig = plt.figure( figsize=( 10.5, 6.2 ) )

    for j in range( phase_speed.shape[0] ):
        plt.scatter(
            wnum, phase_speed[ j, : ],
            s=3, c="blue", alpha=0.3
        )

    plt.scatter(
        wnum, max_c,
        s=6, c="white", edgecolor="black", 
    )
    plt.gca().spines[ "top" ].set_visible( False )
    plt.gca().spines[ "right" ].set_visible( False )

    plt.xlabel( "Zonal Wavenumber", fontsize=18 )
    plt.ylabel( r"Phase Speed [ m s$^{-1}$ ]", fontsize=18 )
    plt.xticks( np.linspace(0, 30, 7), fontsize=16 )
    plt.yticks( fontsize=16 )
    plt.xlim( 0, 30 )
    plt.ylim( 0, None )


    plt.savefig( path_speed, dpi=600, bbox_inches="tight" )
    plt.close()

    fig = plt.figure( figsize=( 10.5, 6.2 ) )

    plt.scatter(
        wnum, max_σ,
        linewidth=2.5, color="black"
    )
    plt.gca().spines[ "top" ].set_visible( False )
    plt.gca().spines[ "right" ].set_visible( False )

    plt.xlabel( "Zonal Wavenumber", fontsize=18 )
    plt.ylabel( r"Growth Rate [ day$^{-1}$ ]", fontsize=18 )
    plt.xticks( np.linspace(0, 30, 7), fontsize=16 )
    plt.yticks( fontsize=16 )
    plt.xlim( 0, 30 )
    plt.ylim( 0, None )
    plt.text( 
        18, 0.8 * np.max( max_σ ),
        f"Max Growth Rate : {np.max( max_σ ):0.3f} day$^{{-1}}$ \nat Wavenumber : {wnum[ np.argmax( max_σ ) ]:0.2f}",
        fontsize=16 )

    plt.savefig( path_growth, dpi=600, bbox_inches="tight" )
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
    skip: int = 1
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
    ax.set_xticklabels(["-40", "-30", "-20", "-10", "0", "10", "20", "30", "40"], fontsize=16)

    ax.set_yticks(np.linspace(0, 14_000, 8))
    ax.set_yticklabels(["0", "2", "4", "6", "8", "10", "12", "14"], fontsize=16)

    ax.set_xlabel("X [100 km]", fontsize=18)
    ax.set_ylabel("Z [km]", fontsize=18)
    ax.minorticks_on()
    ax.set_aspect("auto")

# ---- Create pcolormesh ONCE ----
    pcm = ax.pcolormesh(
        x, z, J[:, :, 0],
        cmap=cmap, vmin=vmin, vmax=vmax,
        shading="auto"
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
        cont_T = curr_ax.contour(x, z, T[:, :, i], levels=T_level, colors="k", linewidths=2, zorder=2)
        cont_w = curr_ax.contour(x, z, w[:, :, i], levels=w_level, colors="seagreen", linewidths=2, zorder=2)

        return (pcm,)  # keep blit simple; contours are re-created

    # blit_ok = not draw_contours  # contours break the benefit of blitting
    ani = FuncAnimation(fig, update, frames=frames, interval=1000 // fps, blit=False)

    # Encoding settings also affect time a lot.
    writer = FFMpegWriter(
        fps=fps,
        bitrate=5_000,
        codec="libx264",
        extra_args=["-pix_fmt", "yuv420p", "-preset", "ultrafast",  "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",]
    )
    ani.save(path, writer=writer, dpi=150)  # dpi=150 usually enough for videos

    plt.close(fig)