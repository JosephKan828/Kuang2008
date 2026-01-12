from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")  # Use non-interactive backend
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

def plot_growth_rate(
    growth_rate: np.ndarray,
    wnum       : np.ndarray,
    path
) -> None:

    # Calculate wavelength
    λ: np.ndarray = 2.0 * np.pi * 4320.0 / wnum

    # Find most unstable mode
    max_idx: np.ndarray = np.argmax( growth_rate, axis=0 )
    # max_σ  : np.ndarray = np.asarray([
    #     growth_rate[max_idx[i], i]
    #     for i in range(wnum.size)
    # ])

    max_σ: np.ndarray = np.take_along_axis(growth_rate, max_idx[None, :], axis=0)[0]

    # Plot  
    x = 40000 / λ

    fig = plt.figure( figsize=( 10.5, 6.2 ) )

    plt.scatter(
        x, max_σ,
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

    plt.savefig( path, dpi=600, bbox_inches="tight" )
    plt.close()

def plot_diagnostics(
    growth_rate : np.ndarray,
    phase_speed : np.ndarray,
    wnum        : np.ndarray,
    path_growth ,
    path_speed
) -> None:

    # Calculate wavelength
    λ: np.ndarray = 2.0 * np.pi * 4320.0 / wnum

    # Find most unstable mode
    max_idx: np.ndarray = np.argmax( growth_rate, axis=0 )
    # max_c  : np.ndarray = np.asarray([
    #     phase_speed[max_idx[i], i]
    #     for i in range(wnum.size)
    # ])
    max_σ: np.ndarray = np.take_along_axis(growth_rate, max_idx[None, :], axis=0)[0]
    max_c: np.ndarray = np.take_along_axis(phase_speed, max_idx[None, :], axis=0)[0]

    # Plot  
    x = 40000 / λ

    fig = plt.figure( figsize=( 10.5, 6.2 ) )

    for j in range( phase_speed.shape[0] ):
        plt.scatter(
            x, phase_speed[ j, : ],
            s=3, c="blue", alpha=0.3
        )

    plt.scatter(
        x, max_c,
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
        x, max_σ,
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
    draw_contours: bool = True
) -> None:
    """
    Scientific animation with full redraw per frame.
    No use of contour collections or artist mutation.
    """

    # ---- Cast to real ----
    J = np.real(J)
    T = np.real(T)
    w = np.real(w)

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

    # Precompute which array updater pcolormesh uses
    # QuadMesh stores flattened array of size nz*nx.
    def _set_pcm(field2d: np.ndarray) -> None:
        pcm.set_array(field2d.ravel())

    def update(i: int):
        nonlocal cont_T, cont_w


        pcm.set_array(J[:, :, i].ravel())
        # _set_pcm(J[:, :, i])

        # if draw_contours:
        #     # remove previous contours (much cheaper than ax.clear())
        #     if cont_T is not None:
        #         for artist in cont_T.collections: # type: ignore
        #             artist.remove()
        #     if cont_w is not None:
        #         for artist in cont_w.collections: # type: ignore
        #             artist.remove()

        #     cont_T = ax.contour(x, z, T[:, :, i], levels=T_level, colors="k", linewidths=1.8)
        #     cont_w = ax.contour(x, z, w[:, :, i], levels=w_level, colors="seagreen", linewidths=1.8)

        # ax.set_title(f"{title}  |  frame {i+1}/{nt}", fontsize=20)

        # # blit=True works well if draw_contours=False (otherwise contours recreate artists each frame)
        # return (pcm,)
    
        # Remove old contours (NO .collections usage)
        if cont_T is not None:
            cont_T.remove()
        if cont_w is not None:
            cont_w.remove()

        # Draw new contours
        cont_T = ax.contour(x, z, T[:, :, i], levels=T_level, colors="k", linewidths=2, zorder=2)
        cont_w = ax.contour(x, z, w[:, :, i], levels=w_level, colors="seagreen", linewidths=2, zorder=2)

        return (pcm,)  # keep blit simple; contours are re-created

    # blit_ok = not draw_contours  # contours break the benefit of blitting
    ani = FuncAnimation(fig, update, frames=frames, interval=1000 // fps, blit=False)

    # Encoding settings also affect time a lot.
    writer = FFMpegWriter(
        fps=fps,
        bitrate=5_000,
        codec="libx264",
        extra_args=["-pix_fmt", "yuv420p", "-preset", "veryfast",  "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",]
    )
    ani.save(path, writer=writer, dpi=150)  # dpi=150 usually enough for videos

    plt.close(fig)