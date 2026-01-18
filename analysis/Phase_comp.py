# ========================================
# Phase analysis
# ========================================

# import package
import h5py
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path

def main() -> None:
    # ----------
    # Load data
    # ----------

    # Input path of output
    INPUT_PATH: Path = Path( __file__ ).parent.parent / "output"

    # load no_rad data
    states: dict[ str, np.ndarray ] = {}

    with h5py.File(INPUT_PATH / "no_rad" / "state.h5", "r") as f:
        states[ "no_rad" ] = f.get( "state" )[ ... ]       #type: ignore
        time: np.ndarray   = f.get( "time" )[ ... ]        #type: ignore
        wnum: np.ndarray   = f.get( "wavenumber" )[ ... ]  #type: ignore
        vars: np.ndarray   = f.get( "variables" )[ ... ]   #type: ignore

    target_λ: float = 8640.0
    target_k: float = 2.0 * np.pi * 4320.0 / target_λ
    kidx   : int   = np.argmin( np.abs( wnum - target_k ) ).astype( int )

    states[ "no_rad" ] = states["no_rad"][ kidx, ... ]

    # load qt_rad, scaling=0.001
    with h5py.File(INPUT_PATH / "qt_rad" / "rad_scaling=0.001" / "state.h5", "r") as f:
        states[ "qt_rad/scaling=0.001" ] = f.get( "state" )[ kidx, ... ] #type: ignore

    # load qt_rad, scaling=0.005
    with h5py.File(INPUT_PATH / "qt_rad" / "rad_scaling=0.005" / "state.h5", "r") as f:
        states[ "qt_rad/scaling=0.005" ] = f.get( "state" )[ kidx, ... ] #type: ignore

    # load qt_rad, scaling=0.01
    with h5py.File(INPUT_PATH / "qt_rad" / "rad_scaling=0.01" / "state.h5", "r") as f:
        states[ "qt_rad/scaling=0.01" ] = f.get( "state" )[ kidx, ... ] #type: ignore
    
    # ---------------------------------
    # Visualize time evolution of phase
    # ---------------------------------
    
    # Plot for different variables
    fig = plt.figure( figsize=(10, 6) )
    for i, var in enumerate( vars ):
        plt.plot( time, states[ "no_rad" ][i].real, label=f"{var}" )
    plt.xlabel( "Time (day)" )
    plt.show()









if __name__ == "__main__":
    main()