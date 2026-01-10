# ==================================================================================
# Compare growth rate between no_rad and qt_rad experiments with rad_scaling = 0.005
# ==================================================================================

import sys
import h5py
import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt

sys.path.append(
    "/home/b11209013/2025_Research/Kuang2008/post/lib"
)

import Diagnostics #type: ignore

def main(
        ROOT: Path,
        case1: str,
        case2: str,
        rad_scaling: str
) -> None:
    
    # ----------------------------------
    # Determine paths
    # ----------------------------------
    
    INPUT: Path  = ROOT / "output"
    OUTPUT: Path = ROOT / "analysis" / "figures"
    
    OUTPUT.mkdir( parents=True, exist_ok=True )

    # ----------------------------------
    # load data
    # ----------------------------------

    fnames: list[ Path ] = []

    for case in [ case1, case2 ]:
        if case == "no_rad":
            fnames.append( INPUT / case / "optrs.h5" )
        else:
            fnames.append( INPUT / case / f"rad_scaling={rad_scaling}" / "optrs.h5" )

    # Load wavenumbers
    with h5py.File( INPUT / case1 / "state.h5", "r" ) as f:
        wnum: np.ndarray = f.get( "wavenumber" )[ ... ]

    # Load operators
    optrs: dict[ str, np.ndarray ] = {
        case1: h5py.File( fnames[0], "r" )[ "operators" ][ ... ],
        case2: h5py.File( fnames[1], "r" )[ "operators" ][ ... ]
    }
    
    # ----------------------------------
    # Compute growth rates
    # ----------------------------------

    σ: dict[ str, np.ndarray ] = {
        case: Diagnostics.modal_growth_rate( optrs[case], wnum )
        for case in [ case1, case2 ]
    }

    # ----------------------------------
    # Compute maximum growth rate
    # ----------------------------------

    max_idxs: dict[ str, np.ndarray ] = {
        case: np.argmax( σ[case], axis=0 )
        for case in [ case1, case2 ]
    }

    max_σs: dict[ str, np.ndarray ] = {
        case: np.take_along_axis( σ[case], max_idxs[case][ None, : ], axis=0 )[0]
        for case in [ case1, case2 ]
    }

    # ----------------------------------
    # Plot growth rates
    # ----------------------------------

    # set axis
    λ       : np.ndarray = 2.0 * np.pi * 4320.0 / wnum
    k_nondim: np.ndarray = 40000 / λ

    # color list
    clist: list[str] = [ "#085993", "#009C24" ]
    llist: list[str] = [ "no radiation", "with radiation" ]

    fig = plt.figure( figsize=( 10.5, 6.2 ) )

    for i, case in enumerate( [ case1, case2 ] ):
        plt.scatter(
            k_nondim, max_σs[ case ],
            linewidth=2.5, color=clist[i], label=llist[i]
        )
        plt.gca().spines[ "top" ].set_visible( False )
        plt.gca().spines[ "right" ].set_visible( False )

    plt.xlabel( "Zonal wavenumber", fontsize=18 )
    plt.ylabel( r"Growth Rate [ day$^{-1}$ ]", fontsize=18 )
    plt.xticks( np.linspace(0, 30, 7), fontsize=16 )
    plt.yticks( fontsize=16 )
    plt.xlim( 0, 30 )
    plt.ylim( 0, None )
    plt.legend( fontsize=16, frameon=False )

    plt.tight_layout()
    plt.savefig( OUTPUT / f"rad_growth_rate_comp.png", dpi=600, bbox_inches="tight" )
    plt.close()


if __name__ == "__main__":
    # experiment parameters
    case1: str = "no_rad"
    case2: str = "qt_rad"
    rad_scaling: str = "0.005"
    
    # setup input and output paths
    ROOT = Path( __file__ ).resolve().parents[ 1 ]

    main( ROOT, case1, case2, rad_scaling )