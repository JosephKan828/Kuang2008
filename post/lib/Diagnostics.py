import numpy as np

def modal_growth_rate(
    optrs: np.ndarray,
    k    : np.ndarray
) -> np.ndarray:

    Nk: int = int( len( k ) )
    Nv: int = 6

    growth: np.ndarray = np.zeros( ( Nv, Nk ), dtype=np.float64 )

    for j in range( Nk ):
        eigvals, eigvecs = np.linalg.eig( optrs[ :, :, j ].T )

        σ = np.real( eigvals )

        growth[ :, j ] = σ

    return growth

def phase_speed(
    optrs: np.ndarray,
    k    : np.ndarray
) -> np.ndarray:

    Nk: int = int( len( k ) )
    Nv: int = 6

    length_scale: float = 4320000.0
    time_scale  : float = 86400.0

    speed: np.ndarray = np.zeros( ( Nv, Nk ), dtype=np.float64 )

    for j in range( Nk ):
        eigvals, eigvecs = np.linalg.eig( optrs[ :, :, j ].T )

        c = -np.imag( eigvals ) / k[ j ] * ( length_scale / time_scale )

        speed[ :, j ] = c

    return speed

def convert_state(
        state: np.ndarray
) -> np.ndarray:

    # original state shape: ( Nv, Nk, Nt )
    γq = 0.7

    L: np.ndarray = state[ -1, :, : ]                                    # low-level heating
    U: np.ndarray = L + 0.7*( state[ -2, :, : ] - 1.5*state[ 2, :, : ] ) # upper-level heating

    J1: np.ndarray = ( L + U )[ None, :, : ]
    J2: np.ndarray = ( L - U )[ None, :, : ]

    state_new = np.concatenate( [ state[:4], J1, J2 ], axis=0 ).astype( np.complex64, copy=False )

    return state_new