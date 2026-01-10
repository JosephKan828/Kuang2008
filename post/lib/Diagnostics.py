import numpy as np

def modal_growth_rate(
    optrs: np.ndarray,
    k    : np.ndarray
) -> np.ndarray:

    Nk: np.int64 = int( len( k ) )
    Nv: np.int64 = 6

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

    Nk: np.int64 = int( len( k ) )
    Nv: np.int64 = 6

    length_scale: np.float64 = 4320000.0
    time_scale  : np.float64 = 86400.0

    speed: np.ndarray = np.zeros( ( Nv, Nk ), dtype=np.float64 )

    for j in range( Nk ):
        eigvals, eigvecs = np.linalg.eig( optrs[ :, :, j ].T )

        c = -np.imag( eigvals ) / k[ j ] * ( length_scale / time_scale )

        speed[ :, j ] = c

    return speed





