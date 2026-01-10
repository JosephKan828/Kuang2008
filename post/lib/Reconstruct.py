import numpy as np

def _compute_matmul(
    pcs  : np.ndarray,
    basis: np.ndarray
) -> np.ndarray:

    # make pcs and basis into numpy array
    basis_arr: np.ndarray = np.asarray( basis ).reshape( -1 )
    pcs_arr  : np.ndarray = np.asarray( pcs )
    
    return np.einsum( "z,xt->zxt", basis_arr, pcs_arr, optimize=True )

def reconstruct(
    kidx  : np.int64,
    state : np.ndarray,
    basis : np.ndarray,
    invmat:  np.ndarray,
) -> tuple[ np.ndarray, np.ndarray ]:

    FourierBasis: np.ndarray = invmat[ :, kidx ][ None, :, None ]   # ( Nt, )
    state_kidx  : np.ndarray = state[ :, kidx, : ][ :, None, : ]

    state_pc    : np.ndarray = FourierBasis * state_kidx

    prof        : np.ndarray = np.einsum( "vz,vxt->vzxt", basis, state_pc, optimize=True )

    return prof
    # state1      : np.ndarray = state[0][ kidx, : ] # ( Nx, )
    # state2      : np.ndarray = state[1][ kidx, : ] # ( Nx, )


    # # convert coefficient to horizontal PC
    # state1_pc : np.ndarray = FourierBasis[ :, None ] * state1[ None, : ]  # Shape: ( Nx, Nt )
    # state2_pc : np.ndarray = FourierBasis[ :, None ] * state2[ None, : ]

    # # reconstruct profile of data
    # basis1 : np.ndarray = np.asarray( basis[0] ).squeeze()
    # basis2 : np.ndarray = np.asarray( basis[1] ).squeeze()

    # prof1 : np.ndarray = _compute_matmul( state1_pc, basis1 )
    # prof2 : np.ndarray = _compute_matmul( state2_pc, basis2 )

    # return prof1, prof2

