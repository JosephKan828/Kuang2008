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
) -> np.ndarray:

    FourierBasis: np.ndarray = invmat[ :, kidx ][ None, :, None ]   # ( Nt, )
    state_kidx  : np.ndarray = state[ :, kidx, : ][ :, None, : ]

    state_pc    : np.ndarray = FourierBasis * state_kidx

    prof        : np.ndarray = np.einsum( "vz,vxt->vzxt", basis, state_pc, optimize=True )

    return prof