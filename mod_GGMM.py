"""Reconstructed version of Brian's file to construct GGMMs"""
import numpy as np
from qutip import Qobj

def construct_ggmm_sub(dim):
    ''' Constructs Generalized Gell-Mann Matrices for a given dimension Hilbert space 'dim'
        
        Returns:
            ggmm_list (list): A list of Qobj matrices representing the Generalized Gell-Mann matrices.
    '''
    matrices = []
    # Symmetric Gell-Mann matrices (off-diagonal)
    for i in range(dim):
        for j in range(i + 1, dim):
            mat = np.zeros((dim, dim), dtype=complex)
            mat[i, j] = 1
            mat[j, i] = 1
            matrices.append(Qobj(mat))

    # Anti-symmetric Gell-Mann matrices (off-diagonal)
    for i in range(dim):
        for j in range(i + 1, dim):
            mat = np.zeros((dim, dim), dtype=complex)
            mat[i, j] = -1j
            mat[j, i] = 1j
            matrices.append(Qobj(mat))

    # Diagonal Gell-Mann matrices
    for k in range(1, dim):
        mat = np.zeros((dim, dim), dtype=complex)
        for i in range(k):
            mat[i, i] = 1
        mat[k, k] = -k
        mat /= np.sqrt(k * (k + 1))
        matrices.append(Qobj(mat))

    return matrices

# Example usage
# dim = 8
# ggmm_list = construct_ggmm_sub(dim)
# for i, matrix in enumerate(ggmm_list):
#     print(f"Gell-Mann Matrix {i+1}:\n{matrix}\n")


def test_GGMMs(dim):
    # Test for dimension 2
    print('Testing GGMMs...')
    matrices = construct_ggmm_sub(dim)
    expected_count = dim**2 - 1
    assert len(matrices) == expected_count, f"Expected {expected_count} matrices, got {len(matrices)}"

    for i, mat in enumerate(matrices):
        # Check if matrix is traceless
        trace = mat.tr()
        assert np.isclose(trace, 0, atol=1e-10), f"Matrix {i+1} is not traceless, trace = {trace}"

        # Check if the matrix is Hermitian
        hermitian_check = mat == mat.dag()
        assert hermitian_check, f"Matrix {i+1} is not Hermitian."

    print("All tests passed.")

# Run the test function
test_GGMMs(8)

# print('check particular mat')
# mat = construct_ggmm_sub(8)[0]
# hermitian_check = (mat == mat.dag())
# print(hermitian_check)
