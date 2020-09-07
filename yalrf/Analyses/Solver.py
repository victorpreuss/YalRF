import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import numpy as np

def solve_linear(A, z, is_sparse=False):
    if is_sparse == True:
        An = scipy.sparse.csc_matrix(A)
        lu = scipy.sparse.linalg.splu(An)
        x = lu.solve(z)
    else:
        lu, piv = scipy.linalg.lu_factor(A)
        x = scipy.linalg.lu_solve((lu, piv), z)

    return x, not np.isnan(np.sum(x))