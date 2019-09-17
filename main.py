"""Matrix Operation library -
Students:
    Joseph Asaf Gardin: 204614911
    Guy Cohen: 307992958"""


import numpy as np
from operations import determinant_fast, inverse, pseudo_inverse, adjoint_matrix


if __name__ == '__main__':
    mat = np.array([[5, 4, 3, 2, 1],
                    [4, 3, 2, 1, 5],
                    [3, 2, 9, 5, 4],
                    [2, 1, 5, 4, 3],
                    [1, 2, 3, 4, 5]], dtype=np.float64)

    # Non square matrix
    # mat = np.array([[1, 3, 5],
    #                 [1, 3, 1],
    #                 [4, 3, 9],
    #                 [5, 2, 0]], dtype=np.float64)

    print 'Determinant is:\n {}'.format(determinant_fast(mat))
    print
    print 'Adjoint Matrix is:\n {}'.format(adjoint_matrix(mat))
    print
    print 'Inverse Matrix is:\n {}'.format(inverse(mat))
    print
    print 'Pseudo-Inverse Matrix is:\n {}'.format(pseudo_inverse(mat))
    print
