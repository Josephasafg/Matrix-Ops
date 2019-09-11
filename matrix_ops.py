import numpy as np


def transposeMatrix(m):
    return map(list,zip(*m))


def getMatrixMinor(m,i,j):
    r = []
    for row in m[i+1:]:
        r.append(row[:j] + row[j+1:])

    print
    # return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]


def minor(Matrix, i, j):
    """A minor of the matrix

    This function returns the minor given by striking out row i and
    column j of the matrix.
    """
    # Verify parameters.

    if i < 0 or i >= Matrix.shape[0]:
        raise ValueError("Row value %d is out of range" % i)
    if j < 0 or j >= Matrix.shape[1]:
        raise ValueError("Column value %d is out of range" % j)
    # Create the output matrix.
    m = np.zeros((Matrix.shape[0] - 1, Matrix.shape[1] - 1))
    # Loop through the matrix, skipping over the row and column specified
    # by i and j.
    minor_row = minor_col = 0
    for self_row in xrange(Matrix.shape[0]):
        if not self_row == i:  # Skip row i.
            for self_col in xrange(Matrix.shape[1]):
                if not self_col == j:  # Skip column j.
                    m[(minor_row, minor_col)] = Matrix[(self_row, self_col)]
                    minor_col += 1
            minor_col = 0
            minor_row += 1
    return m


def getMatrixDeternminant(m):
    #base case for 2x2 matrix
    if len(m) == 2:
        return m[0][0]*m[1][1]-m[0][1]*m[1][0]

    determinant = 0
    for c in range(len(m)):
        determinant += ((-1)**c)*m[0][c]*getMatrixDeternminant(minor(m,0,c))
    return determinant


def getMatrixInverse(m):
    determinant = getMatrixDeternminant(m)
    #special case for 2x2 matrix:
    if len(m) == 2:
        return [[m[1][1]/determinant, -1*m[0][1]/determinant],
                [-1*m[1][0]/determinant, m[0][0]/determinant]]

    #find matrix of cofactors
    cofactors = []
    for r in range(len(m)):
        cofactorRow = []
        for c in range(len(m)):
            mino = minor(m,r,c)
            cofactorRow.append(((-1)**(r+c)) * getMatrixDeternminant(mino))
        cofactors.append(cofactorRow)
    cofactors = transposeMatrix(cofactors)
    for r in range(len(cofactors)):
        for c in range(len(cofactors)):
            cofactors[r][c] = cofactors[r][c]/determinant
    return cofactors


# def inverse(matrix):
#     work_matrix = copy_matrix(matrix)
#     identity_mat_c = copy_matrix(get_identity_mat(size=work_matrix.shape[0]))
#
#     n = len(work_matrix)
#
#     focus_diagonal = 0
#     fdScaler = 1. / work_matrix[focus_diagonal][focus_diagonal]
#
#     for j in range(n):
#         work_matrix[focus_diagonal][j] = fdScaler * work_matrix[focus_diagonal][j]
#         identity_mat_c[focus_diagonal][j] = fdScaler * identity_mat_c[focus_diagonal][j]
#
#     work_matrix = np.around(work_matrix, decimals=3)
#     identity_mat_c = np.around(identity_mat_c, decimals=3)
#
#     n = len(matrix)
#     indices = list(range(n))
#
#     for i in indices[0:focus_diagonal] + indices[focus_diagonal + 1:]:
#         current_row_scalar = work_matrix[i][focus_diagonal]
#         for j in range(n):
#             work_matrix[i][j] = work_matrix[i][j] - current_row_scalar * work_matrix[focus_diagonal][j]
#             identity_mat_c[i][j] = identity_mat_c[i][j] - current_row_scalar * identity_mat_c[focus_diagonal][j]
#
#     work_matrix = np.around(work_matrix, decimals=3)
#     identity_mat_c = np.around(identity_mat_c, decimals=3)
#
#     indices = list(range(n))
#     for fd in range(1, n):
#         fdScaler = 1.0 / work_matrix[fd][fd]
#         for j in range(n):
#             work_matrix[fd][j] *= fdScaler
#             identity_mat_c[fd][j] *= fdScaler
#
#         identity_mat_c[identity_mat_c == 0] = 0
#         work_matrix[work_matrix == 0] = 0
#         work_matrix = np.around(work_matrix, decimals=3)
#         identity_mat_c = np.around(identity_mat_c, decimals=3)
#
#         for i in indices[:fd] + indices[fd + 1:]:
#             current_row_scalar = work_matrix[i][fd]
#             for j in range(n):
#                 work_matrix[i][j] = work_matrix[i][j] - current_row_scalar * work_matrix[fd][j]
#                 identity_mat_c[i][j] = identity_mat_c[i][j] - current_row_scalar * identity_mat_c[fd][j]
#
#             identity_mat_c[identity_mat_c == 0] = 0
#             work_matrix[work_matrix == 0] = 0
#             work_matrix = np.around(work_matrix, decimals=3)
#             identity_mat_c = np.around(identity_mat_c, decimals=3)
#
#     return identity_mat_c


mat = np.array([[5, 4, 3, 2, 1],
                [4, 3, 2, 1, 5],
                [3, 2, 9, 5, 4],
                [2, 1, 5, 4, 3],
                [1, 2, 3, 4, 5]], dtype=np.float64)

# mat = np.array([[1, 4, 7],
#                 [3, 0, 5],
#                 [-1, 9, 11]], dtype=np.float64)


m = getMatrixInverse(mat)
print m