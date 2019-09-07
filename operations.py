import math
import numpy as np

identity_mat = np.array([[1, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0],
                         [0, 0, 1, 0, 0],
                         [0, 0, 0, 1, 0],
                         [0, 0, 0, 0, 1]], dtype=np.float64)


def get_identity_mat(size=2):
    identity = np.zeros((size, size), dtype=np.float64)
    for i in range(size):
        for j in range(size):
            if i == j:
                identity[i, j] = 1
            else:
                identity[i, j] = 0

    return identity


def pseudo_inverse(matrix):
    c_mat = copy_matrix(matrix)
    transposed_mat = transpose(c_mat)  # transposing Matrix

    c_transposed_mat = copy_matrix(transposed_mat)

    mul_mat = np.matmul(c_transposed_mat, matrix)  # multiplying transposed mat and original mat
    # mul_mat = np.matmul(matrix, c_transposed_mat)
    inverse_mat = inverse(mul_mat)

    pseudo_mat = np.matmul(inverse_mat, transposed_mat)
    # pseudo_mat = np.matmul(transposed_mat, inverse_mat)

    return round_decimals(pseudo_mat)


def round_decimals(matrix):
    return np.around(matrix, decimals=2)


def copy_matrix(matrix):
    return np.copy(matrix)


def determinant_fast(mat):
    n = len(mat)
    matrix = copy_matrix(mat)

    if n == 2:
        return (matrix[0, 0] * matrix[1, 1]) - (matrix[1, 0] * matrix[0, 1])

    for focus_diagonal in range(n):
        for i in range(focus_diagonal + 1, n):
            if matrix[focus_diagonal][focus_diagonal] == 0:
                matrix[focus_diagonal][focus_diagonal] = 0.0

            try:
                scalar = np.float64(matrix[i][focus_diagonal] / matrix[focus_diagonal][focus_diagonal])
                if math.isnan(scalar):
                    scalar = 0
            except Exception as e:
                print e

            for j in range(n):
                matrix[i][j] = np.float64(matrix[i][j] - (scalar * matrix[focus_diagonal][j]))

    product = 1.0
    for i in range(n):
        product *= matrix[i][i]

    return product


def matrix_minor(matrix, i, j):
    # Verify parameters.

    if i < 0 or i >= matrix.shape[0]:
        raise ValueError("Row value %d is out of range" % i)
    if j < 0 or j >= matrix.shape[1]:
        raise ValueError("Column value %d is out of range" % j)
    # Create the output matrix.
    m = np.zeros((matrix.shape[0] - 1, matrix.shape[1] - 1))
    # Loop through the matrix, skipping over the row and column specified
    # by i and j.
    minor_row = minor_col = 0
    for self_row in xrange(matrix.shape[0]):
        if not self_row == i:  # Skip row i.
            for self_col in xrange(matrix.shape[1]):
                if not self_col == j:  # Skip column j.
                    m[(minor_row, minor_col)] = matrix[(self_row, self_col)]
                    minor_col += 1
            minor_col = 0
            minor_row += 1
    return m


def adjoint_matrix(matrix):
    mat_size = len(matrix)

    adj_mat = np.zeros((mat_size, mat_size))
    for i in range(mat_size):
        for j in range(mat_size):
            co_factor = matrix_minor(matrix, i, j)
            mat = np.asarray(co_factor, dtype=np.float64)
            adj_mat[i, j] = ((-1) ** (i + j)) * determinant_fast(mat)

    return transpose(adj_mat)


def validate_inverse(matrix_a, matrix_b):
    validate_mat = np.matmul(matrix_a, matrix_b)
    validate_mat = np.around(validate_mat, decimals=0)
    return np.equal(validate_mat, get_identity_mat(size=validate_mat.shape[0]))
    # return np.equal(validate_mat, identity_mat)


def inverse(matrix):
    determinant = determinant_fast(m)

    if len(matrix) == 2:
        return [[m[1][1]/determinant, -1 * matrix[0][1]/determinant],
                [-1*m[1][0]/determinant, matrix[0][0]/determinant]]

    cofactors = []
    for row in range(len(matrix)):
        cofactor_row = []
        for col in range(len(m)):
            minor = matrix_minor(matrix, row, col)
            cofactor_row.append(((-1)**(row + col)) * determinant_fast(minor))

        cofactors.append(cofactor_row)
    cofactors = transpose(cofactors)

    for row in range(len(cofactors)):
        for col in range(len(cofactors)):
            cofactors[row][col] = cofactors[row][col] / determinant

    return cofactors


def transpose(matrix):
    mat = copy_matrix(matrix)
    transposed_mat = np.zeros((mat.shape[1], mat.shape[0]))

    for i in range(len(mat)):
        for j in range(len(mat[0])):
            transposed_mat[j, i] = mat[i, j]

    return transposed_mat


# mat = np.array([[5, 4, 3, 2, 1],
#                 [4, 3, 2, 1, 5],
#                 [3, 2, 9, 5, 4],
#                 [2, 1, 5, 4, 3],
#                 [1, 2, 3, 4, 5]], dtype=np.float64)


# mat = np.array([[1, 2],
#                 [3, 4],
#                 [5, 6]])
# mat = np.array([[3, 2, 1, 2],
#                 [7, 5, 2, 5],
#                 [0, 0, 9, 4],
#                 [0, 0, 11, 5]], dtype=np.float64)

mat = np.array([[2, 5, 7],
                [6, 3, 4],
                [5, -2, -3]], dtype=np.float64)

# mat = np.array([[1, 1, 1, 1],
#                 [5, 7, 7, 9]], dtype=np.float64)


m = inverse(mat)
print m
# inv = inverse(mat)
# print validate_inverse(mat, inv)

# mat = [[1, 2, 3], [0, 1, 5], [5, 6, 0]]

# mat = np.array([[1, 2, 3], [0, 4, 5], [1, 0, 6]], dtype=np.float64)
# print determinant_fast(mat)
# mat = np.array([[3, -2, 4],
#                 [2, -4, 5],
#                 [1, 8, 2]], dtype=np.float64)
