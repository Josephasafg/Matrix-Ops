# coding=utf-8
import math
import numpy as np
from numpy.linalg import LinAlgError


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
    if isinstance(matrix, np.ndarray):
        try:
            c_mat = copy_matrix(matrix)
            if matrix.shape[0] == matrix.shape[1]:
                if determinant_fast(c_mat) != 0:
                    return inverse(c_mat)

            if is_row_linear_dependant(c_mat):
                return pseudo_inverse_one(c_mat)
            else:
                transposed_mat = transpose(c_mat)
                c_transposed_mat = copy_matrix(transposed_mat)
                mul_mat = np.matmul(c_transposed_mat, matrix)
                inverse_mat = inverse(mul_mat)

                pseudo_mat = np.matmul(inverse_mat, transposed_mat)
                return round_decimals(pseudo_mat)  # We've decided to round the result of the elements for a cleaner
                # output. Of course it can be removed for a more accurate result
                # for example -0.09796705 will now be -0.1
        except ValueError as e:
            print e


"""if rows of a matrix are linear dependant we use the following formula 𝐴^𝑇*(𝐴^𝑇*𝐴)^(−1) for pseudo inverse"""
def pseudo_inverse_one(matrix):
    c_mat = copy_matrix(matrix)
    transposed_mat = transpose(c_mat)
    multiplied_mat = np.matmul(transposed_mat, c_mat)
    inversed_mat = inverse(multiplied_mat)
    return np.matmul(transposed_mat, inversed_mat)


def round_decimals(matrix):
    return np.around(matrix, decimals=2)


def copy_matrix(matrix):
    return np.copy(matrix)


#  By upper triangular calculation
def determinant_fast(mat):
    if mat.__class__ == np.ndarray:
        if mat.shape[0] != mat.shape[1]:
            raise ValueError('Matrix must be a square matrix of nXn size')

        n = len(mat)
        matrix = copy_matrix(mat)

        if n == 2:
            return (matrix[0, 0] * matrix[1, 1]) - (matrix[1, 0] * matrix[0, 1])

        sign_change = -1
        sign_change_amount = 0
        for focus_diagonal in range(n):
            for i in range(focus_diagonal + 1, n):
                if matrix[focus_diagonal][focus_diagonal] == 0:
                    matrix[[focus_diagonal, focus_diagonal + 1]] = matrix[[focus_diagonal+1, focus_diagonal]]
                    sign_change_amount += 1
                    continue

                try:
                    scalar = np.float64(matrix[i][focus_diagonal] / matrix[focus_diagonal][focus_diagonal])
                    if math.isnan(scalar):
                        scalar = 0

                    for j in range(n):
                        matrix[i][j] = np.float64(matrix[i][j] - (scalar * matrix[focus_diagonal][j]))
                except Exception:
                    raise ZeroDivisionError('Cant divide in zero')

        diag_product = 1.0
        for i in range(n):
            diag_product *= matrix[i][i]

        sign_change = sign_change ** sign_change_amount
        return round_decimals(diag_product) * sign_change
    else:
        raise ValueError('Argument must be of type np.ndarray and of nXn size')


def matrix_minor(matrix, i, j):
    # Verify parameters.

    if i < 0 or i >= matrix.shape[0]:
        raise ValueError("Row value %d is out of range" % i)
    if j < 0 or j >= matrix.shape[1]:
        raise ValueError("Column value %d is out of range" % j)

    # Create the output matrix.
    output_mat = np.zeros((matrix.shape[0] - 1, matrix.shape[1] - 1))
    # Loop through the matrix, skipping over the row and column specified
    # by i and j.
    minor_row = minor_col = 0
    for self_row in xrange(matrix.shape[0]):
        if not self_row == i:  # Skip row i.
            for self_col in xrange(matrix.shape[1]):
                if not self_col == j:  # Skip column j.
                    output_mat[(minor_row, minor_col)] = matrix[(self_row, self_col)]
                    minor_col += 1
            minor_col = 0
            minor_row += 1
    return output_mat


def adjoint_matrix(matrix):
    if isinstance(matrix, np.ndarray):
        if matrix.shape[0] != matrix.shape[1]:
            raise ValueError('Matrix must be a square matrix of nXn size')

        mat_size = len(matrix)

        adj_mat = np.zeros((mat_size, mat_size))
        for i in range(mat_size):
            for j in range(mat_size):
                co_factor = matrix_minor(matrix, i, j)
                adj_mat[i, j] = ((-1) ** (i + j)) * determinant_fast(co_factor)

        return transpose(adj_mat)


def validate_inverse(matrix_a, matrix_b):
    if matrix_a.shape[1] != matrix_b.shape[0]:
        raise ValueError('Can\'t multiply matrices due to wrong sizes.')

    validate_mat = np.matmul(matrix_a, matrix_b)
    validate_mat = np.around(validate_mat, decimals=0)
    return np.equal(validate_mat, get_identity_mat(size=validate_mat.shape[0]))


def inverse(matrix):
    if isinstance(matrix, np.ndarray):
        if matrix.shape[0] != matrix.shape[1]:
            raise ValueError('Matrix must be a square matrix of nXn size')

    m = copy_matrix(matrix)
    determinant = determinant_fast(m)

    if determinant == 0:
        raise ValueError('Matrix is not invertible')

    if len(matrix) == 2:
        return [[m[1][1]/determinant, -1 * matrix[0][1]/determinant],
                [-1*m[1][0]/determinant, matrix[0][0]/determinant]]

    adjoint_mat = []
    for row in range(len(matrix)):
        adjoint_row = []
        for col in range(len(m)):
            minor = matrix_minor(matrix, row, col)
            adjoint_row.append(((-1)**(row + col)) * determinant_fast(minor))

        adjoint_mat.append(adjoint_row)
    cofactors = transpose(adjoint_mat)

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


def is_row_linear_dependant(matrix):
    try:
        row, vectors = np.linalg.eig(matrix.T)
        if not matrix[row == 0, :].any():
            return False
    except LinAlgError:
        return False

    return True
