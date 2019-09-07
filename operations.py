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

            scalar = np.float64(matrix[i][focus_diagonal] / matrix[focus_diagonal][focus_diagonal])

            for j in range(n):
                matrix[i][j] = np.float64(matrix[i][j] - (scalar * matrix[focus_diagonal][j]))

    product = 1.0
    for i in range(n):
        product *= matrix[i][i]

    return product


def get_matrix_minor(m, i, j):
    result = [row[:j] + row[j + 1:] for row in (m[:i] + m[i + 1:])]
    return result


def adjoint_matrix(matrix):
    mat_size = len(matrix)

    adj_mat = np.zeros((mat_size, mat_size))
    for i in range(mat_size):
        for j in range(mat_size):
            co_factor = get_matrix_minor(matrix, i, j)
            mat = np.asarray(co_factor, dtype=np.float64)
            adj_mat[i, j] = ((-1) ** (i + j)) * determinant_fast(mat)

    return transpose(adj_mat)


def validate_inverse(matrix_a, matrix_b):
    validate_mat = np.matmul(matrix_a, matrix_b)
    validate_mat = np.around(validate_mat, decimals=0)
    return np.equal(validate_mat, identity_mat)


def inverse(matrix):
    work_matrix = copy_matrix(matrix)
    identity_mat_c = copy_matrix(get_identity_mat(size=work_matrix.shape[0]))

    n = len(work_matrix)

    focus_diagonal = 0
    fdScaler = 1. / work_matrix[focus_diagonal][focus_diagonal]

    for j in range(n):
        work_matrix[focus_diagonal][j] = fdScaler * work_matrix[focus_diagonal][j]
        identity_mat_c[focus_diagonal][j] = fdScaler * identity_mat_c[focus_diagonal][j]

    work_matrix = np.around(work_matrix, decimals=3)
    identity_mat_c = np.around(identity_mat_c, decimals=3)

    n = len(matrix)
    indices = list(range(n))

    for i in indices[0:focus_diagonal] + indices[focus_diagonal + 1:]:
        current_row_scalar = work_matrix[i][focus_diagonal]
        for j in range(n):
            work_matrix[i][j] = work_matrix[i][j] - current_row_scalar * work_matrix[focus_diagonal][j]
            identity_mat_c[i][j] = identity_mat_c[i][j] - current_row_scalar * identity_mat_c[focus_diagonal][j]
    work_matrix = np.around(work_matrix, decimals=3)
    identity_mat_c = np.around(identity_mat_c, decimals=3)

    indices = list(range(n))
    for fd in range(1, n):
        fdScaler = 1.0 / work_matrix[fd][fd]
        for j in range(n):
            work_matrix[fd][j] *= fdScaler
            identity_mat_c[fd][j] *= fdScaler

        identity_mat_c[identity_mat_c == 0] = 0
        work_matrix[work_matrix == 0] = 0
        work_matrix = np.around(work_matrix, decimals=3)
        identity_mat_c = np.around(identity_mat_c, decimals=3)
        for i in indices[:fd] + indices[fd + 1:]:
            current_row_scalar = work_matrix[i][fd]
            for j in range(n):
                work_matrix[i][j] = work_matrix[i][j] - current_row_scalar * work_matrix[fd][j]
                identity_mat_c[i][j] = identity_mat_c[i][j] - current_row_scalar * identity_mat_c[fd][j]

            identity_mat_c[identity_mat_c == 0] = 0
            work_matrix[work_matrix == 0] = 0
            work_matrix = np.around(work_matrix, decimals=3)
            identity_mat_c = np.around(identity_mat_c, decimals=3)

    return identity_mat_c


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

mat = np.array([[1, 1, 1, 1],
                [5, 7, 7, 9]], dtype=np.float64)

print pseudo_inverse(mat)
# inv = inverse(mat)
# print validate_inverse(mat, inv)

# mat = [[1, 2, 3], [0, 1, 5], [5, 6, 0]]

# mat = np.array([[1, 2, 3], [0, 4, 5], [1, 0, 6]], dtype=np.float64)
# print determinant_fast(mat)
# mat = np.array([[3, -2, 4],
#                 [2, -4, 5],
#                 [1, 8, 2]], dtype=np.float64)
