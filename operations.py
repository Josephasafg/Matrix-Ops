import numpy as np

identity_mat = np.array([[1, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0],
                         [0, 0, 1, 0, 0],
                         [0, 0, 0, 1, 0],
                         [0, 0, 0, 0, 1]], dtype=np.float64)


def copy_matrix(matrix):
    return np.copy(matrix)


def determinant_fast(A):
    # Section 1: Establish n parameter and copy A
    n = len(A)
    matrix = copy_matrix(A)
    if n == 2:
        return (matrix[0, 0] * matrix[1, 1]) - (matrix[1, 0] * matrix[0, 1])

    # Section 2: Row ops on A to get in upper triangle form
    for fd in range(n):  # A) fd stands for focus diagonal
        for i in range(fd + 1, n):  # B) only use rows below fd row
            if matrix[fd][fd] == 0:  # C) if diagonal is zero ...
                matrix[fd][fd] = 0.0  # change to ~zero
            # D) cr stands for "current row"
            crScaler = np.float64(matrix[i][fd] / matrix[fd][fd])
            # E) cr - crScaler * fdRow, one element at a time
            for j in range(n):
                matrix[i][j] = np.float64(matrix[i][j] - (crScaler * matrix[fd][j]))

    # Section 3: Once matrix is in upper triangle form ...
    product = 1.0
    for i in range(n):
        # ... product of diagonals is determinant
        product *= matrix[i][i]

    return product


def get_matrix_minor(m, i, j):
    result = [row[:j] + row[j + 1:] for row in (m[:i] + m[i + 1:])]
    return result


def adjoint_matrix(matrix):
    mat_size = len(matrix)

    # mat = np.zeros((0,0))
    adj_mat = np.zeros((mat_size, mat_size))
    for i in range(mat_size):
        for j in range(mat_size):
            cofactor = get_matrix_minor(matrix, i, j)
            mat = np.asarray(cofactor, dtype=np.float64)
            adj_mat[i, j] = ((-1) ** (i + j)) * determinant_fast(mat)

    return np.transpose(adj_mat)  #: TODO: write transpose function


def validate_inverse(matrix_a, matrix_b):
    validate_mat = np.matmul(matrix_a, matrix_b)
    validate_mat = np.around(validate_mat, decimals=0)
    return np.equal(validate_mat, identity_mat)


def inverse(matrix):
    work_matrix = matrix.copy()
    identity_mat_c = identity_mat.copy()
    n = len(work_matrix)

    fd = 0  # fd stands for focus diagonal OR the current diagonal
    fdScaler = 1. / work_matrix[fd][fd]

    for j in range(n):  # using j to indicate cycling thru columns
        work_matrix[fd][j] = fdScaler * work_matrix[fd][j]
        identity_mat_c[fd][j] = fdScaler * identity_mat_c[fd][j]

    work_matrix = np.around(work_matrix, decimals=3)
    identity_mat_c = np.around(identity_mat_c, decimals=3)

    n = len(matrix)
    indices = list(range(n))

    for i in indices[0:fd] + indices[fd + 1:]:  # *** skip row with fd in it.
        crScaler = work_matrix[i][fd]  # cr stands for "current row".
        for j in range(n):  # cr - crScaler * fdRow, but one element at a time.
            work_matrix[i][j] = work_matrix[i][j] - crScaler * work_matrix[fd][j]
            identity_mat_c[i][j] = identity_mat_c[i][j] - crScaler * identity_mat_c[fd][j]
    work_matrix = np.around(work_matrix, decimals=3)
    identity_mat_c = np.around(identity_mat_c, decimals=3)

    indices = list(range(n))  # to allow flexible row referencing ***
    # We've already run for fd = 0, now let's run for fd = 1 to the last fd
    for fd in range(1, n):  # fd stands for focus diagonal
        fdScaler = 1.0 / work_matrix[fd][fd]
        # Fidentity_matRST: scale fd row with fd inverse.
        for j in range(n):  # Use j to indicate column looping.
            work_matrix[fd][j] *= fdScaler
            identity_mat_c[fd][j] *= fdScaler

        identity_mat_c[identity_mat_c == 0] = 0
        work_matrix[work_matrix == 0] = 0
        work_matrix = np.around(work_matrix, decimals=3)
        identity_mat_c = np.around(identity_mat_c, decimals=3)
        # SECOND: operate on all rows except fd row.
        for i in indices[:fd] + indices[fd + 1:]:  # *** skip row with fd in it.
            crScaler = work_matrix[i][fd]  # cr stands for "current row".
            for j in range(n):  # cr - crScaler * fdRow, but one element at a time.
                work_matrix[i][j] = work_matrix[i][j] - crScaler * work_matrix[fd][j]
                identity_mat_c[i][j] = identity_mat_c[i][j] - crScaler * identity_mat_c[fd][j]

            identity_mat_c[identity_mat_c == 0] = 0
            work_matrix[work_matrix == 0] = 0
            work_matrix = np.around(work_matrix, decimals=3)
            identity_mat_c = np.around(identity_mat_c, decimals=3)

    return identity_mat_c


mat = np.array([[5, 4, 3, 2, 1],
                [4, 3, 2, 1, 5],
                [3, 2, 9, 5, 4],
                [2, 1, 5, 4, 3],
                [1, 2, 3, 4, 5]], dtype=np.float64)

inv = inverse(mat)
print validate_inverse(mat, inv)

# mat = [[1, 2, 3], [0, 1, 5], [5, 6, 0]]

# mat = np.array([[1, 2, 3], [0, 4, 5], [1, 0, 6]], dtype=np.float64)
# print determinant_fast(mat)
# mat = np.array([[3, -2, 4],
#                 [2, -4, 5],
#                 [1, 8, 2]], dtype=np.float64)
