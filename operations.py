import numpy as np


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


def inverse(matrix):
    pass


# mat = np.array([[3, -2, 4],
#                 [2, -4, 5],
#                 [1, 8, 2]], dtype=np.float64)

A = np.array([[5, 4, 3, 2, 1], [4, 3, 2, 1, 5], [3, 2, 9, 5, 4], [2, 1, 5, 4, 3], [1, 2, 3, 4, 5]], dtype=np.float64)
I = np.array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]], dtype=np.float64)

AM = A.copy()
IM = I.copy()
n = len(AM)

fd = 0  # fd stands for focus diagonal OR the current diagonal
fdScaler = 1. / AM[fd][fd]

for j in range(n):  # using j to indicate cycling thru columns
    AM[fd][j] = fdScaler * AM[fd][j]
    IM[fd][j] = fdScaler * IM[fd][j]

AM = np.around(AM, decimals=3)
IM = np.around(IM, decimals=3)

n = len(A)
indices = list(range(n))

for i in indices[0:fd] + indices[fd + 1:]:  # *** skip row with fd in it.
    crScaler = AM[i][fd]  # cr stands for "current row".
    for j in range(n):  # cr - crScaler * fdRow, but one element at a time.
        AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
        IM[i][j] = IM[i][j] - crScaler * IM[fd][j]
AM = np.around(AM, decimals=3)
IM = np.around(IM, decimals=3)

indices = list(range(n))  # to allow flexible row referencing ***
# We've already run for fd = 0, now let's run for fd = 1 to the last fd
for fd in range(1, n):  # fd stands for focus diagonal
    fdScaler = 1.0 / AM[fd][fd]
    # FIRST: scale fd row with fd inverse.
    for j in range(n):  # Use j to indicate column looping.
        AM[fd][j] *= fdScaler
        IM[fd][j] *= fdScaler

    IM[IM == 0] = 0
    AM[AM == 0] = 0
    AM = np.around(AM, decimals=3)
    IM = np.around(IM, decimals=3)
    # SECOND: operate on all rows except fd row.
    for i in indices[:fd] + indices[fd + 1:]:  # *** skip row with fd in it.
        crScaler = AM[i][fd]  # cr stands for "current row".
        for j in range(n):  # cr - crScaler * fdRow, but one element at a time.
            AM[i][j] = AM[i][j] - crScaler * AM[fd][j]
            IM[i][j] = IM[i][j] - crScaler * IM[fd][j]

        IM[IM == 0] = 0
        AM[AM == 0] = 0
        AM = np.around(AM, decimals=3)
        IM = np.around(IM, decimals=3)

C = np.matmul(A, IM)
C = np.around(C, decimals=0)
print C

# mat = [[1, 2, 3], [0, 1, 5], [5, 6, 0]]

# mat = np.array([[1, 2, 3], [0, 4, 5], [1, 0, 6]], dtype=np.float64)
# print determinant_fast(mat)
