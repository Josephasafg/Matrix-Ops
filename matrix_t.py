import types
import operator

"""Linear Algebra Matrix Class

The Matrix class is an implementation of a linear algebra matrix.  
Arithmetic operations, trace, determinant, and minors are defined for it.  This
is a lightweight alternative to a numerical Python package for people who need
to do basic linear algebra.

Vectors are implemented as 1xN and Nx1 matricies.  There is no separate vector
class.  This implementation enforces the distinction between row and column
vectors.

Indexing is zero-based, i.e. the upper left-hand corner of a matrix is element
(0,0), not element (1,1).

Matricies are stored as a list of lists, where the top level lists are the rows
and the sub-lists are the columns.  Because of the way Python handles list
references, you have be careful when copying matrix objects.  If you have a
matrix a, assign b=a, and then change values in b, you will change values in a
as well.  Matrix copying should be done with copy.deepcopy.

This implementation has no memory-saving optimization for sparse matricies.  A
derived class may implement a more sophisticated storage method by overriding 
the __getitem__ and __setitem__ functions.

Determinants are taken by expanding by minors on the top row.  The private 
functions supplied for expansion by minors are more generic than what is needed
by this implementation.  They may be used by a derived class that wishes to do
more efficient expansion of sparse matricies.

By default, Matrix elements are members of the complex field, but if you want
to perform linear algebra on something other than numbers you may redefine
Matrix.null_element, Matrix.identity_element, and Matrix.inverse_element and 
override the is_scalar_element function.

References:
	George Arfken, "Mathematical Methods for Physicists", 3rd ed. San Diego:
Academic Press Inc. (1985)
"""

__author__ = "Bill McNeill <billmcn@speakeasy.net>"
__version__ = "1.0"


class Matrix_Error(Exception):
    """Abstract parent for all matrix exceptions
	"""
    pass


class Matrix_Arithmetic_Error(Matrix_Error):
    """Incorrect dimensions for arithmetic

	This exception is thrown when you try to add or multiply matricies of
	incompatible sizes.
	"""

    def __init__(self, a, b, operation):
        self.a = a
        self.b = b
        self.operation = operation

    def __str__(self):
        return "Cannot %s a %dx%d and a %dx%d matrix" % \
               (self.operation, \
                self.a.rows(), self.a.cols(), \
                self.b.rows(), self.b.cols())


class Matrix_Multiplication_Error(Matrix_Arithmetic_Error):
    """Thrown when you try to multiply matricies of incompatible dimensions.

	This exception is also thrown when you try to right-multiply a row vector or
	left-multiply a column vector.
	"""

    def __init__(self, a, b):
        Matrix_Arithmetic_Error.__init__(self, a, b, "multiply")


class Matrix_Addition_Error(Matrix_Arithmetic_Error):
    """Thrown when you try to add matricies of incompatible dimensions.
	"""

    def __init__(self, a, b):
        Matrix_Arithmetic_Error.__init__(self, a, b, "add")


class Square_Error(Matrix_Error):
    """Square-matrix only

	This exception is thrown when you try to calculate a function that is only
	defined for square matricies on a non-square matrix.
	"""

    def __init__(self, func):
        self.func = func

    def __str__(self):
        return "%s only defined for square matricies." % self.func


class Trace_Error(Square_Error):
    """Thrown when you try to get the trace of a non-square matrix.
	"""

    def __init__(self):
        Square_Error.__init__(self, "The trace is")


class Minor_Error(Square_Error):
    """Thrown when you try to take a minor of a non-square matrix.
	"""

    def __init__(self):
        Square_Error.__init__(self, "Minors are")


class Determinant_Error(Square_Error):
    """Thrown when you try to take the determinant of a non-square matrix.
	"""

    def __init__(self):
        Square_Error.__init__(self, "The determinant is")


class Matrix:
    """A linear algebra matrix

	This class defines a generic matrix and the basic matrix operations from
	linear algebra.  An instance of this class is a single matrix with
	particular values.
	"""
    null_element = 0
    identity_element = 1
    inverse_element = -1

    def __init__(self, *args):
        pass


    def __str__(self):
        s = ""
        for row in self.m:
            s += "%s\n" % row
        return s

    def __getitem__(self, (row, col)):
        return self.m[row][col]

    def __setitem__(self, (row, col), value):
        self.m[row][col] = value

    def rows(self):
        return len(self.m)

    def cols(self):
        return len(self.m[0])

    def row(self, i):
        return self.m[i]

    def col(self, j):
        r = [row[j] for row in self.m]
        return r

    def scalar_multiply(self, scalar):
        r = []
        for row in self.m:
            r.append(map(lambda x: x * scalar, row))
        return Matrix(r)

    def is_square(self):
        return self.rows() == self.cols()

    def transpose(self):
        r = []
        for col in xrange(self.cols()):
            r.append(self.col(col))
        return Matrix(r)

    def trace(self):
        """The trace of the matrix
        """
        if not self.is_square():
            raise Trace_Error()
        t = 0
        for i in xrange(self.rows()):
            t += self[(i, i)]
        return t

    def determinant(self):
        if not self.is_square():
            raise Determinant_Error()
        # Calculate 2x2 determinants directly.
        if self.rows() == 2:
            return self[(0, 0)] * self[(1, 1)] - self[(0, 1)] * self[(1, 0)]
        # Expand by minors for larger matricies.
        return self.expand_by_minors_on_row(0)

    def expand_by_minors_on_row(self, row):
        """Calculates the determinant by expansion of minors

        This function returns the determinant of the matrix by doing an
        expansion of minors on the specified row.
        """

        assert (row < self.rows())
        d = 0
        for col in xrange(self.cols()):
            # Note: the () around -1 are needed.  Otherwise you get -(1**col).
            d += (-1) ** (row + col) * \
                 self[(row, col)] * self.minor(row, col).determinant()
        return d

    def minor(self, i, j):
        # Verify parameters.
        if not self.is_square():
            raise Minor_Error()
        if i < 0 or i >= self.rows():
            raise ValueError("Row value %d is out of range" % i)
        if j < 0 or j >= self.cols():
            raise ValueError("Column value %d is out of range" % j)
        # Create the output matrix.
        m = Matrix(self.rows() - 1, self.cols() - 1)
        # Loop through the matrix, skipping over the row and column specified
        # by i and j.
        minor_row = minor_col = 0
        for self_row in xrange(self.rows()):
            if not self_row == i:  # Skip row i.
                for self_col in xrange(self.cols()):
                    if not self_col == j:  # Skip column j.
                        m[(minor_row, minor_col)] = self[(self_row, self_col)]
                        minor_col += 1
                minor_col = 0
                minor_row += 1
        return m


def unit_matrix(n):
    m = Matrix(n)
    for i in xrange(m.rows()):
        m[(i, i)] = m.identity_element
    return m
