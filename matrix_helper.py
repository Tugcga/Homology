from fractions import Fraction
from math import gcd
import functools


MATRIX_INCONSISTENT_DIMENSIONS = 0
MATRIX_NOT_CHAIN_MAPS = 1


def get_size(matrix):
    return (len(matrix), len(matrix[0]))


def matrix_to_string(matrix):
    str_array = []
    for row_index in range(len(matrix)):
        row_array = []
        for v in matrix[row_index]:
            row_array.append(str(v))
        str_array.append("|" + " ".join(row_array) + "|")
    return "\n".join(str_array)


def multiply_matrices(m1, m2, size1, size2):
    to_return = [[0 for j in range(size2[1])] for i in range(size1[0])]
    for i in range(size1[0]):
        for j in range(size2[1]):
            for k in range(size1[1]):
                to_return[i][j] += m1[i][k] * m2[k][j]
    return to_return


def build_constant_matrix(dimension, value=0, is_tuples=False):
    '''Generate matrix of the size (dimension[0], dimension[1]) with elements are equal to value

    is_tuples = Treu means that result should be tuples, not lists
    '''
    if is_tuples:
        return tuple(tuple(value for j in range(dimension[1])) for i in range(dimension[0]))
    else:
        return [[value for j in range(dimension[1])] for i in range(dimension[0])]


def copy_matrix(matrix):
    return [[v for v in row] for row in matrix]


def add_columns_to_matrix(matrix, array):
    to_return = []
    for row_index in range(len(matrix)):
        new_row = [v for v in matrix[row_index]]
        for a in range(len(array)):
            new_row.append(array[a][row_index])
        to_return.append(new_row)
    return to_return


def add_rows_to_matrix(matrix, array):
    to_return = copy_matrix(matrix)
    for a in array:
        to_return.append([v for v in a])
    return to_return


def build_transpose(matrix):
    return tuple(tuple(matrix[j][i] for j in range(len(matrix))) for i in range(len(matrix[0])))


def orders_to_string(orders):
    '''Present tuple of group orders as direct sum of Z_p and Z.

    For example, if orders = (1, 4, 0), then the method returns Z_4 + Z.
    '''
    str_array = []
    for v in orders:
        if v == 0:
            str_array.append("Z")
        elif v > 1:
            str_array.append("Z_" + str(v))
    if len(str_array) > 0:
        return " + ".join(str_array)
    else:
        return "0"


def int_lcm(a, b):
    return a * b // gcd(a, b)


def lcm_n(*args):
    return functools.reduce(int_lcm, args)


def get_matrix_kernel(matrix):
    '''Return the basis of the matrix kernel.

    The result is an array of tuples. Each tuple has the same length as the count of columns in the original matrix
    '''
    # make a copy of the matrix and convert it values to rational fractions
    m = [[Fraction(v, 1) for v in r] for r in matrix]
    # start gauss process, the number of steps is min(rows, columns)
    column_index = 0  # the indexes of the upper left corner of the considered submatrix
    row_index = 0
    while column_index < len(m[0]) and row_index < len(m):
        # try to find the first row with non-zero element in the column column_index with (rows start from row_index)
        is_find = False
        i = row_index
        while not is_find:
            if m[i][column_index] != 0:
                is_find = True
            else:
                i += 1
            if i >= len(m):
                is_find = True
        if i >= len(m):
            # no non-zero elements in the column, siply increase it
            column_index += 1
            # and start next iteration
        else:
            # swap i-th row with row_index row
            # do it only for elements in column_index+ columns
            a = m[i][column_index]
            for j in range(column_index, len(m[row_index])):
                m[row_index][j], m[i][j] = m[i][j], m[row_index][j]
                # divide all elements in the row_index row to the first one
                m[row_index][j] = m[row_index][j] / a
            # zerofy all elements in the column_index colum
            for j in range(len(m)):
                # iterate throw rows
                if j != row_index:
                    a = m[j][column_index]
                    for k in range(column_index, len(m[j])):
                        m[j][k] = m[j][k] - a * m[row_index][k]
            # increase row and column indexes
            row_index += 1
            column_index += 1
    # next we should generate basis of the kernel
    # devide all variables to main and passive
    main_vars = []
    passive_vars = []
    i = 0
    j = 0
    while i < len(m) and j < len(m[0]):
        if m[i][j] != 0:
            main_vars.append(j)
            i += 1
            j += 1
        else:
            passive_vars.append(j)
            j += 1
    for k in range(j, len(m[0])):
        passive_vars.append(k)
    kernel = []
    for i in range(len(passive_vars)):
        var_index = passive_vars[i]
        vector = [0 for j in range(len(m[0]))]
        # calculate var_index value as LCM of all denominators in the corresponding column
        vector[var_index] = lcm_n(*[m[j][var_index].denominator for j in range(len(m))])
        # calculate all main variables
        for j in range(len(main_vars)):
            main_var_index = main_vars[j]
            vector[main_var_index] = -1 * vector[var_index] * m[j][var_index]

        # add vector to the answer
        kernel.append(tuple(int(v) for v in vector))
    return kernel


def reduce(matrix):
    '''Apply Gauss process to the matrix. Use only integer numbers and use row transforms
    '''
    # create copy of the matrix
    def is_zero_row(row):
        for v in row:
            if v != 0:
                return False
        return True

    m = [[v for v in row] for row in matrix]
    height = len(m)  # the number of rows in the matrix
    width = len(m[0])  # the number of columns in the matrix
    row_index = 0
    column_index = 0
    while row_index < height and column_index < width:
        # find the minimum number (by absolute value) in the column column_index in the rows after row_index
        c_min = None
        c_min_index = -1
        for i in range(row_index, height):
            if m[i][column_index] != 0:
                v = abs(m[i][column_index])
                if c_min is None:
                    c_min = v
                    c_min_index = i
                else:
                    if v < c_min:
                        c_min = v
                        c_min_index = i
        if c_min_index == -1:
            # there are no non-zero elements in the column, start next step
            column_index += 1
        else:
            if c_min_index != row_index:
                # switch two rows
                for k in range(column_index, width):
                    m[c_min_index][k], m[row_index][k] = m[row_index][k], m[c_min_index][k]
            if m[row_index][column_index] < 0:
                for k in range(column_index, width):
                    m[row_index][k] *= -1
            # decrease values in the column
            a = abs(m[row_index][column_index])
            for i in range(height):
                # iterate by all rows
                if i != row_index:
                    b = abs(m[i][column_index])
                    if b >= a:
                        # apply transform to the row i
                        c = (-1 if m[row_index][column_index] * m[i][column_index] > 0 else 1) * (b // a)
                        for j in range(column_index, width):
                            m[i][j] = m[i][j] + c * m[row_index][j]
            # check, are all elements in the column now = 0
            is_zeros = True
            is_finish = False
            i = row_index + 1
            while not is_finish:
                if i >= height:
                    is_finish = True
                else:
                    if m[i][column_index] != 0:
                        is_zeros = False
                        is_finish = True
                    i += 1
            if is_zeros:
                row_index += 1
                column_index += 1
    # remove zero rows
    to_return = []

    for i in range(height):
        if not is_zero_row(m[i]):
            to_return.append(tuple(v for v in m[i]))
    return to_return


def get_error_string(error_code):
    if error_code == MATRIX_INCONSISTENT_DIMENSIONS:
        return "Inconsistend dimensions of matrices."
    elif error_code == MATRIX_NOT_CHAIN_MAPS:
        return "Matrix pair is not present chain map."
    else:
        return "Unknown error."
