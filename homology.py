import random


class MatrixProcessor(object):
    def __init__(self, matrix_first, matrix_second):
        '''matrix_first is matrix for [d_n], matrix_second for [d_{n + 1}]
        '''
        self._m1 = [[v for v in row] for row in matrix_first]
        self._m2 = [[v for v in row] for row in matrix_second]
        (self._correct, self._incorrect_mode) = self.check_correct()
        if self._correct is False:
            print("Input matrices are incorrect. " + ("M_1 * M_2 \\neq 0." if self._incorrect_mode == 1 else "Inconsistend dimensions."))
        self._m1_base = len(self._m1[0])
        self._transform = [[1 if i == j else 0 for j in range(self._m1_base)] for i in range(self._m1_base)]
        self._reducted = False  # turn to True after reduction pair of matrices
        self._data_diagonal = None  # contains values in the main diagonal (with shifted rows) of the M2 after reduction

    def get_is_correct(self):
        return self._correct

    def get_first_matrix(self):
        return self._m1

    def get_second_matrix(self):
        return self._m2

    def _get_size(self, matrix):
        return (len(matrix), len(matrix[0]))

    def _multiply_matrices(self, m1, m2, size1, size2):
        to_return = [[0 for j in range(size2[1])] for i in range(size1[0])]
        for i in range(size1[0]):
            for j in range(size2[1]):
                for k in range(size1[1]):
                    to_return[i][j] += m1[i][k] * m2[k][j]
        return to_return

    def check_correct(self):
        '''Is matrices, presented in processor, are boundary maps, i.e. M1 * M2 = 0

        Return pair (bool, int), where bool is True for correct matrix pair, int is key for incorrect result
        key = -1: correct result
        key = 0: wrong dimensions
        key = 1: M1 * M2 is not equal to zero matrix
        '''
        # check dimensions
        m1_dim = self._get_size(self._m1)
        m2_dim = self._get_size(self._m2)
        if m1_dim[1] == m2_dim[0]:
            # make multiplication
            mult = self._multiply_matrices(self._m1, self._m2, m1_dim, m2_dim)
            # is all elements are zero
            s = self._get_size(mult)
            for i in range(s[0]):
                for j in range(s[1]):
                    if mult[i][j] != 0:
                        return (False, 1)
            return (True, -1)
        else:
            return (False, 0)

    def _is_row_column_zero(self, matrix, dim, min_row, min_column):
        '''Return True if all elements in min_row and min_column are zero
        '''
        for i in range(min_column + 1, dim[1]):
            if matrix[min_row][i] != 0:
                return False
        for i in range(min_row + 1, dim[0]):
            if matrix[i][min_column] != 0:
                return False
        return True

    def _get_minimal_coordinates(self, matrix, dim, min_row, min_colum):
        to_return = (-1, -1)
        min_value = 0
        for i in range(min_row, dim[0]):
            for j in range(min_colum, dim[1]):
                v = abs(matrix[i][j])
                if v != 0:
                    if (min_value == 0 and v != 0) or (min_value != 0 and v < min_value):
                        min_value = v
                        to_return = (i, j)
        return to_return

    def _switch_transfrom_value(self, i, j, switch_i, switch_j):
        if i == j:
            if i != switch_i and j != switch_j:
                return 1
            else:
                return 0
        else:
            if (i == switch_i and j == switch_j) or (i == switch_j and j == switch_i):
                return 1
            else:
                return 0

    def _negate_transform_value(self, i, j, negate_i):
        if i == j:
            if i == negate_i:
                return -1
            else:
                return 1
        else:
            return 0

    def _add_transform_value(self, i, j, from_i, to_j, coefficient):
        if i == j:
            return 1
        else:
            if i == from_i and j == to_j:
                return coefficient
            else:
                return 0

    def _switch_columns(self, matrix, dim, min_row, column_i, column_j, second_matrix=None, second_dim=(0, 0)):
        '''If second_matrix is not None, then switch row of it
        '''
        for k in range(min_row, dim[0]):
            c = matrix[k][column_i]
            matrix[k][column_i] = matrix[k][column_j]
            matrix[k][column_j] = c
        if second_matrix is not None:
            for k in range(second_dim[1]):
                c = second_matrix[column_i][k]
                second_matrix[column_i][k] = second_matrix[column_j][k]
                second_matrix[column_j][k] = c
            # if second_matrix is not None, then we switch columns of the M1 matrix
            self._transform = self._multiply_matrices(self._transform, [[self._switch_transfrom_value(i, j, column_i, column_j) for j in range(self._m1_base)] for i in range(self._m1_base)], (self._m1_base, self._m1_base), (self._m1_base, self._m1_base))

    def _switch_row(self, matrix, dim, min_column, row_i, row_j, second_matrix=None, second_dim=(0, 0)):
        '''If second matrix is not None, switch columns of it
        '''
        for k in range(min_column, dim[1]):
            c = matrix[row_i][k]
            matrix[row_i][k] = matrix[row_j][k]
            matrix[row_j][k] = c
        if second_matrix is not None:
            for k in range(second_dim[0]):
                c = second_matrix[k][row_i]
                second_matrix[k][row_i] = second_matrix[k][row_j]
                second_matrix[k][row_j] = c
            self._transform = self._multiply_matrices(self._transform, [[self._switch_transfrom_value(i, j, row_i, row_j) for j in range(self._m1_base)] for i in range(self._m1_base)], (self._m1_base, self._m1_base), (self._m1_base, self._m1_base))

    def _negate_row(self, matrix, dim, min_column, row_i, secondary_matrix=None, secondary_dim=(0, 0)):
        for k in range(min_column, dim[1]):
            matrix[row_i][k] *= -1
        if secondary_matrix is not None:
            for k in range(secondary_dim[0]):
                secondary_matrix[k][row_i] *= -1

    def _negate_column(self, matrix, dim, min_row, column_i, secondary_matrix=None, secondary_dim=(0, 0)):
        for k in range(min_row, dim[0]):
            matrix[k][column_i] *= -1
        if secondary_matrix is not None:
            for k in range(secondary_dim[1]):
                secondary_matrix[column_i][k] *= -1
            self._transform = self._multiply_matrices(self._transform, [[self._negate_transform_value(i, j, column_i) for j in range(self._m1_base)] for i in range(self._m1_base)], (self._m1_base, self._m1_base), (self._m1_base, self._m1_base))

    def _add_row(self, matrix, dim, min_column, row_i, row_j, coefficient, write_transform=False):
        '''Add row_i to the row_j with coefficient c
        If write_transform=True, then we should save transformation in transform matrix, because we make it fro M2 and change the basis of the M1
        '''
        for k in range(min_column, dim[1]):
            matrix[row_j][k] += coefficient * matrix[row_i][k]
        if write_transform:
            self._transform = self._multiply_matrices(self._transform, [[self._add_transform_value(i, j, row_j, row_i, -1*coefficient) for j in range(self._m1_base)] for i in range(self._m1_base)], (self._m1_base, self._m1_base), (self._m1_base, self._m1_base))

    def _add_column(self, matrix, dim, min_row, column_i, column_j, coefficient, secondary_matrix=None, secondary_dim=(0, 0)):
        '''Add column_i to the column_j with coefficient c
        For secondary matrix subtract from i-th row the j-th row with coefficient c
        '''
        for k in range(min_row, dim[0]):
            matrix[k][column_j] += coefficient*matrix[k][column_i]
        if secondary_matrix is not None:
            for k in range(secondary_dim[1]):
                secondary_matrix[column_i][k] -= coefficient*secondary_matrix[column_j][k]
            self._transform = self._multiply_matrices(self._transform, [[self._add_transform_value(i, j, column_i, column_j, coefficient) for j in range(self._m1_base)] for i in range(self._m1_base)], (self._m1_base, self._m1_base), (self._m1_base, self._m1_base))

    def _get_diagonal_length(self, matrix, dim):
        to_return = 0
        size = min(dim[0], dim[1])
        for k in range(size):
            if matrix[k][k] != 0:
                to_return += 1
            else:
                return to_return
        return to_return

    def make_reduction(self):
        '''Apply reduction process to matrix pair for calculating homology group

        We make the following transformations:
        1. Integer transform of rows of the M1
        2. Integer transform of columns of the M2
        3. Simultaneous transfrom of columns of the M1 and rows of the M2.
        '''
        m1_dim = self._get_size(self._m1)
        m2_dim = self._get_size(self._m2)
        # 1. Reduce M1. Done it in several steps. Each step zerofy one row and colum
        steps_count = min(m1_dim[0], m1_dim[1])
        for step in range(steps_count):
            # for each step consider submatrix of the M1 from (step, step) element
            # check if we should zerofy row or colum
            while not self._is_row_column_zero(self._m1, m1_dim, step, step):
                # step-row or step-column are not zero, we should zerofy it
                # find minimal non-zero value in submatrix
                pos = self._get_minimal_coordinates(self._m1, m1_dim, step, step)
                # transpose (step, step)-element with pos-element, if we need it
                if pos[0] != step:
                    self._switch_row(self._m1, m1_dim, step, step, pos[0])
                if pos[1] != step:
                    self._switch_columns(self._m1, m1_dim, step, step, pos[1], self._m2, m2_dim)
                # reduce rows of the M1 (in fact step column), without effect on columns of the M2
                # make corner elements positive
                if self._m1[step][step] < 0:
                    self._negate_row(self._m1, m1_dim, step, step)
                a = self._m1[step][step]
                for k in range(step + 1, m1_dim[0]):
                    v = self._m1[k][step]
                    if v != 0:
                        c = abs(v) // a
                        # reduce k-th row by step-row with coefficient -c
                        self._add_row(self._m1, m1_dim, step, step, k, c if v < 0 else -1*c)
                # reduce step row of the M1 (by summand columns) with effect on rows of the M2
                for k in range(step + 1, m1_dim[1]):
                    v = self._m1[step][k]
                    if v != 0:
                        c = abs(v) // a
                        self._add_column(self._m1, m1_dim, step, step, k, c if v < 0 else -1*c, self._m2, m2_dim)
        # 2. Next reduce M2 without any effects on the M1
        # Find the number of non-zero elements in the diagonal of the M1
        shift = self._get_diagonal_length(self._m1, m1_dim)
        # start from row shift + 1 and column 0
        steps_count = min(m2_dim[1], m2_dim[0] - shift)
        for step in range(steps_count):
            # repeat the same process
            while not self._is_row_column_zero(self._m2, m2_dim, step + shift, step):
                pos = self._get_minimal_coordinates(self._m2, m2_dim, step + shift, step)
                if pos[0] != step + shift:
                    self._switch_row(self._m2, m2_dim, step, step + shift, pos[0])
                if pos[1] != step:
                    self._switch_columns(self._m2, m2_dim, step + shift, step, pos[1])
                if self._m2[step + shift][step] < 0:
                    self._negate_row(self._m2, m2_dim, step, step + shift)
                a = self._m2[step + shift][step]
                for k in range(step + shift + 1, m2_dim[0]):
                    v = self._m2[k][step]
                    if v != 0:
                        c = abs(v) // a
                        self._add_row(self._m2, m2_dim, step, step + shift, k, c if v < 0 else -1*c, True)
                for k in range(step + 1, m2_dim[1]):
                    v = self._m2[step + shift][k]
                    if v != 0:
                        c = abs(v) // a
                        self._add_column(self._m2, m2_dim, step + shift, step, k, c if v < 0 else -1*c)
        # last - make all non-zero elements of the M2 positive
        data_array = []
        for i in range(steps_count):
            if self._m2[shift + i][i] < 0:
                self._negate_column(self._m2, m2_dim, shift + i, i)
            data_array.append(self._m2[shift + i][i])
        self._reducted = True
        self._data_diagonal = tuple(data_array)

    def matrix_to_string(self, matrix):
        str_array = []
        for row_index in range(len(matrix)):
            row_array = []
            for v in matrix[row_index]:
                row_array.append(str(v))
            str_array.append("|" + " ".join(row_array) + "|")
        return "\n".join(str_array)

    def __repr__(self):
        l = 0
        for r in self._m1:
            l = max(l, len(str(r)) - 2 - 2*(len(r) - 1))
        return "\n".join([self.matrix_to_string(self._m1), "-"*(2*l + 1), self.matrix_to_string(self._m2)])


class ChainComplex(object):
    def __init__(self):
        pass


def get_matrix_pair():
    dim1 = (3, 3)
    dim2 = (3, 2)
    limit = 2
    is_find = False
    while not is_find:
        # generate ranom matrix
        matrix1 = [[random.randint(-1*limit, limit) for j in range(dim1[1])] for i in range(dim1[0])]
        matrix2 = [[random.randint(-1*limit, limit) for j in range(dim2[1])] for i in range(dim2[0])]
        # check correctness
        mp = MatrixProcessor(matrix1, matrix2)
        if mp.get_is_correct():
            is_find = True
            return (matrix1, matrix2)

if __name__ == "__main__":
    '''mp = MatrixProcessor([[1, 1, 1, 1], [-1, -1, -1, -1]], [[1, 1, 1], [1, -1, -1], [-1, -1, 1], [-1, 1, -1]])
    if mp.get_is_correct():
        mp.make_reduction()
        print(mp)
        print(mp.matrix_to_string(mp._transform))'''
    cc = ChainComplex()
    
