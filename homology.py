from fractions import Fraction
from math import gcd
from functools import reduce


MATRIX_INCONSISTENT_DIMENSIONS = 0
MATRIX_NOT_CHAIN_MAPS = 1


class MatrixHelper(object):
    def __init__(self):
        pass

    def get_size(self, matrix):
        return (len(matrix), len(matrix[0]))

    def matrix_to_string(self, matrix):
        str_array = []
        for row_index in range(len(matrix)):
            row_array = []
            for v in matrix[row_index]:
                row_array.append(str(v))
            str_array.append("|" + " ".join(row_array) + "|")
        return "\n".join(str_array)

    def multiply_matrices(self, m1, m2, size1, size2):
        to_return = [[0 for j in range(size2[1])] for i in range(size1[0])]
        for i in range(size1[0]):
            for j in range(size2[1]):
                for k in range(size1[1]):
                    to_return[i][j] += m1[i][k] * m2[k][j]
        return to_return

    def build_constant_matrix(self, dimension, value=0, is_tuples=False):
        '''Generate matrix of the size (dimension[0], dimension[1]) with elements are equal to value

        is_tuples = Treu means that result should be tuples, not lists
        '''
        if is_tuples:
            return tuple(tuple(value for j in range(dimension[1])) for i in range(dimension[0]))
        else:
            return [[value for j in range(dimension[1])] for i in range(dimension[0])]

    def copy_matrix(self, matrix):
        return [[v for v in row] for row in matrix]

    def add_columns_to_matrix(self, matrix, array):
        to_return = []
        for row_index in range(len(matrix)):
            new_row = [v for v in matrix[row_index]]
            for a in range(len(array)):
                new_row.append(array[a][row_index])
            to_return.append(new_row)
        return to_return

    def add_rows_to_matrix(self, matrix, array):
        to_return = self.copy_matrix(matrix)
        for a in array:
            to_return.append([v for v in a])
        return to_return

    def build_transpose(self, matrix):
        return tuple(tuple(matrix[j][i] for j in range(len(matrix))) for i in range(len(matrix[0])))

    def orders_to_string(self, orders):
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

    def lcm(self, a, b):
        return a * b // gcd(a, b)

    def lcm_n(self, *args):
        return reduce(self.lcm, args)

    def get_matrix_kernel(self, matrix):
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
            vector[var_index] = self.lcm_n(*[m[j][var_index].denominator for j in range(len(m))])
            # calculate all main variables
            for j in range(len(main_vars)):
                main_var_index = main_vars[j]
                vector[main_var_index] = -1 * vector[var_index] * m[j][var_index]

            # add vector to the answer
            kernel.append(tuple(int(v) for v in vector))
        return kernel

    def get_error_string(self, error_code):
        if error_code == MATRIX_INCONSISTENT_DIMENSIONS:
            return "Inconsistend dimensions of matrices."
        elif error_code == MATRIX_NOT_CHAIN_MAPS:
            return "Matrix pair is not present chain map."
        else:
            return "Unknown error."


class MatrixProcessor(object):
    def __init__(self, matrix_first, matrix_second, log_error=True):
        '''matrix_first is matrix for [d_n], matrix_second for [d_{n + 1}]
        '''
        self._helper = MatrixHelper()
        self._m1 = [[v for v in row] for row in matrix_first]
        self._m2 = [[v for v in row] for row in matrix_second]
        (self._correct, self._incorrect_code) = self.check_correct()
        if self._correct is False and log_error:
            print(self._helper.get_error_string(self._incorrect_code))
        self._m1_base = len(self._m1[0])
        self._transform = [[1 if i == j else 0 for j in range(self._m1_base)] for i in range(self._m1_base)]
        self._shift = 0  # shift is equal to the number of bases elements, which are not in ker(M1)
        self._reducted = False  # turn to True after reduction pair of matrices
        self._data_diagonal = None  # contains values in the main diagonal (with shifted rows) of the M2 after reduction

    def get_is_correct(self):
        return self._correct

    def get_incorrect_code(self):
        return self._incorrect_code

    def get_is_reducted(self):
        return self._reducted

    def get_first_matrix(self):
        return self._m1

    def get_second_matrix(self):
        return self._m2

    def get_group_orders(self):
        return self._data_diagonal

    def get_transform_column(self, index):
        return tuple(self._transform[i][index] for i in range(len(self._transform)))

    def get_nonkernel_count(self):
        return self._shift

    def check_correct(self):
        '''Is matrices, presented in processor, are boundary maps, i.e. M1 * M2 = 0

        Return pair (bool, int), where bool is True for correct matrix pair, int is key for incorrect result
        key = -1: correct result
        key = 0: wrong dimensions
        key = 1: M1 * M2 is not equal to zero matrix
        '''
        # check dimensions
        m1_dim = self._helper.get_size(self._m1)
        m2_dim = self._helper.get_size(self._m2)
        if m1_dim[1] == m2_dim[0]:
            # make multiplication
            mult = self._helper.multiply_matrices(self._m1, self._m2, m1_dim, m2_dim)
            # is all elements are zero
            s = self._helper.get_size(mult)
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
            self._transform = self._helper.multiply_matrices(self._transform, [[self._switch_transfrom_value(i, j, column_i, column_j) for j in range(self._m1_base)] for i in range(self._m1_base)], (self._m1_base, self._m1_base), (self._m1_base, self._m1_base))

    def _switch_row(self, matrix, dim, min_column, row_i, row_j, second_matrix=None, second_dim=(0, 0), save_to_transform=False):
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
        if second_matrix is not None or save_to_transform:
            self._transform = self._helper.multiply_matrices(self._transform, [[self._switch_transfrom_value(i, j, row_i, row_j) for j in range(self._m1_base)] for i in range(self._m1_base)], (self._m1_base, self._m1_base), (self._m1_base, self._m1_base))

    def _negate_row(self, matrix, dim, min_column, row_i, secondary_matrix=None, secondary_dim=(0, 0), save_to_transform=False):
        for k in range(min_column, dim[1]):
            matrix[row_i][k] *= -1
        if secondary_matrix is not None:
            for k in range(secondary_dim[0]):
                secondary_matrix[k][row_i] *= -1
        if save_to_transform:
            self._transform = self._helper.multiply_matrices(self._transform, [[self._negate_transform_value(i, j, row_i) for j in range(self._m1_base)] for i in range(self._m1_base)], (self._m1_base, self._m1_base), (self._m1_base, self._m1_base))

    def _negate_column(self, matrix, dim, min_row, column_i, secondary_matrix=None, secondary_dim=(0, 0)):
        for k in range(min_row, dim[0]):
            matrix[k][column_i] *= -1
        if secondary_matrix is not None:
            for k in range(secondary_dim[1]):
                secondary_matrix[column_i][k] *= -1
            self._transform = self._helper.multiply_matrices(self._transform, [[self._negate_transform_value(i, j, column_i) for j in range(self._m1_base)] for i in range(self._m1_base)], (self._m1_base, self._m1_base), (self._m1_base, self._m1_base))

    def _add_row(self, matrix, dim, min_column, row_i, row_j, coefficient, save_to_transform=False):
        '''Add row_i to the row_j with coefficient c
        If save_to_transform=True, then we should save transformation in transform matrix, because we make it fro M2 and change the basis of the M1
        '''
        for k in range(min_column, dim[1]):
            matrix[row_j][k] += coefficient * matrix[row_i][k]
        if save_to_transform:
            self._transform = self._helper.multiply_matrices(self._transform, [[self._add_transform_value(i, j, row_j, row_i, -1*coefficient) for j in range(self._m1_base)] for i in range(self._m1_base)], (self._m1_base, self._m1_base), (self._m1_base, self._m1_base))

    def _add_column(self, matrix, dim, min_row, column_i, column_j, coefficient, secondary_matrix=None, secondary_dim=(0, 0)):
        '''Add column_i to the column_j with coefficient c
        For secondary matrix subtract from i-th row the j-th row with coefficient c
        '''
        for k in range(min_row, dim[0]):
            matrix[k][column_j] += coefficient*matrix[k][column_i]
        if secondary_matrix is not None:
            for k in range(secondary_dim[1]):
                secondary_matrix[column_i][k] -= coefficient*secondary_matrix[column_j][k]
            self._transform = self._helper.multiply_matrices(self._transform, [[self._add_transform_value(i, j, column_i, column_j, coefficient) for j in range(self._m1_base)] for i in range(self._m1_base)], (self._m1_base, self._m1_base), (self._m1_base, self._m1_base))

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
        m1_dim = self._helper.get_size(self._m1)
        m2_dim = self._helper.get_size(self._m2)
        # 1. Reduce M1. Done it in several steps. Each step zerofy one row and colum
        steps_count = min(m1_dim[0], m1_dim[1])
        for step in range(steps_count):
            # for each step consider submatrix of the M1 from (step, step) element
            # check if there are non-zero elements in reduction segment
            pos = self._get_minimal_coordinates(self._m1, m1_dim, step, step)
            should_do = pos != (-1, -1)
            while should_do:
                # there is non-zero element in the sumbatrix
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
                pos = self._get_minimal_coordinates(self._m1, m1_dim, step, step)
                should_do = pos == (-1, -1)
        # 2. Next reduce M2 without any effects on the M1
        # Find the number of non-zero elements in the diagonal of the M1
        shift = self._get_diagonal_length(self._m1, m1_dim)
        self._shift = shift
        # start from row shift + 1 and column 0
        steps_count = min(m2_dim[1], m2_dim[0] - shift)
        for step in range(steps_count):
            # repeat the same process
            while not self._is_row_column_zero(self._m2, m2_dim, step + shift, step):
                pos = self._get_minimal_coordinates(self._m2, m2_dim, step + shift, step)
                if pos[0] != step + shift:
                    self._switch_row(self._m2, m2_dim, step, step + shift, pos[0], save_to_transform=True)
                if pos[1] != step:
                    self._switch_columns(self._m2, m2_dim, step + shift, step, pos[1])
                if self._m2[step + shift][step] < 0:
                    self._negate_row(self._m2, m2_dim, step, step + shift, save_to_transform=True)
                a = self._m2[step + shift][step]
                for k in range(step + shift + 1, m2_dim[0]):
                    v = self._m2[k][step]
                    if v != 0:
                        c = abs(v) // a
                        self._add_row(self._m2, m2_dim, step, step + shift, k, c if v < 0 else -1*c, save_to_transform=True)
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
        # add to data_array remaining zero rows
        for i in range(shift + steps_count, m2_dim[0]):
            data_array.append(0)
        self._reducted = True
        self._data_diagonal = tuple(data_array)

    def __repr__(self):
        l = 0
        for r in self._m1:
            l = max(l, len(str(r)) - 2 - 2*(len(r) - 1))
        return "\n".join([self._helper.matrix_to_string(self._m1), "-"*(2*l + 1), self._helper.matrix_to_string(self._m2)])


class ChainComplex(object):
    def __init__(self, is_cochain=False):
        self._helper = MatrixHelper()
        self._is_cochain = is_cochain  # if True, then each map from n to n + 1, if False, then map from n no n - 1
        self._dimensions = {}  # key - group index, value - it dimension. If there is no key, then this group is zero
        self._boundaries = {}  # key - group index FROM, value - matrix (in the form of tupples) for this map
        self._factors = {}  # key - group index, value - array of generators of the factor subspace
        self._homology = {}  # key - group index, value - link to MatrixProcessor with reducted matrixes. If there is no key, the this homology froup are not calculated yet

    def generate_cochain_complex(self):
        '''Generate cochain complex which is dual to the gieven chain complex. If the current complex is cochain, then None is returned
        '''
        if self._is_cochain is True:
            print("Complex is cochain, nothing to generate.")
            return None
        else:
            return self._get_transpose_complex()

    def _get_transpose_complex(self):
        cc = ChainComplex(is_cochain=not self._is_cochain)  # create cochain complex
        # copy groups dimensions
        for g_index in self._dimensions.keys():
            cc.add_group(g_index, self._dimensions[g_index])
        # find maximal and minimal index
        indexes = list(self._dimensions.keys())
        indexes.sort()
        min_index = indexes[0]
        max_index = indexes[len(indexes) - 1]
        # copy boundaries, get it from previous group
        for index in range(min_index, max_index):
            if index + 1 in self._boundaries.keys():
                cc.add_boundary_map(index, self._helper.build_transpose(self._boundaries[index + 1]))
        # copy factors
        for f_index in self._factors.keys():
            cc.add_factor(f_index, [f for f in self._factors[f_index]])
        return cc

    def generate_chain_complex(self):
        '''Generate chain complex which is dual to the given cochain complex. If the current complex is chain, then None is returned
        '''
        if self._is_cochain is False:
            print("Complex is chain, nothing to generate")
            return None
        else:
            return self._get_transpose_complex()

    def add_group(self, index, dimension):
        '''Add group to the chain complex with specific index and specific dimension. If the group with this index already presented in the chain complex, then it will be overrided and computed homology will be droped.
        '''
        self._dimensions[index] = dimension
        if index in self._homology.keys():
            self._homology.pop(index)
        if index in self._boundaries.keys():
            self._boundaries.pop(index)
        if self._is_cochain and (index - 1) in self._boundaries.keys():
            self._boundaries.pop(index - 1)
        if self._is_cochain is False and (index + 1) in self._boundaries.keys():
            self._boundaries.pop(index + 1)

    def add_factor(self, index, factor_generators):
        '''Add to the chain group C_{index} factor-subgroup, generated by vectors in the list factor_generators

        List of factor_generators should be list of tuples of length dim (c_{index}). Each generator will be added as a column to the scond matrix for calculating homology group
        '''
        if index not in self._dimensions.keys():
            print("Set the dimension of the group C_" + str(index) + " at first.")
        else:
            # check is dimensions of factor generators are correct
            is_correct = True
            for gen in factor_generators:
                if len(gen) != self._dimensions[index]:
                    print("Generator " + str(gen) + " must be dimension " + str(self._dimensions[index]) + ".")
                    is_correct = False
            if is_correct:
                if index not in self._factors.keys():
                    self._factors[index] = []
                for gen in factor_generators:
                    self._factors[index].append(gen)
            else:
                print("Failed to add generators to the subgroup of the C_" + str(index))

    def _get_matrix_size(self, matrix):
        return (len(matrix), len(matrix[0]))

    def add_boundary_map(self, index, matrix):
        '''Add the boundary map to the chain complex as matrix.

        index is the index of the group FROM the map is acts. Matrix will be added if from- and to-groups are non-trivial.
        '''
        if index in self._boundaries.keys():
            self._boundaries.pop(index)
        if index in self._homology.keys():
            self._homology.pop(index)
        dims = self._get_matrix_size(matrix)
        if index in self._dimensions:
            from_dim = self._dimensions[index]
            to_dim = 0
            if self._is_cochain:
                if index + 1 in self._dimensions:
                    to_dim = self._dimensions[index + 1]
            else:
                if index - 1 in self._dimensions:
                    to_dim = self._dimensions[index - 1]
            if to_dim == 0:
                print("Does not need specific boundary map from the group C_" + str(index) + " to the trivial group C_" + str(index + 1 if self._is_cochain else index - 1))
            else:
                if (to_dim, from_dim) == dims:
                    self._boundaries[index] = tuple(tuple(matrix[i][j] for j in range(len(matrix[i]))) for i in range(len(matrix)))
                else:
                    print("Wrong dimension of the boundary map. It should acts from " + str(from_dim) + "-dim to " + str(to_dim) + "-dim.")
        else:
            print("Group C_" + str(index) + " is zero, does not need any specific boundary map.")

    def get_dimension(self, index):
        '''Return dimension of the group with index
        '''
        if index in self._dimensions.keys():
            return self._dimensions[index]
        else:
            return 0

    def get_boundary_map(self, index):
        if index in self._boundaries.keys():
            return self._boundaries[index]
        else:
            from_dim = self.get_dimension(index)
            to_dim = self.get_dimension(index + 1 if self._is_cochain else index - 1)
            return self._helper.build_constant_matrix((to_dim, from_dim), is_tuples=True)

    def get_boundary_map_kernel(self, index):
        '''Return an array of vectors, which form the basis of the boundary map kernel
        '''
        if index in self._boundaries.keys():
            helper = MatrixHelper()
            return helper.get_matrix_kernel(self._boundaries[index])
        else:
            return [tuple(1 if i == j else 0 for j in range(self.get_dimension(index))) for i in range(self.get_dimension(index))]

    def calculate_homology(self, index, force_recalculate=False):
        if (index in self._homology and self._homology[index].get_is_reducted() is False) or (index not in self._homology) or (force_recalculate is True):
            current_dim = self.get_dimension(index)
            to_dim = self.get_dimension(index + 1 if self._is_cochain else index - 1)
            prev_dim = self.get_dimension(index - 1 if self._is_cochain else index + 1)
            if current_dim > 0 and to_dim > 0 and prev_dim > 0:
                matrix_first = self._helper.copy_matrix(self.get_boundary_map(index))
                matrix_second = self._helper.copy_matrix(self.get_boundary_map(index - 1 if self._is_cochain else index + 1))
                if index in self._factors.keys():
                    if self._is_cochain:
                        # add factors as rows to the first matrix
                        matrix_first = self._helper.add_rows_to_matrix(matrix_first, self._factors[index])
                    else:
                        # add factors as columns for matrix_second
                        matrix_second = self._helper.add_columns_to_matrix(matrix_second, self._factors[index])
                mp = MatrixProcessor(matrix_first, matrix_second, log_error=False)
                self._homology[index] = mp
                if mp.get_is_correct():
                    mp.make_reduction()
                    return True
                else:
                    print(self._helper.get_error_string(mp.get_incorrect_code()))
                    return False
            else:
                if current_dim == 0:  # homology is trivial, because the group C_i = 0
                    mp = MatrixProcessor([[0]], [[1]])
                    self._homology[index] = mp
                    mp.make_reduction()
                    return True
                else:
                    if to_dim == 0 and prev_dim > 0:  # ->C->C->0->
                        matrix_second = self._helper.copy_matrix(self.get_boundary_map(index - 1 if self._is_cochain else index + 1))
                        matrix_first = self._helper.build_constant_matrix((1, len(matrix_second)), 0)
                        if index in self._factors.keys():
                            if self._is_cochain:
                                pass
                            else:
                                matrix_second = self._helper.add_columns_to_matrix(matrix_second, self._factors[index])
                        mp = MatrixProcessor(matrix_first, matrix_second)
                        self._homology[index] = mp
                        mp.make_reduction()
                        return True
                    elif to_dim > 0 and prev_dim == 0:  # ->0->C->C->
                        matrix_first = self.get_boundary_map(index)
                        matrix_second = self._helper.build_constant_matrix((len(matrix_first[0]), 1), 0)
                        if index in self._factors.keys():
                            if self._is_cochain:
                                pass
                            else:
                                matrix_second = self._helper.add_columns_to_matrix([[] for i in range(self._dimensions[index])], self._factors[index])
                        mp = MatrixProcessor(matrix_first, matrix_second)
                        self._homology[index] = mp
                        mp.make_reduction()
                        return True
                    else:  # ->0->C->0->
                        matrix_first = self._helper.build_constant_matrix((1, current_dim), 0)
                        matrix_second = self._helper.build_constant_matrix((current_dim, current_dim), 0)
                        if index in self._factors.keys():
                            if self._is_cochain:
                                pass
                            else:
                                matrix_second = self._helper.add_columns_to_matrix(matrix_second, self._factors[index])
                        mp = MatrixProcessor(matrix_first, matrix_second)
                        self._homology[index] = mp
                        mp.make_reduction()
                        return True
        else:
            return True

    def get_homology(self, index):
        '''Return homology group as tuple (a1, a2, ...)

        It means that H = Z_{a1} + Z_{a2} + ...
        Z_0 = Z
        Z_1 = {e} - trivial group
        '''
        is_success = False
        if (index not in self._homology.keys()) or (index in self._homology.keys() and self._homology[index].get_is_reducted() is False):
            # no calculations, try to do it
            is_success = self.calculate_homology(index)
        else:
            is_success = True
        if is_success:
            orders = self._homology[index].get_group_orders()
            if len(orders) == 0:
                return (1,)
            else:
                return self._homology[index].get_group_orders()
        else:
            print("Data for calculating homology group H_" + str(index) + " incorrect. Break calculations.")
            return None

    def get_generator(self, group_index, generator_index):
        '''Return the generator of the summand, which corresponds to generator_index, of the homology group H_{group_index}

        The result is a tuple in default basis of the chain group C_{group_index}, and None if homology is not calculated yet or generator_index is invalid
        '''
        is_success = False
        if (group_index not in self._homology.keys()) or (group_index in self._homology.keys() and self._homology[group_index].get_is_reducted() is False):
            # no calculations, try to do it
            is_success = self.calculate_homology(group_index)
        else:
            is_success = True
        if is_success:
            orders = self._homology[group_index].get_group_orders()
            if len(orders) > 0:
                shift = self._homology[group_index].get_nonkernel_count()
                if generator_index < len(orders):
                    return self._homology[group_index].get_transform_column(generator_index + shift)
                else:
                    print("Invalid generator index. Orders are " + str(orders) + " but you request generator with index " + str(generator_index) + ".")
            else:
                return tuple(0 for i in range(self._dimensions[group_index]))
        else:
            print("Data for calculating homology group H_" + str(group_index) + " incorrect. Break calculations.")
            return None

    def __repr__(self):
        str_array = []
        write_dim = False
        # get minimal and maximal index
        dim_keys = list(self._dimensions.keys())
        dim_keys.sort()
        min_index = dim_keys[0]
        max_index = dim_keys[len(dim_keys) - 1]
        # create iterator
        r = range(min_index, max_index + 1) if self._is_cochain is True else range(max_index, min_index - 1, -1)
        str_array.append("0")
        for i in r:
            if i in self._dimensions.keys():
                str_array.append("C_" + str(i) + (":dim=" + str(self._dimensions[i]) if write_dim else ""))
            else:
                str_array.append("0")
        str_array.append("0")
        return "->".join(str_array)


def example_basic():
    print("Example basic output:")
    helper = MatrixHelper()  # helper for matrix outputs
    # 1. Create chain complex
    cc = ChainComplex()
    # 2. Add chain groups by setting their indexes and dimensions
    cc.add_group(3, 1)
    cc.add_group(2, 3)
    cc.add_group(1, 4)
    cc.add_group(0, 2)
    # 3. Set boundary maps
    cc.add_boundary_map(3, [[0], [0], [0]])
    cc.add_boundary_map(2, [[1, 1, 1], [1, -1, -1], [-1, -1, 1], [-1, 1, -1]])
    cc.add_boundary_map(1, [[1, 1, 1, 1], [-1, -1, -1, -1]])
    # 4. Output scheme of the chain complex
    print(cc)
    # 5. For each group index
    for i in range(4):
        # 5.1. Calculate homology group
        orders = cc.get_homology(i)
        # 5.2. Output the isomorphism class of the group
        print("H_" + str(i) + " = " + helper.orders_to_string(orders))
        # 5.3. Output generators of summands in original chain group basis
        for v_index in range(len(orders)):
            print("\t order " + str(orders[v_index]) + ": " + str(cc.get_generator(i, v_index)))


def example_convert_complex():
    print("Example of conversion output:")
    helper = MatrixHelper()
    # 1. Create chain complex
    cc = ChainComplex()
    # consinder 0->C_2->C_1->0, where C_2 is 4-dim and C_1 is 3-dim
    # 2. Add groups
    cc.add_group(2, 4)
    cc.add_group(1, 3)
    # 3. Add boundary map from C_2 to C_1
    cc.add_boundary_map(2, [[1, -1, 2, 1], [2, -1, 3, -3], [1, -1, 1, -1]])
    # 4. Build corresponding cochain complex
    ccc = cc.generate_cochain_complex()
    # 5. calculate cohomology groups
    h1 = ccc.get_homology(1)
    h2 = ccc.get_homology(2)
    print("H^1 = " + helper.orders_to_string(h1))
    print("H^2 = " + helper.orders_to_string(h2))
    # 6. Show the generater cocycle for the H^2
    print(h2)
    for v_index in range(len(h2)):
        print("H^2 generator " + str(v_index) + ": " + str(ccc.get_generator(2, v_index)) + " has order " + str(h2[v_index]))


if __name__ == "__main__":
    example_basic()
    example_convert_complex()
