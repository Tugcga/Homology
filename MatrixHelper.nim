import math
import bigints
import nimpy

proc get_size(matrix: seq[seq[int]]): (int, int) {.exportpy.} =
    return (len(matrix), len(matrix[0]))

proc matrix_to_string(matrix: seq[seq[int]]): string {.exportpy.} =
    result = ""
    for i, row in matrix:
        result.add("|")
        for j, v in matrix[i]:
            result.add($v)  # $... - to string
            if j < len(matrix[i]) - 1:
                result.add(" ")
        if i < len(matrix) - 1:
            result.add("|\n")
        else:
            result.add("|")
    return result

proc multiply_matrices(a: seq[seq[int]], b: seq[seq[int]]): seq[seq[int]] {.exportpy.} =
    var (a_r, a_c) = get_size(a)
    var (b_r, b_c) = get_size(b)
    if a_c == b_r:
        for i in 0..<a_r:
            var row: seq[int]
            for j in 0..<b_c:
                var value: int = 0
                for k in 0..<a_c:
                    value += a[i][k] * b[k][j]
                row.add(value)
            result.add(row)
    return result

proc build_constant_matrix(dimension: (int, int), value: int): seq[seq[int]] {.exportpy.} =
    for i in 0..<dimension[0]:
        var row: seq[int]
        for j in 0..<dimension[1]:
            row.add(value)
        result.add(row)
    return result

proc copy_matrix(matrix: seq[seq[int]]): seq[seq[int]] {.exportpy.} = 
    return matrix

proc add_columns_to_matrix(matrix: seq[seq[int]], columns: seq[seq[int]]): seq[seq[int]] {.exportpy.} = 
    for i in 0..<len(matrix):
        var row: seq[int]
        for j in 0..<len(matrix[i]):
            row.add(matrix[i][j])
        for j in 0..<len(columns):
            row.add(columns[j][i])
        result.add(row)
    return result

proc add_rows_to_matrix(matrix: seq[seq[int]], rows: seq[seq[int]]): seq[seq[int]] {.exportpy.} = 
    for i in 0..<len(matrix):
        var row: seq[int]
        for j in 0..<len(matrix[i]):
            row.add(matrix[i][j])
        result.add(row)
    for i in 0..<len(rows):
        var row: seq[int]
        for j in 0..<len(rows[i]):
            row.add(rows[i][j])
        result.add(row)
    return result

proc build_transpose(matrix: seq[seq[int]]): seq[seq[int]] {.exportpy.} =
    for i in 0..<len(matrix[0]):
        var row: seq[int]
        for j in 0..<len(matrix):
            row.add(matrix[j][i])
        result.add(row)
    return result

proc orders_to_string(orders: openArray[int]): string {.exportpy.} =
    if len(orders) > 0:
        for i, v in orders:
            if v == 0:
                result.add("Z")
            else:
                result.add("Z_" & $v)
            if i < len(orders) - 1:
                result.add(" + ")
    else:
        result = "0"
    return result

proc int_lcm(a: int, b: int): int {.exportpy.} = 
    return lcm(a, b)

proc int_gcd(a: int, b: int): int {.exportpy.} =
    return gcd(a, b)

proc lcm_n(values: openArray[int]): int {.exportpy.} =
    if len(values) == 0:
        return 0
    elif len(values) == 1:
        return values[0]
    else:
        result = lcm(values[0], values[1])
        for i in 2..<len(values):
            result = lcm(result, values[i])
    return result

proc get_matrix_kernel(matrix: seq[seq[int]]): seq[seq[int]] {.exportpy.} =
    return result

proc is_zero_row(row: openArray[int]): bool =
    for v in row:
        if v != 0:
            return false
    return true

proc is_zero_row(row: openArray[BigInt]): bool =
    for v in row:
        if v != 0:
            return false
    return true

proc is_zero_row(row: openArray[int64]): bool =
    for v in row:
        if v != 0:
            return false
    return true

proc abs_bigint(value: BigInt): BigInt =
    if value > 0:
        return value
    else:
        return initBigInt(0) - value

proc reduce(matrix: seq[seq[int]]): seq[seq[int64]] {.exportpy.} =
    # make copy of the matrix with values as BigInt
    # var m = copy_matrix(matrix)
    var m: seq[seq[int64]]
    for i in 0..<len(matrix):
        var row: seq[int64]
        for j in 0..<len(matrix[i]):
            row.add(int64(matrix[i][j]))
        m.add(row)

    var height: int = len(m)  # the number of rows in the matrix
    var width: int = len(m[0])  # the number of columns in the matrix
    var row_index: int = 0
    var column_index: int = 0
    while row_index < height and column_index < width:
        # find the minimum number (by absolute value) in the column column_index in the rows after row_index
        var c_min: int64
        var c_min_index: int = -1
        for i in row_index..<height:
            if m[i][column_index] != 0:
                var v: int64 = abs(m[i][column_index])
                if c_min_index == -1:
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
                for k in column_index..<width:
                    swap(m[c_min_index][k], m[row_index][k])
            if m[row_index][column_index] < 0:
                for k in column_index..<width:
                    m[row_index][k] = -1 * m[row_index][k]
            # decrease values in the column
            var a: int64 = abs(m[row_index][column_index])
            for i in 0..<height:
                # iterate by all rows
                if i != row_index:
                    var b: int64 = abs(m[i][column_index])
                    if b >= a:
                        # apply transform to the row i
                        if m[row_index][column_index] * m[i][column_index] > 0:
                            var c: int64 =  -1 * (b div a)
                            for j in column_index..<width:
                                m[i][j] = m[i][j] + c * m[row_index][j]
                        else:
                            var c: int64 = b div a
                            for j in column_index..<width:
                                m[i][j] = m[i][j] + c * m[row_index][j]
            # check, are all elements in the column now = 0
            var is_zeros: bool = true
            var is_finish: bool = false
            var i: int = row_index + 1
            while not is_finish:
                if i >= height:
                    is_finish = true
                else:
                    if m[i][column_index] != 0:
                        is_zeros = false
                        is_finish = true
                    i += 1
            if is_zeros:
                row_index += 1
                column_index += 1
    # remove zero rows
    for i in 0..<height:
        if not is_zero_row(m[i]):
            var row: seq[int64]
            for v in m[i]:
                row.add(v)
            result.add(row)
    #[for i in 0..<height:
        if not is_zero_row(m[i]):
            var row: seq[int]
            for v in m[i]:
                row.add(int(v.limbs[0]))  # <-- here we return not actual value, but it some part in BigInt representation
            result.add(row)]#
    return result


proc main(): void =
    var matrix = @[@[1, 2, 3, 4], @[2, 3, 4, 5], @[3, 4, 5, 6]]
    echo(get_size(matrix))

    echo(matrix_to_string(matrix))

    var a = @[@[1, 2, 3], @[2, 3, 4]]
    var b = @[@[1, 2], @[2, 1], @[1, 1]]

    var c = multiply_matrices(a, b)
    echo(matrix_to_string(c))

    var d = build_constant_matrix((2, 2), 1)
    echo(matrix_to_string(d))
    echo(copy_matrix(d))
    echo(add_columns_to_matrix(d, @[@[1, 2], @[2, 3]]))
    echo(add_rows_to_matrix(d, @[@[0, 0], @[1, 2]]))

    echo(build_transpose(a))

    echo(orders_to_string([2, 1, 0]))

    echo(int_lcm(6, 4), ", ", int_gcd(10, 12), ", ", lcm_n([2, 4, 6, 12]))

    var q = @[@[4, -9], @[-3, -10]]
    echo(reduce(q))

    var w = @[@[1, 2], @[3, 4]]
    var e = copy_matrix(w)
    e[0][0] = -1
    w[1][1] = -5
    echo(matrix_to_string(e))
    echo(matrix_to_string(w))
    var s = @[@[-3, -1, 6, 3, -6], @[-3, 3, 10, -3, 8], @[3, -6, -9, -3, 3], @[9, 10, -1, -4, -2], @[-10, -3, 7, 2, 10]]
    echo(reduce(s))

# main()
