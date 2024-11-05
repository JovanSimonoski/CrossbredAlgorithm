import copy
import math
import os
import time
from sage.all import Matrix, GF, vector
from functools import cmp_to_key
from itertools import combinations_with_replacement, combinations

variables_str = []


class Monomial:
    """
    A class representing a monomial.
    Each monomial in constructed of list of variables(strings) and has a degree equal to the number of variables.
    """

    def __init__(self, variables):
        if not variables:
            self.variables = ['1']
        else:
            self.variables = copy.deepcopy(list(variables))

        if len(variables) > 1:
            sort_monomial_variables(self.variables)

        self.degree = len(variables)

    def __get_as_tuple__(self):
        return tuple(self.variables)

    def __str__(self):
        return ''.join(variable for variable in self.variables)

    def __eq__(self, other):
        if str(self) == str(other):
            return True
        return False

    def __hash__(self):
        return hash((tuple(self.variables), self.degree))


def generate_monomials_types(degree, k, num_variables):
    """
    Generates all the necessary lists of different monomial types with 'num_variables' variables
    and max. degree 'degree'.

        -monomials[list of Monomial objects with max. degree of 2, sorted in grevlex]

        -monomials_degree_d[list of Monomial objects with max. degree of d, sorted in grevlex]

        -monomials_fukuoka_mq_challenge[list of Monomial objects with max. degree of d, sorted in grevlex,
            based on the Fukuoka MQ Challenge format]

        -sorted_monomials_deg_k[monomials_degree_d sorted by deg_k]

    Parameters:
        degree(int): The max. degree of the monomials.
        k(int): The 'k' parameter
        num_variables(int): The number of variables.

    Returns:
        List[Monomial]: The generated monomials
    """

    monomials = generate_monomials(num_variables, 2)
    monomials_degree_d = generate_monomials(num_variables, degree)
    monomials_fukuoka_mq_challenge = generate_monomials_fukuoka_format(num_variables, 2)
    sorted_monomials_deg_k = sort_deg_k_grevlex(copy.deepcopy(monomials_degree_d), k)
    return monomials, monomials_degree_d, monomials_fukuoka_mq_challenge, sorted_monomials_deg_k


def generate_monomials(num_variables, degree):
    """
    Generates all the monomials with 'num_variables' variables and max. degree 'degree',
    having that x^2 = x based on the GF(2). The monomials are sorted in grevlex.

    Parameters:
        num_variables(int): The number of variables.
        degree(int): The max. degree of the monomials.

    Returns:
        List[Monomial]: The generated monomials
    """

    var = []
    for i in range(1, num_variables + 1):
        var.append('x' + str(i))

    global variables_str
    variables_str = copy.deepcopy(var)
    variables_str.append('1')

    monomials_tuples = []

    for d in range(degree + 1):
        for combination in combinations(var, d):
            monomials_tuples.append(combination)

    sort_grevlex(monomials_tuples)

    monomials_list = []
    for mon in monomials_tuples:
        m_tmp_list = []
        for variable in mon:
            m_tmp_list.append(variable)
        monomials_list.append(m_tmp_list)

    monomials_obj = []
    for mon in monomials_list:
        monomials_obj.append(Monomial(mon))

    return monomials_obj


def generate_monomials_fukuoka_format(num_variables, degree):
    """
    Generates all the monomials with 'num_variables' variables and max. degree 'degree'
    based on the Fukuoka MQ Challenge format. The monomials are sorted in grevlex.

    Parameters:
        num_variables(int): The number of variables.
        degree(int): The max. degree of the monomials.

    Returns:
        List[Monomial]: The generated monomials
    """

    var = []
    for i in range(1, num_variables + 1):
        var.append('x' + str(i))

    monomials_tuples = []

    for degree in range(degree + 1):
        for comb in combinations_with_replacement(var, degree):
            monomials_tuples.append(comb)

    sort_grevlex(monomials_tuples)

    monomials_list = []
    for mon in monomials_tuples:
        m_tmp_list = []
        for variable in mon:
            m_tmp_list.append(variable)
        monomials_list.append(m_tmp_list)

    monomials_new = []
    # Only works for degree 2
    for monomial in monomials_list:
        if len(monomial) > 1 and monomial[0] == monomial[1]:
            monomials_new.append([monomial[0]])
        else:
            monomials_new.append(monomial)

    monomials_obj = []
    for mon in monomials_new:
        monomials_obj.append(Monomial(mon))

    return monomials_obj


def sort_deg_k_grevlex(monomials, k):
    """
    Sorts monomials by deg_k.

    Parameters:
        monomials(List[Monomial]): The monomials we want to sort.
        k(int): The 'k' parameter.

    Returns:
        List[Monomial]: The sorted monomials
    """

    grevlex = cmp_to_key(sort_monomial_grevlex)

    one = False
    if Monomial([]) in monomials:
        one = True
        monomials.remove(Monomial([]))

    sorted_monomials_deg_k = copy.deepcopy(monomials)
    sorted_monomials_deg_k = copy.deepcopy(sorted(sorted_monomials_deg_k,
                                                  key=lambda monomial: (
                                                      -deg_k(monomial, k), grevlex(monomial.__get_as_tuple__()))
                                                  ))

    if one:
        sorted_monomials_deg_k.append(Monomial([]))

    return sorted_monomials_deg_k


def sort_monomial_variables(variables):
    """
    Sorts the variables of a monomial in grevlex.

    Parameters:
        variables(List[string]): The list of variables
    """

    variables.sort(key=lambda variable: int(variable[1:]))


def sort_grevlex(list_to_sort):
    """
    Sorts a list in grevlex.

    Parameters:
        list_to_sort(List[tuple(string)]): The list of monomials(represented as tuples of strings-variables).
    """

    one = False
    if list_to_sort[0] == tuple(''):
        one = True
        list_to_sort.remove(list_to_sort[0])

    list_to_sort.sort(key=cmp_to_key(sort_monomial_grevlex))

    if one:
        list_to_sort.append(tuple(''))


def sort_monomial_grevlex(monomial_1, monomial_2):
    """
    Compares two monomials.

    Parameters:
        monomial_1(tuple(string)): The first monomial. Each variable is a string in a tuple.
        monomial_2(tuple(string)): The second monomial. Each variable is a string in a tuple.

    Returns:
        int: 1 if monomial_1 > monomial_2, -1 if monomial_1 < monomial_2, 0 otherwise
    """

    len_mon_1 = len(monomial_1)
    len_mon_2 = len(monomial_2)

    if len_mon_1 != len_mon_2:
        return 1 if len_mon_1 < len_mon_2 else -1

    for i in range(len_mon_1 - 1, -1, -1):
        deg1 = int(monomial_1[i][1:])
        deg2 = int(monomial_2[i][1:])
        if deg1 != deg2:
            return 1 if deg1 > deg2 else -1
    return 0


def deg_k(monomial, k):
    """
    Calculates the deg_k degree of a monomial

    Parameters:
        monomial(Monomial): The monomial we want to get the degree of
        k(int): The 'k' parameter

    Returns:
        int: The deg_k degree of the monomial
    """

    degree = 0

    for v in monomial.variables:
        if v == '1':
            return degree
        if int(v[1:]) <= k:
            degree += 1
        else:
            return degree

    return degree


def add_leading_zeros(equations, len_monomials):
    """
    Adds leading zeros to the equations depending on the desired length.

    Parameters:
        equations(List[List[int]]): The list of equations. Each equation is represented as a list of integers.
        len_monomials(int): The desired length of the equations
    """

    tmp_zeros = [0] * len_monomials
    for e in equations:
        e.reverse()
        e.extend(tmp_zeros)
        e.reverse()


def generate_variables(num_variables):
    """
    Generates all variables corresponding to the number of variables.

    Parameters:
        num_variables(int): The number of variables.

    Returns:
        List[string]: List of variables
    """

    variables = []
    for i in range(num_variables):
        variables.append(f'x{i + 1}')
    return variables


def format_equations_fukuoka(equations, monomials_fukuoka_mq_challenge):
    """
    Formats the given equations which are in Fukuoka MQ Challenge format, that means including monomials like: x^2.
    This function removes those columns and recalculates the equation if the values of those columns were '1'.

    Parameters:
        equations(List[List[int]]): The equations we are formatting. Each equation is represented as a list of integers.
        monomials_fukuoka_mq_challenge(List[Monomial]): The monomials corresponding to the Fukuoka MQ Challenge format
    """

    len_monomials_fukuoka_challenge = len(monomials_fukuoka_mq_challenge)
    positions_dict = {}
    indexes_to_remove = []
    for i in range(len_monomials_fukuoka_challenge - 2, - 1, -1):
        if len(monomials_fukuoka_mq_challenge[i].variables) == 1:
            mon = str(monomials_fukuoka_mq_challenge[i])
            if mon not in positions_dict:
                positions_dict.update({f'{mon}': i})
            else:
                positions_dict.update({f'{mon}': (positions_dict[mon], i)})
                indexes_to_remove.append(i)
    for e in equations:
        for key in positions_dict.keys():
            index_to_remove = positions_dict[key][1]
            index_to_sustain = positions_dict[key][0]
            if e[index_to_remove] == 1:
                e[index_to_sustain] += 1
                e[index_to_sustain] %= 2
        for index in indexes_to_remove:
            del e[index]


def check_consistency(system_polynomials):
    """
    Checks if a given system of polynomials is consistent or not.

    Parameters:
        system_polynomials(List[List[int]]): The system we're trying to check

    Returns:
        Boolean: True if the system is consistent, otherwise False
    """

    for s in system_polynomials:
        if all(x == 0 for x in s[:-1]) and s[-1] == 1:
            return False
    return True


def solve_linear_system(k, solution, system_polynomials, num_variables, answer):
    """
    Tries to solve a given linear system.

    Parameters:
        k(int): The 'k' parameter.
        solution(List[int]): The current solution (values of the previously fixed variables).
        system_polynomials(List[List[int]]): The system we're trying to solve
        num_variables(int): The number of variables.
        answer(List[int]): The actual answer to the example,
            derived from the Fukuoka MQ Challenge input (the answer file).

    Returns:
        None: If no solution is found or the system is checked as inconsistent.
        If solution is found, it doesn't return anything. It prints the solution
            and exits the program with exit code '0'.

    Raises:
        ValueError: If system is unsolvable.
    """

    if not check_consistency(system_polynomials):
        return

    linear_system = []
    for s in system_polynomials:
        linear_system.append(copy.copy(s[len(s) - num_variables - 1: len(s) - num_variables + k - 1]))

    sage_linear_system = Matrix(GF(2), linear_system)

    constant_column = []
    for s in system_polynomials:
        constant_column.append(s[-1])
    sage_vector = vector(GF(2), constant_column)

    try:
        kernel = sage_linear_system.solve_right(sage_vector)
        s = []
        for value in kernel:
            s.append(value)
        s.extend(solution[::-1])

        if all(val == 0 for val in s):
            return
        if s == answer:
            end_time = time.time()
            path = f'../runs_data/run_n{num_variables}_k{k}'

            with open(f'{path}/starting_time.txt', 'r') as f:
                start_time = float(f.read())
            if os.path.exists(f'{path}/starting_time.txt'):
                os.remove(f'{path}/starting_time.txt')

            with open(f'{path}/info.txt', 'a') as f:
                f.write('\n' + str(s) + '\n')
                f.write(str(answer) + '\n')
                f.write(f'\nFinished in {round((end_time - start_time) / 60, 4)} minutes.\n')
            exit(0)
    except ValueError:
        pass

    return


def construct_dictionaries(monomials, sorted_monomials_deg_k, num_variables):
    """
    Constructs dictionaries that will later be used for different purposes mainly to
    speed up the process of the matrix construction and sorting.

        -mon_bin_dict{monomial(string): binary representation of the monomial(List[int])}:

        -index_bin_dict{index(int): binary representation of the monomial corresponding
            to that index in deg_k sorted equations(List[int])}

        -bin_index_dict{the reverse of index_bin_dict}

        -default_to_deg_k_index_dict{index in the grevlex sorted equation(int): corresponding
            index in the deg_k sorted equations(int)}

        -deg_k_to_default_index_dict{the reverse of default_to_deg_k_index_dict}

    Parameters:
        monomials(List(Monomial)): List of all the monomials with max. degree 2 sorted in grevlex.
        sorted_monomials_deg_k(List[Monomial]): List of all the monomials of max. degree D sorted by deg_k.
        num_variables(int): Number of variables.

    Returns:
        All the constructed dictionaries
    """

    len_monomials = len(monomials)
    mon_bin_dict = {}
    for mon in monomials[:-1]:
        bin_rep = [0] * num_variables
        for var in mon.variables:
            bin_rep[int(var[1:]) - 1] = 1
        bin_rep = ''.join(str(bit) for bit in bin_rep)
        mon_bin_dict.update({str(mon): bin_rep})

    mon_index_dict = {}
    index_bin_dict = {}
    bin_index_dict = {}
    for i in range(len_monomials - 1):
        mon_index_dict.update({str(monomials[i]): i})
        index_bin_dict.update({i: mon_bin_dict[str(monomials[i])]})
        bin_index_dict.update({mon_bin_dict[str(monomials[i])]: i})
    mon_index_dict.update({str(monomials[-1]): len_monomials - 1})

    default_to_deg_k_index_dict = {}
    deg_k_to_default_index_dict = {}
    for i in range(len_monomials):
        mon_deg_k = str(sorted_monomials_deg_k[i])
        default_index = mon_index_dict[mon_deg_k]
        if default_index == i:
            continue
        default_to_deg_k_index_dict.update({default_index: i})
        deg_k_to_default_index_dict.update({i: default_index})

    return bin_index_dict, default_to_deg_k_index_dict, deg_k_to_default_index_dict, index_bin_dict, mon_bin_dict


def sort_matrix_columns_by_dictionary(indexes_dictionary, matrix):
    """
    Sorts a given matrix by the key:value pairs it the given indexes dictionary.

    Parameters:
        indexes_dictionary(dictionary): Dictionary of keys of indexes and values of corresponding indexes
            by whom the positions of the values are switched and the matrix is being sorted.
        matrix(List[List[int]]): List of lists of integer, the matrix we want to sort.

    Returns:
        List[List[int]]: The sorted matrix.
    """

    sorted_matrix = []
    for m in matrix:
        m_sorted = [0] * len(m)
        for i in range(len(m)):
            if i in indexes_dictionary.keys():
                new_index = indexes_dictionary[i]
                m_sorted[new_index] = m[i]
            else:
                m_sorted[i] = m[i]
        sorted_matrix.append(m_sorted)

    return sorted_matrix


def construct_first_sub_matrix_and_save_to_file(equations, monomials, num_variables, index_bin_dict, mon_bin_dict,
                                                bin_index_dict, indexes_dictionary, sorted_monomials_deg_k, k):
    """
    Constructs the first sub matrix which is a sub matrix of the original Macaulay matrix
    whose rows correspond to products u*f with deg_k(u) >= d-1. Then, save it to a file instead of keeping in
    operational memory,
    (d is hard coded as 1 for the entire implementation which means that the first sub matrix is the actual original
    Macaulay matrix but for the sake of efficiency we don't store both)

    Parameters:
        equations(List[List[int]]): List of lists of integers. Each equation is represented as a list of integers.
        monomials(List(Monomial)): List of all the monomials with max. degree 2 sorted in grevlex.
        num_variables(int): Number of variables.
        index_bin_dict(dictionary): Dictionary with keys of indexes and values of binary representation
            of the corresponding monomial (that is represented by the index/position in the equation/column).
        mon_bin_dict(dictionary): Dictionary with keys of monomials(string) and values of the binary representation
            of the monomial.
        bin_index_dict(dictionary): Dictionary with keys of binary representations of monomials and values of index
            of the corresponding monomial in an equation.
        indexes_dictionary(dictionary): Dictionary of keys of indexes and values of corresponding indexes
            by whom the positions of the values are switched and the matrix is being sorted.
        sorted_monomials_deg_k(List[Monomial]): List of all the monomials of max. degree D sorted by deg_k.
        k (int): The 'k' parameter

    Returns:
        List[List[int]]: The first sub matrix.

    """

    path = f'../runs_data/run_n{num_variables}_k{k}'
    len_monomials = len(monomials)
    max_degree_u = monomials[0].degree - 2

    tmp_counter = 0
    for i in range(len(sorted_monomials_deg_k) - 1, -1, -1):
        if deg_k(sorted_monomials_deg_k[i], k) > 1:
            break
        tmp_counter += 1

    linear_separator = len(sorted_monomials_deg_k) - tmp_counter

    with open(f'{path}/macaulay_matrix_1.txt', 'a') as sub_matrix_1_file, open(f'{path}/macaulay_matrix_2.txt',
                                                                               'a') as sub_matrix_2_file:
        for mon in monomials[:-1]:
            if mon.degree > max_degree_u:
                continue
            mon_bin = mon_bin_dict[str(mon)]
            for e in equations:
                line = [0] * len_monomials
                for i in range(len_monomials - 1):
                    if e[i] == 1:
                        result = ''
                        current_mon_binary = index_bin_dict[i]
                        for j in range(num_variables):
                            result += str(int(mon_bin[j]) | int(current_mon_binary[j]))
                        index = bin_index_dict[result]
                        line[index] += 1
                        line[index] %= 2
                if e[-1] == 1:
                    index = bin_index_dict[mon_bin]
                    line[index] += 1
                    line[index] %= 2

                line_sorted = [0] * len(line)
                for i in range(len(line)):
                    if i in indexes_dictionary.keys():
                        new_index = indexes_dictionary[i]
                        line_sorted[new_index] = line[i]
                    else:
                        line_sorted[i] = line[i]

                sub_matrix_2_file.write(str(line_sorted[:linear_separator]) + '\n')
                sub_matrix_1_file.write(str(line_sorted[linear_separator:]) + '\n')

        for e in equations:
            line_sorted = [0] * len(e)
            for i in range(len(e)):
                if i in indexes_dictionary.keys():
                    new_index = indexes_dictionary[i]
                    line_sorted[new_index] = e[i]
                else:
                    line_sorted[i] = e[i]

            sub_matrix_2_file.write(str(line_sorted[:linear_separator]) + '\n')
            sub_matrix_1_file.write(str(line_sorted[linear_separator:]) + '\n')


def get_size_of_mm(n, degree, m):
    """
    Calculates the size of a Macaulay matrix for given number of values n, degree D and number of equation m.

    Parameters:
        n(int): Number of variables.
        degree(int): Degree of the Macaulay matrix.
        m(int): Number of equations.
    """

    num_cols = 0
    for i in range(degree, -1, -1):
        num_cols += math.comb(n, i)

    num_rows = (n + 1) * m
    for i in range(2, degree - 1):
        print(f'i={i}')
        num_rows += math.comb(n, i) * m

    print(f'-----------\n'
          f'Number of Rows x Columns:{num_rows}x{num_cols}'
          f'\n-----------')


def rref(matrix):
    """
    Function that transforms a binary matrix into Reduced Row Echelon Form

    Parameters:
        matrix(List[List[int]]): The matrix that we want to transform
    """

    rows, cols = len(matrix), len(matrix[0])
    lead = 0
    for r in range(rows):
        if lead >= cols:
            return matrix
        i = r
        while matrix[i][lead] == 0:
            i += 1
            if i == rows:
                i = r
                lead += 1
                if cols == lead:
                    return matrix
        matrix[i], matrix[r] = matrix[r], matrix[i]

        for i in range(rows):
            if i != r and matrix[i][lead] == 1:
                matrix[i] = [(iv ^ rv) for iv, rv in zip(matrix[i], matrix[r])]
        lead += 1


def read_sub_matrix(sub_matrix, num, path):
    """
    Function that reads a sub matrix from a file.

    Parameters:
        sub_matrix(): The sub matrix
        num(int): The number of the sub matrix
        path(String): Path to the current working directory
    """

    with open(f'{path}/sub_mm_{num}.txt', 'r') as file:
        for line in file:
            line = line.strip().strip('[]')
            line = line.strip().strip(',')
            row = [int(x) for x in line.split(', ')]
            sub_matrix.append(row)


def write_sub_matrix(sub_matrix, num, path):
    """
    Function that writes a sub matrix to a file.

    Parameters:
        sub_matrix(): The sub matrix
        num(int): The number of the sub matrix
        path(String): Path to the current working directory
    """

    with open(f'{path}/sub_mm_{num}.txt', 'w') as file:
        for row in sub_matrix:
            file.write(str(row) + '\n')


def transpose_matrix_in_iterations(matrix_file, path):
    """
    Function that applies one matrix to another one using Gaussian Elimination

    Parameters:
        matrix_file(String): Name of the file where the matrix is stored
        path(String): Path to the current working directory
    """

    temp_files = []
    with open(f'{path}/{matrix_file}', 'r') as f:
        for row_index, line in enumerate(f):
            line = line.strip().strip('[]')
            row = [int(x) for x in line.split(',')]

            temp_files.append(f'{path}/tmp_{row_index}')

            with open(f'{path}/tmp_{row_index}', 'w') as tmp_file:
                for element in row:
                    tmp_file.write(str(element) + ', \n')

    with open(f'{path}/transposed_mm.txt', 'w') as out_f:
        temp_files = [open(f'{path}/tmp_{i}', 'r') for i in range(len(temp_files))]
        for col_line in zip(*temp_files):
            out_f.write(' '.join(line.strip() for line in col_line) + '\n')

    for file in temp_files:
        file.close()
        os.remove(file.name)


def apply_matrix(matrix_to_apply, matrix_to_be_applied_on):
    """
    Function that applies one matrix to another one using Gaussian Elimination

    Parameters:
        matrix_to_apply(List[List[int]]): The matrix we are applying
        matrix_to_be_applied_on(List[List[int]]): The matrix that we want to apply the changes to
    """

    for row in matrix_to_apply:
        pivot = -1
        for p in range(len(row)):
            if row[p] == 1:
                pivot = p
                break
        if pivot != -1:
            for i in range(len(matrix_to_be_applied_on)):
                if matrix_to_be_applied_on[i][pivot] == 1:
                    matrix_to_be_applied_on[i] = [row[q] ^ matrix_to_be_applied_on[i][q] for q in range(len(row))]


def find_and_swap_missing_pivot(reduced_sub_matrix, other_sub_matrix, starting_column=0):
    """
    Finds and swaps missing pivot rows from one sub matrix to another

    Parameters:
        reduced_sub_matrix(List[List[int]]): The main sub matrix who has missing pivot rows
        other_sub_matrix(List[List[int]]): The other sub matrix from whom
            we get the rows we need in the main sub matrix
        starting_column(int): Starting column of the sub matrix considering it's not the whole matrix
    """

    rows_reduced = len(reduced_sub_matrix)
    columns_reduced = len(reduced_sub_matrix[0])
    rows_other = len(other_sub_matrix)

    r = 0
    for col in range(starting_column, columns_reduced, 1):
        if r == rows_reduced:
            break
        if reduced_sub_matrix[r][col] == 0:
            for row in range(rows_other):
                if other_sub_matrix[row][col] == 1:
                    tmp = other_sub_matrix[row]
                    other_sub_matrix[row] = reduced_sub_matrix[r]
                    reduced_sub_matrix[r] = tmp
                    break
        r += 1


def check_rref(matrix, starting_column=0):
    """
    Check if a matrix is in a Reduced Row Echelon Form
    * It is used only for the first 2 sub matrices (may not be always correct)

    Parameters:
        matrix(List[List[int]]): The sub matrix
        starting_column(int): Starting column of the sub matrix considering it's not the whole matrix

    Returns:
        bool: True if the matrix is in Reduced Row Echelon Form, false otherwise
    """

    rows = len(matrix)
    columns = len(matrix[0])

    r = 0
    for i in range(starting_column, columns):
        if matrix[r][i] == 0:
            return False
        r += 1
        if r >= rows:
            break

    return True


def get_pivots(sub_matrix, pivots, starting_column):
    """
    Extracts the indexes of the pivots from a sub matrix and adds them to the pivots list

    Parameters:
        sub_matrix(List[List[int]]): The sub matrix
        pivots(List[int]): Indexes of the pivots
        starting_column(int): Starting column of the sub matrix considering it's not the whole matrix
    """

    rows, cols = len(sub_matrix), len(sub_matrix[0])

    for r in range(rows):
        for c in range(starting_column, cols):
            if sub_matrix[r][c] == 1:
                pivots.append(c)
                break


def construct_null_space_basis(null_space_basis, free_vars, pivots, sub_matrix, starting_row=0):
    """
    Constructs the null space basis of a matrix in many iterations using sub matrices

    Parameters:
        null_space_basis(List[List[int]]): The null space basis
        free_vars(List[int]): Indexes of free variables
        pivots(List[int]): Indexes of the pivots
        sub_matrix(List[List[int]]): The sub matrix
        starting_row(int): Starting row of the sub matrix considering it's not the whole matrix
    """

    for index, free_var in enumerate(free_vars):
        null_vector = null_space_basis[index]
        for r in range(len(sub_matrix)):
            if sub_matrix[r][free_var] == 1:
                null_vector[pivots[r + starting_row]] = 1
