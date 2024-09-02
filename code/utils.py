import copy
import math

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
            print('\nGOT THE CORRECT SOLUTION:')
            print(s)
            print(answer)

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


def construct_first_sub_matrix(equations, monomials, num_variables, index_bin_dict, mon_bin_dict, bin_index_dict):
    """
    Constructs the first sub matrix which is a sub matrix of the original Macaulay matrix
    whose rows correspond to products u*f with deg_k(u) >= d-1.
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

    Returns:
        List[List[int]]: The first sub matrix.

    """
    len_monomials = len(monomials)
    mm_sub_matrix_1 = []
    max_degree_u = monomials[0].degree - 2

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
            mm_sub_matrix_1.append(copy.copy(line))

    mm_sub_matrix_1.extend(equations)
    return mm_sub_matrix_1


def construct_second_sub_matrix(k, mm_sub_matrix_1_sorted_deg_k, sorted_monomials_deg_k):
    """
    Extracts the second sub matrix from the first sub matrix.
    The second sub matrix contains the columns of the first sub matrix corresponding to monomials m with deg_k(m) > d.

    Parameters:
        k (int): The 'k' parameter
        mm_sub_matrix_1_sorted_deg_k (List[List[int]]): The first sub matrix represented as
            list of lists of integers (matrix) and sorted by deg_k.
        sorted_monomials_deg_k(List[Monomial]): List of all the monomials of max. degree D sorted by deg_k.

    Returns:
        List[List[int]]: The second sub matrix.
    """
    mm_sub_matrix_2_counter = 0
    for m in sorted_monomials_deg_k:
        if deg_k(m, k) <= 1:  # HARDCODED
            break
        mm_sub_matrix_2_counter += 1

    mm_sub_matrix_2 = []
    for m in mm_sub_matrix_1_sorted_deg_k:
        line = copy.copy(m[:mm_sub_matrix_2_counter])
        mm_sub_matrix_2.append(line)

    return mm_sub_matrix_2


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
