import copy
from sage.all import Matrix, GF, vector
from functools import cmp_to_key
from itertools import combinations_with_replacement, combinations

variables_str = []


class Monomial:
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
    monomials = generate_monomials(num_variables, 2)
    monomials_degree_d = generate_monomials(num_variables, degree)
    monomials_fukuoka_mq_challenge = generate_monomials_fukuoka_format(num_variables, 2)
    sorted_monomials_deg_k = sort_deg_k_grevlex(copy.deepcopy(monomials_degree_d), k)
    return monomials, monomials_degree_d, monomials_fukuoka_mq_challenge, sorted_monomials_deg_k


def generate_monomials(num_var, degree):
    var = []
    for i in range(1, num_var + 1):
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


def generate_monomials_fukuoka_format(num_var, degree):
    var = []
    for i in range(1, num_var + 1):
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
    variables.sort(key=lambda variable: int(variable[1:]))


def sort_grevlex(list_to_sort):
    one = False
    if list_to_sort[0] == tuple(''):
        one = True
        list_to_sort.remove(list_to_sort[0])

    list_to_sort.sort(key=cmp_to_key(sort_monomial_grevlex))

    if one:
        list_to_sort.append(tuple(''))


def sort_monomial_grevlex(monomial_1, monomial_2):
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
    tmp_zeros = [0] * len_monomials
    for e in equations:
        e.reverse()
        e.extend(tmp_zeros)
        e.reverse()


def generate_variables(num_variables):
    variables = []
    for i in range(num_variables):
        variables.append(f'x{i + 1}')
    return variables


def format_equations_fukuoka(equations, monomials_fukuoka_mq_challenge):
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
    for s in system_polynomials:
        if all(x == 0 for x in s[:-1]) and s[-1] == 1:
            return False
    return True


def solve_linear_system(k_parameter, solution, system_polynomials, num_variables, answer):
    if not check_consistency(system_polynomials):
        return

    linear_system = []
    for s in system_polynomials:
        linear_system.append(copy.copy(s[len(s) - num_variables - 1: len(s) - num_variables + k_parameter - 1]))

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
            print('ANSWER:')
            print(answer)

            exit(0)
    except ValueError:
        pass

    return


def construct_dictionaries(monomials, sorted_monomials_deg_k, num_variables):
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

    return mm_sub_matrix_1


def construct_second_sub_matrix(k, mm_sub_matrix_1_sorted_deg_k, sorted_monomials_deg_k):
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
