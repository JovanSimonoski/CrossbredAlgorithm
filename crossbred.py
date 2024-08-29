from sage.all import Matrix, GF, vector
from utils import *


def crossbred(num_variables, degree, k, equations, answer):
    monomials, monomials_degree_d, monomials_fukuoka_mq_challenge, sorted_monomials_deg_k = generate_monomials_types(
        degree, k, num_variables)

    len_monomials = len(monomials)
    len_monomials_degree_d = len(monomials_degree_d)

    """ --------- This code segment removes the x^2 monomials from the Fukuoka MQ Challenge input format --------- """
    format_equations_fukuoka(equations, monomials_fukuoka_mq_challenge)

    """ ---- Constructing the variables ---- """
    variables = generate_variables(num_variables)

    """ ---- Constructing dictionaries that will help for faster execution ---- """
    bin_index_dict, default_to_deg_k_index_dict, deg_k_to_default_index_dict, index_bin_dict, mon_bin_dict = (
        construct_dictionaries(monomials_degree_d, sorted_monomials_deg_k, num_variables))

    """ ---- Adding leading zero's to the equations to match the Macaulay matrix degree ---- """
    add_leading_zeros(equations, len_monomials_degree_d - len_monomials)

    """ ---- Constructing the first sub matrix Mac(F) sorted by deg_k and of degree D ---- """
    print('\nFirst sub matrix generating')
    mm_sub_matrix_1 = construct_first_sub_matrix(equations, monomials_degree_d, num_variables, index_bin_dict,
                                                 mon_bin_dict, bin_index_dict)
    print('DONE')

    print('\nSorting the sub matrix to deg_k')
    mm_sub_matrix_1_sorted_deg_k = sort_matrix_columns_by_dictionary(default_to_deg_k_index_dict, mm_sub_matrix_1)
    print('DONE')

    """ ---- Constructing the second sub matrix M(F) sorted by deg_k and of degree D ---- """
    print('\nSecond sub matrix generating')
    mm_sub_matrix_2 = construct_second_sub_matrix(k, mm_sub_matrix_1_sorted_deg_k, sorted_monomials_deg_k)
    print('DONE')

    """ ---- Finding vectors in the left kernel of the second sub matrix M(F) ----"""
    sage_sub_matrix_1 = Matrix(GF(2), mm_sub_matrix_1_sorted_deg_k)
    sage_sub_matrix_2 = Matrix(GF(2), mm_sub_matrix_2)
    print('\nFinding vectors in the left kernel of the second sub matrix M(F)')
    left_kernel = sage_sub_matrix_2.left_kernel()
    print('DONE')

    """ ---- Calculating the polynomials corresponding to vector_i * Mac(F) (first sub matrix) ----"""
    print(
        '\nCalculating the polynomials corresponding to vector_i * Mac(F) (first sub matrix) and sorting them to '
        'default')
    p_polynomials = []
    for v in left_kernel.basis():
        p_polynomials.append(list(v * sage_sub_matrix_1))

    p_polynomials_default_sorted = sort_matrix_columns_by_dictionary(deg_k_to_default_index_dict, p_polynomials)
    print('DONE')

    def fast_evaluate(system_polynomials, remaining, k_parameter, solution):
        if remaining == k_parameter:
            return solve_linear_system(k_parameter, solution, system_polynomials, num_variables, answer)

        x_l = variables[remaining - 1]
        x_l_binary = mon_bin_dict[x_l]

        x_l_appear_list = []
        for var in range(len_monomials_degree_d - 1):
            if x_l in monomials_degree_d[var].variables:
                x_l_appear_list.append(var)

        p_0_system = []
        p_1_system = []

        for poly in system_polynomials:
            p_0 = list(copy.copy(poly))
            p_1 = [0] * len_monomials_degree_d
            for var in x_l_appear_list:
                if poly[var] == 1:
                    p_0[var] = 0

                    current_mon = index_bin_dict[var]

                    if current_mon == x_l_binary:
                        new_mon_index = len_monomials_degree_d - 1
                    else:
                        res = ''
                        for bit in range(len(current_mon)):
                            res += str(int(current_mon[bit]) ^ int(x_l_binary[bit]))
                        new_mon_index = bin_index_dict[res]

                    p_1[new_mon_index] += 1
                    p_1[new_mon_index] %= 2

            p_1_result = []
            for p in range(len(p_1)):
                p_1_result.append(int(p_1[p]) ^ int(p_0[p]))

            p_0_system.append(list(p_0))
            p_1_system.append(p_1_result)

        s_1 = copy.copy(solution)
        s_1.append(1)

        s_0 = copy.copy(solution)
        s_0.append(0)

        fast_evaluate(p_0_system, remaining - 1, k_parameter, s_0)
        fast_evaluate(p_1_system, remaining - 1, k_parameter, s_1)

    print('\nCalling Fast Evaluate and waiting for solution')
    fast_evaluate(p_polynomials_default_sorted, num_variables, k, [])
