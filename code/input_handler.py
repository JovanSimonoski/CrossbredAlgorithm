def parse_input(n, m):
    """
    Parsing the input based on the Fukuoka MQ Challenge format and then getting few parameters from the user:
        n(number of variables),
        seed(the seed of the example we want to try - Fukuoka related),
        k(the 'k' parameter that represents how many variables to keep / number of not fixed variables),
        degree(the degree of the Macaulay matrix, represented in the paper as 'D').

    Parameters:
        n(int): Number of variables.
        m(int): Number of equations.

    Returns:
        degree(int): Degree.
        k(int): The 'k' parameter.
        eq(List[List[int]]): Equations from the input.
        answer(List[int]): The answer of the example.
    """
    file_path = f'../example_input/auto_generated/example_n{n}_m{m}'

    input_list = []
    with open(file_path, 'r') as file:
        for line in file:
            input_list.append(line.strip())

    input_separator_index = input_list.index('*********************')

    input_values = input_list[:input_separator_index - 1]
    input_equations = input_list[input_separator_index + 1:]

    values = []
    for line in input_values:
        values.append(line.split(':')[-1].strip())

    k = input('k: ')
    degree = input('D: ')

    values.append(degree)
    values.append(k)

    variables = ['galois_field', 'num_variables', 'num_equations', 'seed', 'order', 'degree', 'k']
    var_val = dict(zip(variables, values))
    eq = [list(map(int, line.strip(' ;').split())) for line in input_equations]

    degree = int(var_val['degree'])
    k = int(var_val['k'])

    answer_tmp = []
    with open(f'{file_path}-answer', 'r') as file:
        for line in file:
            answer_tmp = line.strip()[1:-1]

    answer_tmp = answer_tmp.split(', ')

    answer = []
    for a in answer_tmp:
        answer.append(int(a))

    return degree, k, eq, answer
