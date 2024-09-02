import os
from random import randint


def generate_mq():
    """
    Generates a MQ polynomial systems with given number of variables and number of equations.
    Writes the generated system and solution into files with predefined file paths.

    Returns:
        n(int): Number of variables.
        m(int): Number of equations.
    """
    n = int(input('n: '))
    m = int(input('m: '))

    file_path = f'../example_input/auto_generated/example_n{n}_m{m}'

    if os.path.exists(file_path):
        return n, m

    solution = [randint(0, 1) for _ in range(n)]

    with open(f'{file_path}-answer', 'w') as f:
        f.write(str(solution) + '\n')

    with open(file_path, 'w') as f:
        f.write(
            f"Galois Field : GF(2)\n"
            f"Number of variables (n) : {n}\n"
            f"Number of polynomials (m) : {m}\n"
            f"Seed : 0\n"
            f"Order : "
            f"graded reverse lex order\n"
            f"\n*********************\n")

        result = 0
        for v in reversed(solution):
            result = result << 1 | v

        for i in range(m):
            result = 0
            for j in range(n):
                for k in range(j + 1):
                    r = randint(0, 1)
                    f.write(f'{r} ')

                    if r == 1:
                        result ^= solution[j] & solution[k]

            for j in range(n):
                r = randint(0, 1)
                f.write(f'{r} ')

                if r == 1:
                    result ^= solution[j]

            f.write(str(result) + ";\n")

    return n, m
