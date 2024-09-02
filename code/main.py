from crossbred import crossbred
from input_handler import parse_input

if __name__ == '__main__':
    num_variables, degree, k, equations, answer = parse_input()
    crossbred(num_variables, degree, k, equations, answer)
