from crossbred import crossbred
from input_handler import parse_input
from mq import *

if __name__ == '__main__':
    n, m = generate_mq()
    degree, k, equations, answer = parse_input(n, m)
    crossbred(n, degree, k, equations, answer)
