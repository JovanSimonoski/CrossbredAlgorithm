from crossbred import crossbred
from input_handler import parse_input
from mq import *
import argparse
import time

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Getting the parameters.")
    parser.add_argument("-n", type=int, help="The number of variables")
    parser.add_argument("-m", type=int, help="The number of equations")
    parser.add_argument("-k", type=int, help="The 'k' parameter")
    parser.add_argument("-D", type=int, help="The degree of the Macaulay matrix")
    args = parser.parse_args()

    n = int(args.n)
    m = int(args.m)
    k = int(args.k)
    degree = int(args.D)

    with open('current_time.txt', 'w') as f:
        f.write(str(time.time()))

    generate_mq(n, m)
    equations, answer = parse_input(n, m)
    crossbred(n, degree, k, equations, answer)
