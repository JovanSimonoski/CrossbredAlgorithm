# Crossbred Algorithm Implementation

This is a Python implementation of the crossbred algorithm presented in the paper [_"A crossbred algorithm for solving Boolean polynomial systems"_](https://ia.cr/2017/372) by Antoine Joux and Vanessa Vitse.

## Algorithm Overview

The general steps of the algorithm are given in this pseudo code from the paper:

![Algorithm Pseudo Code](https://github.com/JovanSimonoski/CrossbredAlgorithm/blob/master/pseudo_code.png?raw=true)

## How to Run It
Clone the repo:
```bash
git clone https://github.com/JovanSimonoski/CrossbredAlgorithm.git
cd CrossbredAlgorithm
cd code
```
To execute the algorithm, use the following example command:
```bash
sage main.py -n 10 -m 20 -k 5 -D 2
```
The arguments for the program are: 
  - the number of variables 'n', 
  - the number of equations 'm',
  - the 'k' parameter and
  - the degree of the Macaulay matrix 'D'. 

A random system will be generated for the 'n' and 'm' parameters. 

The program then tries to solve the system and compares the obtained solution with the correct solution of the generated system.

This implementation works only for hard coded value of d = 1.

