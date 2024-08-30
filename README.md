# Crossbred Algorithm Implementation

This is a Python implementation of the crossbred algorithm presented in the paper _"A crossbred algorithm for solving Boolean polynomial systems"_ by Antoine Joux and Vanessa Vitse.

## Algorithm Overview

The general steps of the algorithm are given in this pseudo code from the paper:

![Algorithm Pseudo Code](https://github.com/JovanSimonoski/CrossbredAlgorithm/blob/master/pseudo_code.png?raw=true)

## How to Run It
Clone the repo:
```bash
git clone https://github.com/JovanSimonoski/CrossbredAlgorithm.git
cd CrossbredAlgorithm
```
To execute the algorithm, use the following command:
```bash
sage main.py
```
Next you will be promted to choose: 
  - the number of variables, 
  - the seed,
  - the 'k' parameter and
  - the degree of the Macaulay matrix - 'D'. 

This implementation works only for hard coded value of d = 1.

