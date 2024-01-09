NOMAD is a blackbox optimization software coded in C++. PyNomad is the Python interface to NOMAD.

NOMAD stands for Nonlinear Optimization using the MADS (Mesh Adaptive Direct Search) Algorithm. NOMAD and PyNomad solve constrained optimization problems of blackbox functions in the form:

minimize f(x) 

subject to: nlcon(x) <= 0, lb <= x <= ub and x in R

Once PyNomad is installed, a longer description is obtained by running PyNomad.info(). Examples of PyNomad utilization are provided in https://github.com/bbopt/nomad/tree/master/examples/advanced/library/PyNomad.

**Please cite NOMAD 4 with reference:**

C. Audet, S. Le Digabel, V. Rochon Montplaisir, and C. Tribes.
Algorithm 1027: NOMAD version 4: Nonlinear optimization with the MADS algorithm.
*ACM Transactions on Mathematical Software*.
Volume 48, Issue 3, Article No.: 35, pp 1–22.
https://doi.org/10.1145/3544489
