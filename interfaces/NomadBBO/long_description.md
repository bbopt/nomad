NOMAD is a blackbox optimization software coded in C++. NomadBBO is the Python interface to NOMAD for solving problems.

NOMAD stands for Nonlinear Optimization using the MADS (Mesh Adaptive Direct Search) Algorithm. NomadBBO can solve mixed variables constrained optimization problems of blackbox functions in the form:

minimize f(x) 

subject to: nlcon(x) <= 0, lb <= x <= ub

where: x is a mixed-variable point x = (x^cat, x^int, x^cont)

Once NomadBBO is installed, a longer description is obtained by running NomadBBO.info(). 

**Please cite NOMAD 4 with reference:**

C. Audet, S. Le Digabel, V. Rochon Montplaisir, and C. Tribes.
Algorithm 1027: NOMAD version 4: Nonlinear optimization with the MADS algorithm.
*ACM Transactions on Mathematical Software*.
Volume 48, Issue 3, Article No.: 35, pp 1–22.
https://doi.org/10.1145/3544489
