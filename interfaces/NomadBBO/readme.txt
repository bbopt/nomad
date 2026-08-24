**********************
**********************
NOMADBBO
**********************
**********************

NomadBBO is a Python package for mixed-variable blackbox optimization
built on top of NOMAD.

It supports:

- Continuous variables
- Integer variables
- Categorical variables
- Blackbox objective functions
- Blackbox constraints

NomadBBO provides a Pythonic interface based on:

- ProblemDefinition
- CatMadsOptimizer
- optimize()

instead of requiring users to manually construct NOMAD parameter files.


**********************
INSTALLATION
**********************

pip install NomadBBO


**********************
QUICK START
**********************

from CatMads import *
import NomadBBO

problem = ProblemDefinition(...)
optimizer = CatMadsOptimizer()

result = NomadBBO.optimize(
    problemDefinition=problem,
    optimizer=optimizer,
    budget=1000
)