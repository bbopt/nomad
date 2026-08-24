"""
Public Python API for NomadBBO.
"""

from .NomadBBO_src import optimize

from .python_src import (
    BaseVariable,
    Binary,
    Choice,
    Integer,
    Real,
    ProblemDefinition,
    DataManager,
    MixedDataManager,
    QuantitativeDataManager,
    BaseOptimizer,
    MixedOptimizer,
    QuantitativeOptimizer,
    DOEgenerator,
    LatinHypercubeDOE,
)

__all__ = [
    "optimize",
    "BaseVariable",
    "Binary",
    "Choice",
    "Integer",
    "Real",
    "ProblemDefinition",
    "DataManager",
    "MixedDataManager",
    "QuantitativeDataManager",
    "BaseOptimizer",
    "MixedOptimizer",
    "QuantitativeOptimizer",
    "DOEgenerator",
    "LatinHypercubeDOE",
]

__version__ = "1.0b4"