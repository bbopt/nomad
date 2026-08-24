"""
DOE (Design of Experiment) utilities for CatMads optimization.

This module provides the DOEMixin class that adds initial DOE management
capabilities to CatMadsProblemDefinition.
"""

import math
import random
from abc import ABC, abstractmethod
from numbers import Real as NumberReal
from typing import Dict, Any, List, Optional


def _lhs_integer_column(n_samples: int, lb: int, ub: int, rng: random.Random) -> List[int]:
    """
    Generate a permuted-stratified (LHS-style) column for an integer variable in [lb, ub].

    When n_samples <= number of options, each stratum maps to a distinct integer.
    When n_samples > number of options, full coverage is repeated and the remainder
    is filled with a random subset, then the whole column is shuffled.
    """
    n_options = ub - lb + 1
    if n_options == 1:
        return [lb] * n_samples

    if n_samples <= n_options:
        strata_width = n_options / n_samples
        permutation = list(range(n_samples))
        rng.shuffle(permutation)
        values: List[int] = []
        for k in permutation:
            val = lb + int(k * strata_width + rng.random() * strata_width)
            values.append(min(val, ub))
    else:
        repeats = n_samples // n_options
        remainder = n_samples % n_options
        pool = list(range(lb, ub + 1)) * repeats
        pool += rng.sample(range(lb, ub + 1), remainder)
        rng.shuffle(pool)
        values = pool

    return values


class DOEgenerator(ABC):
    """
    Generic DOE specification.

    Subclasses must implement `generate_raw_points` and return a list of points
    in raw NOMAD numeric format.
    """

    @abstractmethod
    def generate_raw_points(self, problemDefinition, n_samples: int) -> List[List[float]]:
        """Generate DOE points in raw NOMAD format."""

    def apply(self, owner, problemDefinition,  n_samples: int) -> List[Dict[str, Any]]:
        """
        Build validated DOE records using owner validation helpers and store them in owner.initial_doe.
        """
        raw_points = self.generate_raw_points(problemDefinition, n_samples)
        owner.initial_doe = [owner.validateDOEPoint(problemDefinition, point) for point in raw_points]
        return owner.getInitialDOE()


class LatinHypercubeDOE(DOEgenerator):
    """Latin Hypercube DOE."""

    def __init__(self, seed: Optional[int] = 0):
        self.seed = seed

    def generate_raw_points(self, problemDefinition, n_samples) -> List[List[float]]:
        rng = random.Random(self.seed)
        if not isinstance(n_samples, int) or n_samples < 1:
            raise ValueError(f"n_samples must be a positive integer, got {n_samples!r}.")
        self.n_samples = n_samples
        variables = list(problemDefinition.vars.items())
        n_vars = len(variables)

        columns: List[List[float]] = []
        for var_name, var in variables:
            if var.type == "binary":
                n_ones = (self.n_samples + 1) // 2
                n_zeros = self.n_samples - n_ones
                pool = [0.0] * n_zeros + [1.0] * n_ones
                rng.shuffle(pool)
                columns.append(pool)

            elif var.type == "choice":
                n_options = len(var.options)
                col = _lhs_integer_column(self.n_samples, 0, n_options - 1, rng)
                columns.append([float(v) for v in col])

            elif var.type == "integer":
                lb, ub = var.bounds
                col = _lhs_integer_column(self.n_samples, int(lb), int(ub), rng)
                columns.append([float(v) for v in col])

            elif var.type == "real":
                lb, ub = var.bounds
                if lb == ub:
                    columns.append([float(lb)] * self.n_samples)
                else:
                    strata_width = (ub - lb) / self.n_samples
                    permutation = list(range(self.n_samples))
                    rng.shuffle(permutation)
                    col = []
                    for k in permutation:
                        val = lb + (k + rng.random()) * strata_width
                        col.append(max(lb, min(ub, val)))
                    columns.append(col)

            else:
                raise ValueError(f"Unsupported variable type '{var.type}' for variable '{var_name}'.")

        return [[columns[j][i] for j in range(n_vars)] for i in range(self.n_samples)]



class DOEMixin:
    """
    Mixin providing initial-point management utilities.

    The initial design can combine:
    - user-provided initial points, stored in the problem definition as X0s;
    - automatically generated DOE points, typically used to complete the initial design.

    In the optimization interface, the target DOE size is computed from the problem
    and budget. If the user-provided X0s are fewer than this target, the remaining
    points are generated with the selected DOE generator, e.g., LatinHypercubeDOE.

    This mixin provides validation and formatting helpers for manually provided DOE
    points.
    """

    def _normalize_numeric_value(self, value: Any, context: str) -> float:
        """
        Normalize a numeric value to float and reject NaN/Inf.
        """
        if isinstance(value, bool) or not isinstance(value, NumberReal):
            raise TypeError(f"{context} must be numeric, got {type(value).__name__}.")
        num_value = float(value)
        if not math.isfinite(num_value):
            raise ValueError(f"{context} must be finite, got {value}.")
        return num_value

    def _validate_raw_point(self, problemDefinition, x_raw: List[Any]) -> List[float]:
        """
        Validate a DOE point against variable definitions and return normalized numeric values.
        """
        if not isinstance(x_raw, (list, tuple)):
            raise TypeError(f"DOE point must be a list or tuple, got {type(x_raw).__name__}.")

        dimension = problemDefinition.getDimension()
        if len(x_raw) != dimension:
            raise ValueError(f"DOE point has dimension {len(x_raw)}, expected {dimension}.")

        normalized = []
        for index, (var_name, var) in enumerate(problemDefinition.vars.items()):
            raw_value = self._normalize_numeric_value(x_raw[index], f"DOE value for variable '{var_name}'")

            if var.type == "binary":
                if raw_value not in (0.0, 1.0):
                    raise ValueError(f"Binary variable '{var_name}' must be 0 or 1, got {x_raw[index]}.")
                normalized_value = int(raw_value)
            elif var.type == "choice":
                rounded = int(round(raw_value))
                if not math.isclose(raw_value, rounded, abs_tol=1e-12):
                    raise ValueError(
                        f"Choice variable '{var_name}' must use an integer index in [0, {len(var.options) - 1}], got {x_raw[index]}."
                    )
                if rounded < 0 or rounded >= len(var.options):
                    raise ValueError(
                        f"Choice variable '{var_name}' index {rounded} out of bounds [0, {len(var.options) - 1}]."
                    )
                normalized_value = rounded
            elif var.type == "integer":
                rounded = int(round(raw_value))
                if not math.isclose(raw_value, rounded, abs_tol=1e-12):
                    raise ValueError(
                        f"Integer variable '{var_name}' must be an integer value in bounds {var.bounds}, got {x_raw[index]}."
                    )
                lb, ub = var.bounds
                if rounded < lb or rounded > ub:
                    raise ValueError(
                        f"Integer variable '{var_name}' value {rounded} out of bounds [{lb}, {ub}]."
                    )
                normalized_value = rounded
            elif var.type == "real":
                lb, ub = var.bounds
                if raw_value < lb or raw_value > ub:
                    raise ValueError(
                        f"Real variable '{var_name}' value {raw_value} out of bounds [{lb}, {ub}]."
                    )
                normalized_value = raw_value
            else:
                raise ValueError(f"Unsupported variable type '{var.type}' for variable '{var_name}'.")

            normalized.append(float(normalized_value))

        return normalized

    def _validate_bbot_values(self, problemDefinition, y: Optional[List[Any]]) -> Optional[List[float]]:
        """
        Validate optional BBOT evaluations.
        """
        if y is None:
            return None

        if not isinstance(y, (list, tuple)):
            raise TypeError(f"DOE evaluations must be a list or tuple, got {type(y).__name__}.")

        expected = len(problemDefinition.bbot)
        if len(y) != expected:
            raise ValueError(f"DOE evaluations length is {len(y)}, expected {expected} (len(bbot)).")

        normalized_y = []
        for index, value in enumerate(y):
            bbot_name = problemDefinition.bbot[index]
            normalized_y.append(self._normalize_numeric_value(value, f"DOE evaluation y[{index}] for BBOT '{bbot_name}'"))

        return normalized_y

    def validateDOEPoint(self, problemDefinition, x_raw: List[Any], y: Optional[List[Any]] = None) -> Dict[str, Any]:
        """
        Validate and normalize a DOE point.

        Returns a normalized record:
        {
            "x_raw": List[float],
            "x_mixed": Dict[str, Any],
            "y": Optional[List[float]],
        }
        """
        normalized_x = self._validate_raw_point(problemDefinition, x_raw)
        mixed_point = problemDefinition.convertPointToMixedVariableInput(normalized_x)
        normalized_y = self._validate_bbot_values(problemDefinition, y)

        return {
            "x_raw": normalized_x,
            "x_mixed": mixed_point,
            "y": normalized_y,
        }

    def setInitialDOE(self, problemDefinition, doe: List[Dict[str, Any]]) -> None:
        """
        Set the initial DOE after full validation.

        Accepted input formats per entry:
        - {"x": [...], "y": [...]} or {"x_raw": [...], "y": [...]} or {"point": [...], "eval": [...]}.
        - y/eval is optional.
        """
        if not isinstance(doe, list):
            raise TypeError(f"initialDOE must be a list of points, got {type(doe).__name__}.")

        normalized_doe = []
        for index, entry in enumerate(doe):
            if not isinstance(entry, dict):
                raise TypeError(f"initialDOE[{index}] must be a dict, got {type(entry).__name__}.")

            x_raw = entry.get("x")
            if x_raw is None:
                x_raw = entry.get("x_raw")
            if x_raw is None:
                x_raw = entry.get("point")
            if x_raw is None:
                raise ValueError(f"initialDOE[{index}] must define one of keys 'x', 'x_raw' or 'point'.")

            y = entry.get("y")
            if y is None and "eval" in entry:
                y = entry.get("eval")

            normalized_doe.append(self.validateDOEPoint(problemDefinition, x_raw=x_raw, y=y))

        self.initial_doe = normalized_doe

    def addDOEPoint(self, problemDefinition, x_raw: List[Any], y: Optional[List[Any]] = None) -> Dict[str, Any]:
        """
        Validate and append one DOE point.
        """
        record = self.validateDOEPoint(problemDefinition, x_raw=x_raw, y=y)
        self.initial_doe.append(record)
        return record

    def getInitialDOE(self) -> List[Dict[str, Any]]:
        """
        Return a shallow copy of the normalized DOE list.
        """
        return [
            {
                "x_raw": point["x_raw"].copy(),
                "x_mixed": point["x_mixed"].copy(),
                "y": None if point["y"] is None else point["y"].copy(),
            }
            for point in self.initial_doe
        ]

    def createLatinHypercubeDOE(self, problemDefinition, n_samples: int, seed: Optional[int] = None) -> List[Dict[str, Any]]:
        """
        Generate an initial DOE using Latin Hypercube Sampling and set it as self.initial_doe.

        Each variable type is handled appropriately:
        - real:    continuous stratified sampling within [lb, ub].
        - integer: stratified sampling rounded to valid integers within [lb, ub].
        - choice:  stratified sampling over the integer index space [0, n_options - 1].
        - binary:  balanced sampling over {0, 1}.

        Parameters
        ----------
        n_samples : int
            Number of sample points to generate (must be >= 1).
        seed : int, optional
            Random seed for reproducibility.

        Returns
        -------
        List[Dict[str, Any]]
            The generated DOE (same list stored in self.initial_doe).
        """
        generator = LatinHypercubeDOE(seed=seed)
        return generator.apply(self, problemDefinition, n_samples)
