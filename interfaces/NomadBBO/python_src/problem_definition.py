"""
Problem definition for CatMads optimization.

This module contains the main class for defining optimization problems
with categorical and continuous variables.
"""

import math
from numbers import Real as NumberReal
from typing import Dict, Any, List, Optional
from .variables import BaseVariable


class ProblemDefinition:
    """
    Class for defining a optimization problem,
    which can contain categorical, integer, binary and continuous variables.
    
    This class encapsulates the problem variables, blackbox output types,
    initial points,and provides the interface for problem evaluation.
    Initial points can be provided in the problem definition 
    or through the optimizer's DOE or both.
    """
    
    def __init__(
        self,
        vars: Dict[str, BaseVariable],
        bbot: List[str],
        X0s: Optional[List[Any]] = None,
        **kwargs,
    ):
        """
        Initialize a CatMads optimization problem.

        Parameters:
        - vars: Dict[str, BaseVariable]
            A dictionary defining the problem's variables. Each key is the name of the variable,
            and the value must be an instance derived from BaseVariable (e.g., Binary, Integer, Real, Choice).
            Example:
            vars = {
                "b": Binary(),
                "x": Choice(options=["nothing", "multiply"]),
                "y": Integer(bounds=(0, 2)),
                "z": Real(bounds=(0, 5)),
            }

        - bbot: List[str]
            A list of strings defining the black box output types. Only "OBJ" (objective) and 
            "PB" (constraint/problem) are allowed.
            Example: ["OBJ", "PB", "PB"]

        - runParams: Dict[str, Any], optional
            A dictionary of Nomad parameters to customize the optimization.
            If not provided, will check kwargs for 'runParams'. 
            If neither is found, uses default values.
            Example: {"DISPLAY_DEGREE": 2, "MAX_BB_EVAL": 100}

        - X0s: optional
            Optional initial point or list of initial points in raw NOMAD format.
            Each point must match the problem dimension.

        - kwargs: dict
            Any additional parameters for customization.
        """
        # Validate that all variables are instances of BaseVariable
        for var_name, var_instance in vars.items():
            if not isinstance(var_instance, BaseVariable):
                raise TypeError(f"Variable '{var_name}' must be an instance of BaseVariable or its subclasses. "
                               f"Got {type(var_instance).__name__} instead.")
        
        # Validate that bbot is a list of valid strings
        if not isinstance(bbot, list):
            raise TypeError(f"Parameter 'bbot' must be a list of strings, got {type(bbot).__name__} instead.")
        
        allowed_bbot_values = {"OBJ", "PB"}
        for i, value in enumerate(bbot):
            if not isinstance(value, str):
                raise TypeError(f"All elements in 'bbot' must be strings. Element at index {i} is {type(value).__name__}.")
            if value not in allowed_bbot_values:
                raise ValueError(f"Invalid value '{value}' in 'bbot' at index {i}. Only 'OBJ' and 'PB' are allowed.")
        
        self.vars = vars
        self.bbot = bbot
    
        if 'initialDOE' in kwargs:
            raise TypeError("initialDOE now belongs to the Optimizer/CatMadsOptimizer, not CatMadsProblemDefinition.")

        self.X0s: Optional[List[List[float]]] = None
        if X0s is not None:
            self.setX0(X0s)


    def _normalize_x0_point(self, point: List[Any]) -> List[float]:
        """
        Validate one X0 point in raw NOMAD format.
        """
        if not isinstance(point, (list, tuple)):
            raise TypeError(f"Each X0 point must be a list or tuple, got {type(point).__name__}.")

        if len(point) != self.getDimension():
            raise ValueError(f"Each X0 point must have dimension {self.getDimension()}, got {len(point)}.")

        normalized_point = []
        for index, value in enumerate(point):
            if isinstance(value, bool) or not isinstance(value, NumberReal):
                raise TypeError(f"X0[{index}] must be numeric, got {type(value).__name__}.")

            num_value = float(value)
            if not math.isfinite(num_value):
                raise ValueError(f"X0[{index}] must be finite, got {value}.")
            normalized_point.append(num_value)

        self.convertPointToMixedVariableInput(normalized_point)
        return normalized_point

    def setX0(self, X0s: Optional[List[Any]]) -> None:
        """
        Set optional initial points X0.

        Accepted formats:
        - one point: [x1, x2, ..., xn]
        - several points: [[...], [...], ...]
        """
        if X0s is None:
            self.X0s = None
            return

        if not isinstance(X0s, (list, tuple)):
            raise TypeError(f"X0s must be a point or a list of points, got {type(X0s).__name__}.")

        if len(X0s) == 0:
            self.X0s = []
            return

        first_item = X0s[0]
        if isinstance(first_item, (list, tuple)):
            normalized_points = [self._normalize_x0_point(point) for point in X0s]
        else:
            normalized_points = [self._normalize_x0_point(X0s)]

        self.X0s = normalized_points

    def getX0s(self) -> Optional[List[List[float]]]:
        """
        Return a copy of the optional initial points X0s.
        """
        if self.X0s is None:
            return None
        return [point.copy() for point in self.X0s]

    def evaluate(self, X: Dict[str, Any]) -> Dict[str, Any]:
        """
        Evaluate the objective function for a given solution.

        Parameters:
        - X: dict
            A dictionary where keys correspond to variable names and values are their current state.

        Returns:
        - A dictionary with the computed objectives and any other relevant data.
        """
        raise NotImplementedError("Subclass must implement the evaluate method.")
    
    def getVariables(self) -> List[BaseVariable]:
        """
        Get a list of all variables in the problem definition.
        
        Returns:
        - List[BaseVariable]: List of all variable instances.
        """
        return list(self.vars.values())
    
    def getVariableNames(self) -> List[str]:
        """
        Get a list of all variable names in the problem definition.
        
        Returns:
        - List[str]: List of all variable names.
        """
        return list(self.vars.keys())
    
    def getBBOutputType(self) -> List[str]:
        """
        Get the black box output types.
        
        Returns:
        - List[str]: List of black box output types (e.g., ["OBJ", "PB", "PB"]).
        """
        return self.bbot.copy()  # Return a copy to prevent external modification
    
    def getDimension(self) -> int:
        """
        Get the total dimension of the problem (number of variables).
        
        Returns:
        - int: Total number of variables in the problem.
        """
        return len(self.vars)
    
    def getVariableByName(self, name: str) -> BaseVariable:
        """
        Get a specific variable by name.
        
        Parameters:
        - name: str
            The name of the variable to retrieve.
            
        Returns:
        - BaseVariable: The variable instance.
        
        Raises:
        - KeyError: If the variable name is not found.
        """
        if name not in self.vars:
            raise KeyError(f"Variable '{name}' not found in problem definition.")
        return self.vars[name]
    
    def getVariableTypes(self) -> Dict[str, str]:
        """
        Get the types of all variables.
        
        Returns:
        - Dict[str, str]: Dictionary mapping variable names to their types.
        """
        return {name: var.type for name, var in self.vars.items()}
    
    def getVariablesByType(self, var_type: str = None) -> Dict[str, List[tuple]]:
        """
        Get variables grouped by their type, or filter by a specific type.
        
        Parameters:
        - var_type: str, optional
            If specified, return only variables of this type.
            If None, return all variables grouped by type.
        
        Returns:
        - Dict[str, List[tuple]]: Dictionary where keys are variable types 
          and values are lists of (name, variable_instance) tuples.
          If var_type is specified, returns dict with only that type.
        
        Example:
        >>> problem.getVariablesByType()
        {
            'binary': [('x1', Binary()), ('x2', Binary())],
            'choice': [('y1', Choice(['A', 'B'])), ('y2', Choice(['X', 'Y', 'Z']))],
            'integer': [('z1', Integer((0, 10)))],
            'real': [('w1', Real((0.0, 5.0)))]
        }
        
        >>> problem.getVariablesByType('binary')
        {'binary': [('x1', Binary()), ('x2', Binary())]}
        """
        grouped_vars = {}
        
        for name, var in self.vars.items():
            if var.type not in grouped_vars:
                grouped_vars[var.type] = []
            grouped_vars[var.type].append((name, var))
        
        if var_type is not None:
            if var_type not in grouped_vars:
                return {var_type: []}
            return {var_type: grouped_vars[var_type]}
        
        return grouped_vars
    
    def getVariableIndexByName(self, name: str) -> int:
        """
        Get the index (position) of a variable by its name.
        
        Parameters:
        - name: str
            The name of the variable.
            
        Returns:
        - int: The index (0-based) of the variable in the order they were defined.
        
        Raises:
        - KeyError: If the variable name is not found.
        
        Example:
        >>> problem = CatMadsProblemDefinition(vars={"x": Binary(), "y": Choice([...]), "z": Real(...)})
        >>> problem.getVariableIndexByName("y")  # Returns 1
        """
        variable_names = list(self.vars.keys())
        if name not in variable_names:
            raise KeyError(f"Variable '{name}' not found in problem definition.")
        return variable_names.index(name)
    
    def getVariableIndicesByType(self, var_type: str = None) -> Dict[str, List[int]]:
        """
        Get the indices of variables grouped by their type, or filter by a specific type.
        
        Parameters:
        - var_type: str, optional
            If specified, return only indices of variables of this type.
            If None, return indices of all variables grouped by type.
        
        Returns:
        - Dict[str, List[int]]: Dictionary where keys are variable types 
          and values are lists of indices (0-based positions).
          If var_type is specified, returns dict with only that type.
        
        Example:
        >>> problem = CatMadsProblemDefinition(vars={
        ...     "x1": Binary(), "y1": Choice([...]), "x2": Binary(), "z1": Real(...)
        ... })
        >>> problem.getVariableIndicesByType()
        {'binary': [0, 2], 'choice': [1], 'real': [3]}
        
        >>> problem.getVariableIndicesByType('binary')
        {'binary': [0, 2]}
        """
        variable_names = list(self.vars.keys())
        grouped_indices = {}
        
        for index, (name, var) in enumerate(self.vars.items()):
            if var.type not in grouped_indices:
                grouped_indices[var.type] = []
            grouped_indices[var.type].append(index)
        
        if var_type is not None:
            if var_type not in grouped_indices:
                return {var_type: []}
            return {var_type: grouped_indices[var_type]}
        
        return grouped_indices
    
    def getVariableBounds(self, exclude_binary: bool = False) -> Dict[str, tuple]:
        """
        Get the bounds for all variables (excluding binary if specified).
        
        Parameters:
        - exclude_binary: bool
            If True, exclude binary variables from the result.
            
        Returns:
        - Dict[str, tuple]: Dictionary mapping variable names to their bounds.
                           Binary variables return (0, 1).
                           Choice variables return (0, len(options)-1).
        """
        bounds = []
        for name, var in self.vars.items():
            if exclude_binary and var.type == "binary":
                continue
                
            if var.type == "binary":
                bounds.append((0, 1))
            elif var.type == "choice":
                # Choice variables: bounds based on number of options
                bounds.append((var.options[0], var.options[-1]))
            elif hasattr(var, 'bounds'):
                # Integer and Real variables have bounds attribute
                bounds.append(var.bounds)
            else:
                bounds.append((None, None))  # No bounds defined
            
        return bounds


    def getVariableTypeSummary(self) -> Dict[str, Any]:
        idx = self.getVariableIndicesByType()
        cat_idx = idx.get("choice", [])          # categorical = choice only
        bin_idx = idx.get("binary", [])
        quant_idx = idx.get("integer", []) + idx.get("real", [])
        return {
            "nbCat": len(cat_idx),
            "nbBin": len(bin_idx),
            "nbQuant": len(quant_idx),
            "cat_idx": cat_idx,
            "bin_idx": bin_idx,
            "quant_idx": quant_idx,
            "dimension": self.getDimension(),
            "variable_names": self.getVariableNames(),
        }

    def getNbCategoriesPerCategoricalVariable(self) -> List[int]:
        return [len(var.options) for var in self.vars.values() if var.type == "choice"]

    def getQuantitativeBounds(self) -> List[tuple]:
        bounds = []
        for var in self.vars.values():
            if var.type in ("integer", "real"):
                bounds.append(var.bounds)
        return bounds

    def convertPointToMixedVariableInput(self, point, variable_names: List[str] = None) -> Dict[str, Any]:
        """
        Convert a point (list of numeric values) to a dictionary compatible 
        with the evaluate function.
        
        Parameters:
        - point: List-like object with numeric values from NOMAD
        - variable_names: List[str], optional
            Order of variable names. If None, uses the order from self.vars.keys()
        
        Returns:
        - Dict[str, Any]: Dictionary with variable names as keys and properly 
                         converted values for the evaluate function.
        """
        if variable_names is None:
            variable_names = list(self.vars.keys())
        
        if len(point) != len(variable_names):
            raise ValueError(f"Point has {len(point)} values but {len(variable_names)} variables expected.")
        
        result = {}
        
        for i, var_name in enumerate(variable_names):
            if var_name not in self.vars:
                raise ValueError(f"Variable '{var_name}' not found in problem definition.")
            
            var = self.vars[var_name]
            raw_value = point[i]
            
            # Convert based on variable type
            if var.type == "binary":
                # Binary: convert to int (0 or 1)
                # If raw_value is not exactly 0 or 1, return an error
                if raw_value not in (0, 1):
                    raise ValueError(f"Binary variable '{var_name}' has invalid value {raw_value}. Expected 0 or 1.")
                result[var_name] = raw_value
            elif var.type == "choice":
                # Choice: convert index to actual choice value
                index = int(round(raw_value))
                result[var_name] = var.indexToValue(index)
            elif var.type == "integer":
                # Integer: round to nearest integer
                result[var_name] = int(round(raw_value))
            elif var.type == "real":
                # Real: keep as float
                result[var_name] = float(raw_value)
            else:
                # Unknown type: keep as is
                result[var_name] = raw_value
        
        return result

