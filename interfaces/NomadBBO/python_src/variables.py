"""
Variable types for CatMads optimization problems.

This module defines the base variable class and all specific variable types
that can be used in CatMads optimization problems.
"""

from typing import List, Union


class BaseVariable:
    """
    Base class for all variable types.
    This class should be inherited by specific variable classes.
    """
    
    def __init__(self, var_type: str):
        """
        Initialize a base variable.

        Parameters:
        - var_type: str
            A string representing the type of the variable (e.g., 'binary', 'choice').
        """
        self.type = var_type

    def validate(self, value):
        """
        Abstract method for validation. Should be implemented in subclasses to validate values.

        Parameters:
        - value: Any
            Value to be validated.

        Raises:
        - NotImplementedError
            If the method is not implemented by a subclass.
        """
        raise NotImplementedError("Validate method must be implemented in subclass.")
    
    def getBounds(self) -> tuple:
        """
        Get the bounds for this variable. Should be implemented in subclasses.
        
        Returns:
        - tuple: (lower_bound, upper_bound)
        """
        raise NotImplementedError("getBounds method must be implemented in subclass.")
    
    def getInitialValue(self):
        """
        Get a default initial value for this variable.
        
        Returns:
        - Default value appropriate for this variable type.
        """
        raise NotImplementedError("getInitialValue method must be implemented in subclass.")


class Binary(BaseVariable):
    """Binary variable that can take values 0 or 1."""
    
    def __init__(self):
        super().__init__("binary")

    def validate(self, value):
        if value not in [0, 1]:
            raise ValueError(f"Invalid value '{value}' for Binary variable. It must be 0 or 1.")
        return value
    
    def getBounds(self) -> tuple:
        return (0, 1)
    
    # For now, this should not be used in CatMads optimization. 
    # A DOE will provide initial values.
    def getInitialValue(self):
        return 0


class Choice(BaseVariable):
    """Categorical variable that can take values from a predefined list of options."""
    
    def __init__(self, options: List[str]):
        if not isinstance(options, list) or len(options) == 0:
            raise ValueError("Choice options must be a non-empty list.")
        super().__init__("choice")
        self.options = options

    def validate(self, value):
        if value not in self.options:
            raise ValueError(f"Invalid value '{value}' for Choice variable. Allowed options: {self.options}")
        return value
    
    def getBounds(self) -> tuple:
        return (0, len(self.options) - 1)
    
    def indexToValue(self, index: int) -> str:
        """
        Convert an integer index to the corresponding choice value.
        
        Parameters:
        - index: int
            The index (0-based) of the option to select.
            
        Returns:
        - str: The corresponding option value.
        
        Raises:
        - ValueError: If index is out of bounds.
        """
        if not isinstance(index, int) or index < 0 or index >= len(self.options):
            raise ValueError(f"Invalid index {index} for Choice variable. "
                           f"Must be integer between 0 and {len(self.options)-1}.")
        return self.options[index]
    
    def valueToIndex(self, value: str) -> int:
        """
        Convert a choice value to its corresponding integer index.
        
        Parameters:
        - value: str
            The choice value to convert.
            
        Returns:
        - int: The corresponding index (0-based).
        
        Raises:
        - ValueError: If value is not in options.
        """
        if value not in self.options:
            raise ValueError(f"Value '{value}' not found in options: {self.options}")
        return self.options.index(value)
    
    # For now, this should not be used in CatMads optimization. 
    # A DOE will provide initial values.
    def getInitialValue(self):
        return 0  # Index of first option


class Integer(BaseVariable):
    """Integer variable with specified bounds."""
    
    def __init__(self, bounds: Union[tuple, list]):
        if not isinstance(bounds, (tuple, list)) or len(bounds) != 2:
            raise ValueError("Integer bounds must be a tuple or list with two elements.")
        super().__init__("integer")
        self.bounds = bounds

    def validate(self, value):
        if not isinstance(value, int) or not (self.bounds[0] <= value <= self.bounds[1]):
            raise ValueError(f"Invalid value '{value}' for Integer variable. "
                           f"It must be between {self.bounds[0]} and {self.bounds[1]}.")
        return value
    
    def getBounds(self) -> tuple:
        return self.bounds
    
    # For now, this should not be used in CatMads optimization. 
    # A DOE will provide initial values.
    def getInitialValue(self):
        return self.bounds[0]  # Lower bound as initial value


class Real(BaseVariable):
    """Real (continuous) variable with specified bounds."""
    
    def __init__(self, bounds: Union[tuple, list]):
        if not isinstance(bounds, (tuple, list)) or len(bounds) != 2:
            raise ValueError("Real bounds must be a tuple or list with two elements.")
        super().__init__("real")
        self.bounds = bounds

    def validate(self, value):
        if not isinstance(value, (int, float)) or not (self.bounds[0] <= value <= self.bounds[1]):
            raise ValueError(f"Invalid value '{value}' for Real variable. "
                           f"It must be between {self.bounds[0]} and {self.bounds[1]}.")
        return value
    
    def getBounds(self) -> tuple:
        return self.bounds
    
    # For now, this should not be used in CatMads optimization. 
    # A DOE will provide initial values.
    def getInitialValue(self):
        return self.bounds[0]  # Lower bound as initial value