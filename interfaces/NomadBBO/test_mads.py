import sys
import os

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import
from NomadBBO import (
    ProblemDefinition,
    QuantitativeOptimizer,
    Choice,
    Integer,
    Real,
    optimize,
)

class MyMadsProblem(ProblemDefinition):

    def __init__(self, **kwargs):
        vars = {
            # Order of variables for now: categorical, integer and continuous variables
            "x": Integer(bounds=(0, 4)),  # Integer variable
            "i": Integer(bounds=(0, 3)),  # Integer variable
            "y": Integer(bounds=(0, 2)),  # Integer variable
            "z": Real(bounds=(0, 5)),  # Continuous variable
        }
        X0s = [ 0, 0, 0, 0.0]  # Initial point for the optimization (raw format)
        super().__init__(vars=vars,
                        bbot=["OBJ", "PB", "PB"],
                        X0s=X0s,
                        **kwargs)

    def evaluate(self, X):
        i, x, y, z = X["i"], X["x"], X["y"], X["z"]

        # -------------------------
        # Simple base objective
        # -------------------------
        f = 100.0 * i * (y + z)

        #
        # Categories are coded by integers
        #
        
        # -------------------------
        # Category effects (all different,
        # but some are more similar than others)
        # A and B: very similar
        # C: intermediate
        # D and E: more distinct
        # -------------------------
        if x == 0:
            f = 1.00 * f + 10.0
        elif x == 1:
            f = 1.02 * f + 12.0   # very close to A (correlated category)
        elif x == 2:
            f = 0.95 * f + 25.0   # moderately different
        elif x == 3:
            f = 1.15 * f + 5.0    # clearly different behaviour
        elif x == 4:
            f = 0.85 * f + 40.0   # most different

        g1_base = (-z + y) + 0.35
        g2_base = (z / 2.0 - y) + 0.25

        # Category-dependent shifts (A~B similar, others different)
        if x == 0:
            s1, s2 = 0.00, 0.00
        elif x == 1:
            s1, s2 = 0.03, 0.02   # very close to A
        elif x == 2:
            s1, s2 = 0.10, 0.08
        elif x == 3:
            s1, s2 = -0.12, 0.18
        elif x == 4:
            s1, s2 = 0.18, -0.10

        g1 = g1_base + s1
        g2 = g2_base + s2

        return {"OBJ": f, "PB": [g1, g2]}


# Test the imports and problem definition creation
if __name__ == "__main__":
    print("Testing Mads imports and problem definition...")
    
    # Create an instance of the problem
    try:
        problem = MyMadsProblem()
        print("✓ Problem instance created successfully")
        
        # Test the variable types
        print(f"✓ Variables defined: {list(problem.vars.keys())}")
        print(f"✓ Variable types: {problem.getVariableTypes()}")
        print(f"✓ Variable bounds: {problem.getVariableBounds()}")
        
        
        
        # Test evaluation
        nomad_point = [1, 1, 2, 3.5]  # Example NOMAD point
        eval_input = problem.convertPointToMixedVariableInput(nomad_point) # Conversion needed to match the variable order in the problem definition
        print(f"✓ NOMAD point conversion: {nomad_point} → {eval_input}")

        result = problem.evaluate(eval_input)
        print(f"✓ Evaluation result: {result}")
        
        print("\nAll tests passed! Imports are working correctly.")
        
    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()

    print("\n Let's run the optimization")
    # Create an instance of the problem for optimization
    problem_instance = MyMadsProblem()
    
    # By default
    optimizer = QuantitativeOptimizer(isBOSearchUsed=True)
    result = optimize(problem_instance, optimizer, budget=100)


    fmt = ["{} = {}".format(n,v) for (n,v) in result.items()]
    output = "\n".join(fmt)
    print("\nNOMAD results \n" + output + " \n")
