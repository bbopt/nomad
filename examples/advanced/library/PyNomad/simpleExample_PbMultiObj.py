import PyNomad
import sys

# This example of blackbox function is for a single process
# The blackbox output must be put in the EvalPoint passed as argument
def bb(x):
    try:
        x0 = x.get_coord(0)
        x1 = x.get_coord(1)
        f1 = (x0-1.0)*(x0-1.0)+(x0-x1)*(x0-x1)
        f2 = (x0-x1) * (x0-x1) + (x1-3) * (x1-3)
        rawBBO = str(f1) + " " + str(f2)
        x.setBBO(rawBBO.encode("UTF-8"))
    except:
        print("Unexpected eval error", sys.exc_info()[0])
        return 0
    return 1 # 1: success 0: failed evaluation

X0 = [0, 0]
params = ["DIMENSION 2","BB_OUTPUT_TYPE OBJ OBJ", "DMULTIMADS_OPTIMIZATION yes", "DIRECTION_TYPE ORTHO 2N", "MAX_BB_EVAL 400", "X0 * 2.0" , "LOWER_BOUND * -1" ,"UPPER_BOUND * 5" , "DISPLAY_DEGREE 2", "SOLUTION_FILE SOL.txt"]

# NOTES:
# Result returned by PyNomad.optimize() return a single point
# Solution file contains pareto points

result = PyNomad.optimize(bb, X0, [] , [], params)

fmt = ["{} = {}".format(n,v) for (n,v) in result.items()]
output = "\n".join(fmt)
print("\nNOMAD results \n" + output + " \n")
