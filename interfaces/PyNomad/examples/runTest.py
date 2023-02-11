import PyNomad

# This example of blackbox function is for a single process
# The blackbox output must be put in the EvalPoint passed as argument
def bb(x):
    dim = x.size()
    f = sum([x.get_coord(i)**2 for i in range(dim)])    
    x.setBBO(str(f).encode("UTF-8"))
    return 1 # 1: success 0: failed evaluation

x0 = [0.71, 0.51, 0.51]
lb = [-1, -1, -1]
ub=[]

params = ["BB_OUTPUT_TYPE OBJ", "MAX_BB_EVAL 100", "UPPER_BOUND * 1", "DISPLAY_DEGREE 2", "DISPLAY_ALL_EVAL false", "DISPLAY_STATS BBE OBJ"] 

result = PyNomad.optimize(bb, x0, lb, ub, params)

fmt = ["{} = {}".format(n,v) for (n,v) in result.items()]
output = "\n".join(fmt)
print("\nNOMAD results \n" + output + " \n")
