import PyNomad

# This example of blackbox function is for a single process
# The blackbox output must be put in the EvalPoint passed as argument
def bb(x):
    dim = x.size()
    f = sum([x.get_coord(i)**2 for i in range(dim)])    
    x.setBBO(str(f).encode("UTF-8"))
    return 1 # 1: success 0: failed evaluation

# Special custom function called at the end of each Mads Mega Iteration.
# Each mega iteration can return several best feasible points.
# Return True if one of the point matches the F target, False otherwise.
def megaIterEndCallback(block):
    try:
        nbPoints = block.size()
        for k in range(nbPoints):
            if block.get_x(k).getF() < 0.00001:
                print('Target F=0.000001 is matched. Let''s stop Nomad')
                print('X=(',block.get_x(k).displayX(),') BBO=(',block.get_x(k).getBBO(),')')
                return True
    except Exception as e:
        print("Unexpected error in megaIterCallback()", str(e))
    return False

x0 = [0.71, 0.51, 0.51]
lb = [-1, -1, -1]
ub=[]

params = ["BB_OUTPUT_TYPE OBJ", "MAX_BB_EVAL 100", "UPPER_BOUND * 1", "DISPLAY_DEGREE 2", "NM_SEARCH false", "DISPLAY_ALL_EVAL false", "DISPLAY_STATS BBE OBJ"]

# The Custom callback function must be passed to Nomad
PyNomad.setCustomMegaIterEndCallback(megaIterEndCallback)

result = PyNomad.optimize(bb, x0, lb, ub, params)

fmt = ["{} = {}".format(n,v) for (n,v) in result.items()]
output = "\n".join(fmt)
print("\nNOMAD results \n" + output + " \n")
