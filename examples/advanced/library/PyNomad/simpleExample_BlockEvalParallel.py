import PyNomad
import sys

from asyncio import gather, run


# This example is for a block, i.e., a vector of EvalPoints, evaluated asynchronously.

# Evaluation of a single point.
async def bb_single(x):
    try:
        dim = x.size()
        f = sum([x.get_coord(i)**2 for i in range(dim)])
    except Exception as e:
        print("Unexpected error in bb_single()", str(e))
        return "Inf"
    return str(f)

# The blackbox output must be put in each EvalPoint passed as argument (x.setBBO()).
async def async_bb_block(block):
    nbPoints = block.size()
    evals = []
    evalOk = []
    for i in range(nbPoints):
        evalOk.append(False)
        xs = block.get_x(i)
        evals.append(bb_single(xs))
    
    try:
        results = await gather(*evals)
        for k in range(nbPoints):
            x = block.get_x(k)
            x.setBBO(results[k].encode("UTF-8"))
            evalOk[k] = True
    except Exception as e:
        print("Unexpected error in async_bb_block()", str(e))
        return 0
    return evalOk # list where 1 is success, 0 is a failed evaluation

def sync_bb_block(block):
    result = run(async_bb_block(block))
    return result


x0 = [0.71, 0.51, 0.51]
lb = [-1, -1, -1]
ub=[]

params =  ["BB_OUTPUT_TYPE OBJ", "MAX_BB_EVAL 100", "UPPER_BOUND * 1", "DIRECTION_TYPE ORTHO N+1 NEG"]
params += ["DISPLAY_DEGREE 2", "DISPLAY_STATS BBE BLK_SIZE OBJ", "DISPLAY_ALL_EVAL true"]
params += ["MEGA_SEARCH_POLL yes", "BB_MAX_BLOCK_SIZE 4"]

result = PyNomad.optimize(sync_bb_block, x0, lb, ub, params)

fmt = ["{} = {}".format(n,v) for (n,v) in result.items()]
output = "\n".join(fmt)
print("\nNOMAD results \n" + output + " \n")
