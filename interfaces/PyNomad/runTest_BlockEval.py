import PyNomad
import sys

# This example is for a block, i.e., a vector of EvalPoints.
# Evaluations may be conducted in parallel. The parallelism must be managed by the
# user.
# The blackbox output must be put in each EvalPoint passed as argument.
def bb_block(block):
    nbPoints = block.size()
    evalOk = [False for i in range(nbPoints)]
    try:
        for k in range(nbPoints):
            # eval each point
            x = block.get_x(k)
            dim = x.size()
            f = sum([x.get_coord(i)**2 for i in range(dim)])
            x.setBBO(str(f).encode("UTF-8"), B"OBJ")
            evalOk[k] = True
    except Exception as e:
        print("Unexpected error in bb_block()", str(e))
        return 0
    return evalOk # list where 1 is success, 0 is a failed evaluation

x0 = [0.71, 0.51, 0.51]
lb = [-1, -1, -1]
ub=[]

params =  ["BB_OUTPUT_TYPE OBJ", "MAX_BB_EVAL 100", "UPPER_BOUND * 1"]
params += ["DISPLAY_DEGREE 2", "DISPLAY_STATS BBE BLK_SIZE OBJ", "DISPLAY_ALL_EVAL false"]
params += ["NB_THREADS_OPENMP 1", "BB_MAX_BLOCK_SIZE 4"]

x_return, f_return, h_return, nb_evals, nb_iters, stopflag = PyNomad.optimize(bb_block, x0, lb, ub, params)
print ("\n NOMAD outputs \n X_sol={} \n F_sol={} \n H_sol={} \n NB_evals={} \n NB_iters={} \n".format(x_return,f_return,h_return,nb_evals,nb_iters))
