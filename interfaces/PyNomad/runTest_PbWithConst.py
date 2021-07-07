import PyNomad
import sys

# This example of blackbox function is for a single process
# The blackbox output must be put in the EvalPoint passed as argument
def bb(x):
    try:
        dim = x.size()
        f = x.get_coord(4)
        g1 = sum([(x.get_coord(i)-1)**2 for i in range(dim)])-25
        g2 = 25-sum([(x.get_coord(i)+1)**2 for i in range(dim)])
        rawBBO = str(f) + " " + str(g1) + " " + str(g2)
        x.setBBO(rawBBO.encode("UTF-8"))
    except:
        print("Unexpected eval error", sys.exc_info()[0])
        return 0
    return 1 # 1: success 0: failed evaluation

X0 = [0, 0, 0, 0, 0]
params = ["DIMENSION 5","BB_OUTPUT_TYPE OBJ PB EB", "MAX_BB_EVAL 100", "X0 * 0" , "LOWER_BOUND * -6" , "DISPLAY_DEGREE 2", "DISPLAY_ALL_EVAL false", "DISPLAY_STATS BBE OBJ"]

x_return, f_return, h_return, nb_evals, nb_iters, stopflag = PyNomad.optimize(bb, X0, [] , [], params)
print ("\n NOMAD outputs \n X_sol={} \n F_sol={} \n H_sol={} \n NB_evals={} \n NB_iters={} \n".format(x_return,f_return,h_return,nb_evals,nb_iters))
