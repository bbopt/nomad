import PyNomad
import shutil

def bb(x):
    out =[]
    for i in range(len(x)):
        f =[]
        f.append(sum( [ x[i][j] for j in range(len(x[i])) ] ))
        # print(f)
        out.append(f)
    return out


# copy cache file
shutil.copyfile("cache0.txt","cache.txt")

params = ["DISPLAY_DEGREE 2", "DIMENSION 3",  "LOWER_BOUND ( -1 -1 -1 )", "UPPER_BOUND * 1", "BB_OUTPUT_TYPE OBJ", "CACHE_FILE cache.txt","INITIAL_FRAME_SIZE ( 1.0 1.0 1.0 )","MEGA_SEARCH_POLL yes" ]

# Stopping criterions on bb eval and on frame size
MAX_BB_EVAL= 100
MIN_FRAME_SIZE= 1E-4
MAX_ITERATIONS= 10

nbEval= 0
frameSize=100000

iteration= 0

for iteration in range(MAX_ITERATIONS):

    # Increment iteration counter
    iteration = iteration + 1

    if nbEval > MAX_BB_EVAL :
        print("Reached max bb eval")
        break

    if frameSize<MIN_FRAME_SIZE:
        print("Reached min frame size")
        break

    print("\n---------------------")
    print("Iteration ",iteration)

    # Suggest candidates using Mads MegaSearchPoll (a cache file must be provided)
    print("Initial parameters:\n",params)
    candidates = PyNomad.suggest(params)

    # Eval candidates points
    evals=bb(candidates)
    print("Suggested candidates with evals:\n",candidates,evals)
    print("---------------------")

    # Account new evaluations
    nbEval = nbEval + len(candidates)


    # Observe new evals -> update parameters and cache
    updatedParams = PyNomad.observe(params,candidates,evals,"cache.txt")

    # Decode bytes into string
    for i in range(len(updatedParams)):
        updatedParams[i] = updatedParams[i].decode('utf-8')
    for i in range(len(params)):
        if type(params[i]) is bytes:
           params[i] = params[i].decode('utf-8')

    print("Updated parameters by observe:\n",updatedParams)
    print("---------------------")

    # Replace updated params in params OR add if not present
    for i in range(len(updatedParams)):
        split1 = updatedParams[i].split()
        found = False
        for j in range(len(params)):
            split2 = params[j].split()
            if ( split2[0].upper() == split1[0].upper() ):
                params[j] = updatedParams[i]
                found = True
                break;
        if not found:
            params.append(updatedParams[i])

    print("Parameters for next iteration:\n",params)
    print("\n")
