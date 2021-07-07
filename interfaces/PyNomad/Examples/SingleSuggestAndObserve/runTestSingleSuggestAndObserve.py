import PyNomad
import shutil

def bb(x):
    out =[]
    for i in range(len(x)):
        f =[]
        f.append(sum( [ x[i][j] for j in range(len(x[i])) ] ))
        out.append(f)
    return out


# copy cache file
shutil.copyfile("cache0.txt","cache.txt")

params = ["DISPLAY_DEGREE 2", "DIMENSION 3",  "LOWER_BOUND ( -1 -1 -1 )", "UPPER_BOUND * 1", "BB_OUTPUT_TYPE OBJ", "CACHE_FILE cache.txt","MEGA_SEARCH_POLL yes" ]

# Suggest candidates using Mads MegaSearchPoll (a cache file must be provided)
print("Initial parameters:\n",params)
candidates = PyNomad.suggest(params)

# Eval candidates points
evals=bb(candidates)
print("Suggested candidates with evals:\n",candidates,evals)

# Observe new evals -> update parameters and cache
updatedParams = PyNomad.observe(params,candidates,evals,"cache1.txt")

print("Update parameters:\n",updatedParams)
