import PyNomad

import shutil

# Copy cache file
shutil.copyfile("cache0.txt","cache.txt")

params = ["DISPLAY_DEGREE 2", "DIMENSION 3",  "LOWER_BOUND ( -1 -1 -1 )", "UPPER_BOUND * 1", "MAX_BB_EVAL 100", "BB_OUTPUT_TYPE OBJ", "CACHE_FILE cache.txt","MEGA_SEARCH_POLL yes" ]


points=[[0.5,0,0.5]]
evals=[[-1]]

print("Initial parameters:\n",params)

# Observe evaluated points and update cache and parameters
updatedParams = PyNomad.observe(params,points,evals,"cache1.txt")

# Decode bytes into string
for i in range(len(updatedParams)):
    updatedParams[i] = updatedParams[i].decode('utf-8')

print("Modified parameters:\n",updatedParams)
