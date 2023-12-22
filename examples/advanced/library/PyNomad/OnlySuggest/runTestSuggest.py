# import PyNomad
from PyNomad import PyNomadMainStep

#params = ["DISPLAY_DEGREE 2", "DIMENSION 3",  "LOWER_BOUND ( -1 -1 -1 )", "UPPER_BOUND * 1", "BB_OUTPUT_TYPE OBJ", "CACHE_FILE cache.txt", "MEGA_SEARCH_POLL yes" , "SEED 0"]
params = ["DISPLAY_DEGREE 2", "DIMENSION 3",  "LOWER_BOUND ( -1 -1 -1 )", "UPPER_BOUND * 1", "BB_OUTPUT_TYPE OBJ", "CACHE_FILE cache.txt", "LH_EVAL 5" , "SEED 0"]

#candidates = PyNomad.suggest(params)
mainstep = PyNomadMainStep(params)

candidates = mainstep.suggest()

print(candidates)

# To obtain the same candidates the Random Number Generator (static) must be reset
# Otherwise, each suggest call will most likely propose different candidates
# PyNomad.resetRandomNumberGenerator()

#candidates2 = PyNomad.suggest(params)
candidates2 = mainstep.suggest()

print(candidates2)

#print(PyNomad.getRNGState())
