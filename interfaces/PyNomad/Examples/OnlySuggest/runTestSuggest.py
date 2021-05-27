import PyNomad

params = ["DISPLAY_DEGREE 2", "DIMENSION 3",  "LOWER_BOUND ( -1 -1 -1 )", "UPPER_BOUND * 1", "BB_OUTPUT_TYPE OBJ", "CACHE_FILE cache.txt", "MEGA_SEARCH_POLL yes"]

candidates = PyNomad.suggest(params)

print(candidates)
