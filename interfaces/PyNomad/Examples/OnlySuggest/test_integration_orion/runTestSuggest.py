import PyNomad


params = ['DISPLAY_DEGREE 2', 'DIMENSION 2', 'BB_INPUT_TYPE ( R R  )', 'BB_OUTPUT_TYPE OBJ ', 'LOWER_BOUND ( 0.0 0.0  )', 'UPPER_BOUND ( 1.0 1.0  )', 'CACHE_FILE cache.txt', 'MEGA_SEARCH_POLL yes', 'INITIAL_FRAME_SIZE ( 0.02 0.02 )', 'H_MAX_0 inf']

PyNomad.setRNGState("1788766504 1007038235 3073713670")

candidates = PyNomad.suggest(params)

print(candidates)


