import PyNomad


params = ["DISPLAY_DEGREE 2", "DIMENSION 3", "X0 (0.71 0.51 0.51)", "LOWER_BOUND ( -1 -1 -1 )", "UPPER_BOUND * 1", "MAX_BB_EVAL 100", "BB_EXE \"$python bb.py\"" , "BB_OUTPUT_TYPE OBJ"]
# params = ["DISPLAY_DEGREE 2"]


PyNomad.optimizeWithMainStep(params)
