import PyNomad

def objective (point):
    return sum(point.get_coord(index) ** 2 for index in range(point.size()))

def blackbox (eval_point):
    eval_value = objective(eval_point)
    eval_point.setBBO(str(eval_value).encode('utf-8'))
    return True

def test_single ():
    x0 = [ 0.71, 0.51, 0.51 ]
    lb = [ -1, -1, -1 ]
    ub = []

    params = [
        'BB_OUTPUT_TYPE OBJ',
        'UPPER_BOUND * 1',
        'DISPLAY_DEGREE 1',
        'DISPLAY_ALL_EVAL false',
        'DISPLAY_STATS BBE OBJ'
    ] 

    result = PyNomad.optimize(blackbox, x0, lb, ub, params)

    assert result['run_flag'] == 1
    for value in result['x_best']:
        assert abs(value) < 1e-6

