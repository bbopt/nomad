import PyNomad

def objective (point):
    return point.get_coord(4)

def square (point, shift):
    return sum((point.get_coord(index) + shift) ** 2 
        for index in range(point.size()))

def constraint_one (point):
    return square(point, - 1) - 25

def constraint_two (point):
    return 25 - square(point, + 1)

def blackbox (eval_point):
    eval_f0 = objective(eval_point)
    eval_g0 = constraint_one(eval_point)
    eval_g1 = constraint_one(eval_point)

    eval_point.setBBO(f'{eval_f0} {eval_g0} {eval_g1}'.encode('utf-8'))
    return True

def test_single_constraint ():
    x0 = [0, 0, 0, 0, 0]
    params = [
        'DIMENSION 5',
        'BB_OUTPUT_TYPE OBJ EB EB',
        'MAX_BB_EVAL 100',
        'X0 * 0' ,
        'LOWER_BOUND * -6' ,
        'DISPLAY_DEGREE 2',
        'DISPLAY_ALL_EVAL false',
        'DISPLAY_STATS BBE OBJ'
    ]

    result = PyNomad.optimize(blackbox, x0, [], [], params)

    x_result = result['x_best']
    x_expect = [ -0.22, 1.53, 2.66, 1.52, -3.41 ]

    assert result['exit_status'] == 1
    for result_value, expect_value in zip(x_result, x_expect):
        assert abs(result_value - expect_value) < 1e-6

