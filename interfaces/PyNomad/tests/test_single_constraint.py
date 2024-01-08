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
    eval_g1 = constraint_two(eval_point)

    eval_point.setBBO(f'{eval_f0} {eval_g0} {eval_g1}'.encode('utf-8'))
    return True

def test_single_constraint ():
    x0 = [0, 0, 0, 0, 0]
    params = [
        'DIMENSION 5',
        'MAX_BB_EVAL 400',
        'BB_OUTPUT_TYPE OBJ EB EB',
        'X0 * 0' ,
        'LOWER_BOUND * -6' ,
        'DISPLAY_DEGREE 2',
        'DISPLAY_ALL_EVAL false',
        'DISPLAY_STATS BBE OBJ'
    ]

    result = PyNomad.optimize(blackbox, x0, [], [], params)

    f_result = result['f_best']
    f_expect = -4.0

    assert result['run_flag'] == 0
    assert abs(f_expect - f_result) < 1e-1

