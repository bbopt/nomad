#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "NomadStdCInterface.h"

bool moustache_bb(int nb_inputs, double *x, int nb_outputs, double *bb_outputs, bool *count_eval, NomadUserDataPtr data)
{
    bool eval_ok = true;

    double f = -x[0];
    double g = -(fabs(cos(x[0])) + 0.1) * sin(x[0]) + 2;
    double epsilon = 0.05 + 0.05 * (1 - 1 / (1 + fabs(x[0] - 11)));

    // fix bb_outputs
    bb_outputs[0] = f;
    bb_outputs[1] = g - epsilon - x[1];
    bb_outputs[2] = x[1] - g - epsilon;

    *count_eval = true;

    return eval_ok;
}

int solve_moustache_pb()
{
    // fix essential parameters of the blackbox
    int nb_inputs = 2;
    int nb_outputs = 3;
    char type_bb_outputs[] = "OBJ PB PB";

    // the problem will terminate after 5000 evaluations
    int max_bb_eval = 5000;

    // fix lower and upper bounds.
    double lb[] = {0.0, 0.0};
    double ub[] = {20.0, 4.0};

    // create Nomad problem
    NomadProblem nomad_pb = createNomadProblem(moustache_bb,
                                               NULL,
                                               nb_inputs,
                                               nb_outputs);

    // fix parameters without the NOMAD terminology

    // main parameters
    addNomadArrayOfDoubleParam(nomad_pb, "LOWER_BOUND", lb);
    addNomadArrayOfDoubleParam(nomad_pb, "UPPER_BOUND", ub);
    addNomadStringParam(nomad_pb, "BB_OUTPUT_TYPE", type_bb_outputs);

    addNomadValParam(nomad_pb, "MAX_BB_EVAL", max_bb_eval);

    // display options
    addNomadValParam(nomad_pb, "DISPLAY_DEGREE", 2);
    addNomadBoolParam(nomad_pb, "DISPLAY_ALL_EVAL", false);
    addNomadBoolParam(nomad_pb, "DISPLAY_UNSUCCESSFUL", false);

    // run problem
    double x0[2] = {0, 2.0}; // starting point

    NomadResult nomad_result = createNomadResult();
    int run_flag = solveNomadProblem(nomad_result, nomad_pb, 1, x0, NULL);
    printf("Run status: %d\n", run_flag);
    const int nb_solutions = nbSolutionsNomadResult(nomad_result);
    const bool exists_feas = feasibleSolutionsFoundNomadResult(nomad_result);
    printf("The solver has found %d solutions ", nb_solutions);
    if (exists_feas)
    {
        printf("and they are feasible\n");
    }
    else
    {
        printf("and they are infeasible\n");
    }

    freeNomadProblem(nomad_pb);
    freeNomadResult(nomad_result);
    return 0;
}

bool speedreducer_bb(int nb_inputs, double *x, int nb_outputs, double *bb_outputs, bool *count_eval, NomadUserDataPtr data)
{
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];
    double x4 = x[3];
    double x5 = x[4];
    double x6 = x[5];
    double x7 = x[6];

    // cost function
    double A = 3.3333 * x3 * x3 + 14.9334 * x3 - 43.0934;
    double B = x6 * x6 + x7 * x7;
    double C = x6 * x6 * x6 + x7 * x7 * x7;
    double D = x4 * x6 * x6 + x5 * x7 * x7;
    double f = 0.7854 * x1 * x2 * x2 * A - 1.508 * x1 * B + 7.477 * C + 0.7854 * D;

    // Constraints
    double A1 = sqrt(745.0 * x4 * 745.0 * x4 + 16900000.0 * x2 * x2 * x3 * x3);
    double A2 = sqrt(745.0 * x5 * 745.0 * x5 + 157500000.0 * x2 * x2 * x3 * x3);
    double g1 = 27.0 - x1 * x2 * x2 * x3;                          // <= 0
    double g2 = 397.5 - x1 * x2 * x2 * x3 * x3;                    // <= 0
    double g3 = 1.93 * x4 * x4 * x4 - x2 * x6 * x6 * x6 * x6 * x3; // <= 0
    double g4 = 1.93 * x5 * x5 * x5 - x2 * x7 * x7 * x7 * x7 * x3; // <= 0
    double g5 = A1 - 110.0 * x2 * x3 * x6 * x6 * x6;               // <= 0
    double g6 = A2 - 85.0 * x2 * x3 * x7 * x7 * x7;                // <= 0
    double g7 = x2 * x3 - 40.0;                                    // <= 0
    double g8 = 5.0 * x2 - x1;                                     // <= 0
    double g9 = x1 - 12 * x2;                                      // <= 0
    double g10 = 1.9 + 1.5 * x6 - x4;                              // <= 0
    double g11 = 1.9 + 1.1 * x7 - x5;                              // <= 0

    bb_outputs[0] = f;
    bb_outputs[1] = g1;
    bb_outputs[2] = g2;
    bb_outputs[3] = g3;
    bb_outputs[4] = g4;
    bb_outputs[5] = g5;
    bb_outputs[6] = g6;
    bb_outputs[7] = g7;
    bb_outputs[8] = g8;
    bb_outputs[9] = g9;
    bb_outputs[10] = g10;
    bb_outputs[11] = g11;

    *count_eval = true;

    // the evaluation has succeeded
    return true;
}

void speedreducer_bb_block(int block_size, int nb_inputs, double *x,
                            int nb_outputs, double *bb_outputs,
                            bool *count_eval, bool *eval_ok, NomadUserDataPtr data)
{
    double *inputs = malloc(nb_inputs * sizeof(double));
    double *outputs = malloc(nb_outputs * sizeof(double));
    for (int index = 0; index < block_size; ++index)
    {
        for (int i = 0; i < nb_inputs; ++i)
        {
            inputs[i] = x[index * nb_inputs + i];
        }
        // Call the blackbox on each element of the block
        // There could be some applications where it is faster
        // to parallelize the blocks
        eval_ok[index] = speedreducer_bb(nb_inputs, inputs,
                                         nb_outputs, outputs,
                                         &count_eval[index], data);
        for (int i = 0; i < nb_outputs; ++i)
        {
            bb_outputs[index * nb_outputs + i] = outputs[i];
        }
    }
    free(inputs);
    free(outputs);
}

int solve_speedreducer_pb()
{
    // fix essential parameters of the blackbox
    int dim = 7;
    int nb_outputs = 12;
    char type_bb_outputs[] = "OBJ PB PB PB PB PB PB PB PB PB PB PB";

    // the problem will terminate after 5000 evaluations
    int max_bb_eval = 5000;

    // fix lower and upper bounds.
    double lb[7] = {2.6, 0.7, 17, 7.3, 7.3, 2.9, 5.0};
    double ub[7] = {3.6, 0.8, 28, 8.3, 8.3, 3.9, 5.5};

    // create Nomad problem
    NomadProblem nomad_pb = createNomadProblem(speedreducer_bb,
                                               speedreducer_bb_block,
                                               dim,
                                               nb_outputs);

    // fix the parameters without the NOMAD terminology

    // main parameters
    addNomadArrayOfDoubleParam(nomad_pb, "LOWER_BOUND", lb);
    addNomadArrayOfDoubleParam(nomad_pb, "UPPER_BOUND", ub);
    addNomadStringParam(nomad_pb, "BB_OUTPUT_TYPE", type_bb_outputs);

    addNomadValParam(nomad_pb, "MAX_BB_EVAL", max_bb_eval);

    // display options
    addNomadValParam(nomad_pb, "DISPLAY_DEGREE", 2);
    addNomadBoolParam(nomad_pb, "DISPLAY_ALL_EVAL", false);
    addNomadBoolParam(nomad_pb, "DISPLAY_UNSUCCESSFUL", false);

    // set non opportunistic eval
    addNomadBoolParam(nomad_pb, "EVAL_OPPORTUNISTIC", false);

    // Activate evaluation per block
    addNomadValParam(nomad_pb, "BB_MAX_BLOCK_SIZE", 4);

    // run problem
    double x0[7] = {3.000000000000000e+00,
                    7.500000000000000e-01,
                    22,
                    8.000000000000000e+00,
                    8.000000000000000e+00,
                    3.400000000000000e+00,
                    5.250000000000000e+00};

    NomadResult nomad_result = createNomadResult();
    int run_flag = solveNomadProblem(nomad_result,
                                     nomad_pb, 1, x0,
                                     NULL);
    printf("Run status: %d\n", run_flag);
    const int nb_solutions = nbSolutionsNomadResult(nomad_result);
    const bool exists_feas = feasibleSolutionsFoundNomadResult(nomad_result);
    printf("The solver has found %d solutions ", nb_solutions);
    if (exists_feas)
    {
        printf("and they are feasible\n\n");
    }
    else
    {
        printf("and they are infeasible\n\n");
    }

    freeNomadProblem(nomad_pb);
    freeNomadResult(nomad_result);

    return 0;
}

int main()
{
    solve_moustache_pb();
    solve_speedreducer_pb();
    return EXIT_SUCCESS;
}
