#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "NomadStdCInterface.h"

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
bool mustache_bb(int nb_inputs, double *x, int nb_outputs, double *bb_outputs, bool *count_eval, NomadUserDataPtr data)
{
    bool eval_ok = true;

    double f = -x[0];

    double g = -(fabs(cos(x[0])) + 0.1) * sin(x[0]) + 2.0;
    double epsilon = 0.05 + 0.05 * (1.0 - 1.0 / (1.0 + fabs(x[0] - 11.0)));

    double c1 = g - epsilon - x[1];
    double c2 = x[1] - g - epsilon;

    bb_outputs[0] = f;
    bb_outputs[1] = c1;
    bb_outputs[2] = c2;

    *count_eval = true;

    return eval_ok;
}

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main(int argc, char **argv)
{
    // indispensable parameters to create the problem
    int nb_inputs = 2;
    int nb_outputs = 3;

    // create Nomad problem
    NomadProblem nomad_pb = createNomadProblem(mustache_bb,
                                               NULL,
                                               nb_inputs,
                                               nb_outputs);

    double lb[] = {0.0, 0.0};
    double ub[] = {20.0, 4.0};
    addNomadArrayOfDoubleParam(nomad_pb, "LOWER_BOUND", lb);
    addNomadArrayOfDoubleParam(nomad_pb, "UPPER_BOUND", ub);

    addNomadParam(nomad_pb,"BB_OUTPUT_TYPE OBJ PB PB");
    addNomadParam(nomad_pb, "DISPLAY_DEGREE 2");
    addNomadParam(nomad_pb, "DISPLAY_STATS EVAL ( SOL ) OBJ CONS_H H_MAX");
    addNomadParam(nomad_pb, "DISPLAY_ALL_EVAL true");
    addNomadParam(nomad_pb, "MAX_BB_EVAL 1000");

    // Essential parameter for this test: stop as soon as a feasible point is found
    addNomadParam(nomad_pb, "STOP_IF_FEASIBLE true");

    double x0[2] = {0.0, 2.0};

    NomadResult nomad_result = createNomadResult();

    int run_flag = solveNomadProblem(nomad_result, nomad_pb, 1, x0, NULL);
    printf("Run status: %d\n", run_flag);

    bool exists_feas = feasibleSolutionsFoundNomadResult(nomad_result) &&
                       nbSolutionsNomadResult(nomad_result) > 0;

    if (!exists_feas)
    {
        printf("Exception: no feasible solution found with STOP_IF_FEASIBLE\n");
        freeNomadProblem(nomad_pb);
        freeNomadResult(nomad_result);
        return EXIT_FAILURE;
    }

    // If a solution was found, load and display it
    double x_sol[2];
    double outputs_sol[2];

    loadInputSolutionsNomadResult(x_sol, 1, nomad_result);
    loadOutputSolutionsNomadResult(outputs_sol, 1, nomad_result);

    printf("Best feasible solution found: \n");
    printf("x = [ %f %f ]\n", x_sol[0], x_sol[1]);
    printf("f_feas = %f \n", outputs_sol[0]);
    printf("Constraints = [ %f %f ] \n", outputs_sol[1], outputs_sol[2]);

    freeNomadProblem(nomad_pb);
    freeNomadResult(nomad_result);
    return EXIT_SUCCESS;
}