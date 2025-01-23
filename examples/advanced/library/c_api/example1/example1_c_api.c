#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "NomadStdCInterface.h"

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
bool example1_bb(int nb_inputs, double *x, int nb_outputs, double *bb_outputs, bool *count_eval, NomadUserDataPtr data)
{
    bool eval_ok = true;

    // based on G2
    double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, prod1 = 1.0, prod2 = 1.0;

    for (int i = 0; i < nb_inputs; ++i)
    {
        sum1 += pow(cos(x[i]), 4);
        sum2 += x[i];
        sum3 += (i + 1) * x[i] * x[i];
        prod1 *= pow(cos(x[i]), 2);
        if (prod2 != 0.0)
        {
            if (x[i] == 0.0)
            {
                prod2 = 0.0;
            }
            else
            {
                prod2 *= x[i];
            }
        }
    }

    const double g1 = -prod2 + 0.75;
    const double g2 = sum2 - 7.5 * nb_inputs;

    double f = 10 * g1 + 10 * g2;
    if (0.0 != sum3)
    {
        f -= fabs(((sum1 - 2 * prod1) / sqrt(sum3)));
    }

    // Scale
    f *= 1e-5;

    double c2000 = -f - 2000;

    // Fix outputs
    bb_outputs[0] = g1;
    bb_outputs[1] = g2;
    bb_outputs[2] = f;
    bb_outputs[3] = c2000;

    *count_eval = true;

    return eval_ok;
}

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main(int argc, char **argv)
{
    // indispensable parameters to create the problem
    int nb_inputs = 10;
    int nb_outputs = 4;

    // create Nomad problem
    NomadProblem nomad_pb = createNomadProblem(example1_bb,
                                               NULL,
                                               nb_inputs,
                                               nb_outputs);

    double granularity_values[10] = {0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001,
                                     0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001};

    // Fix parameters using NOMAD convention

    // fix important parameters
    addNomadParam(nomad_pb,"BB_OUTPUT_TYPE PB PB OBJ EB");

    // fix some external parameters
    addNomadArrayOfDoubleParam(nomad_pb, "GRANULARITY",  granularity_values);

    addNomadParam(nomad_pb, "DISPLAY_DEGREE 2");
    addNomadParam(nomad_pb, "DISPLAY_STATS EVAL ( SOL ) OBJ CONS_H H_MAX");
    addNomadParam(nomad_pb, "DISPLAY_ALL_EVAL true");
    addNomadParam(nomad_pb, "DISPLAY_UNSUCCESSFUL false");

    // and the number of blackbox allowed
    addNomadParam(nomad_pb, "MAX_BB_EVAL 1000");

    // run problem
    double x0[10] = {7.0, 7.0, 7.0, 7.0, 7.0,
                     7.0, 7.0, 7.0, 7.0, 7.0}; // starting point

    NomadResult nomad_result = createNomadResult();

    int run_flag = solveNomadProblem(nomad_result, nomad_pb, 1, x0, NULL);
    printf("Run status: %d\n", run_flag);

    int nb_solutions = nbSolutionsNomadResult(nomad_result);
    printf("The algorithm has found %d solutions\n", nb_solutions);

    // Get a solution
    double x_sol[10];
    double outputs_sol[4];

    bool exists_feas = feasibleSolutionsFoundNomadResult(nomad_result) &&
                       nb_solutions > 0;
    loadInputSolutionsNomadResult(x_sol, 1, nomad_result);
    loadOutputSolutionsNomadResult(outputs_sol, 1, nomad_result);

    // display found solutions
    if (exists_feas)
    {
        printf("Best feasible solution found: \n");
        printf("x_feas = [ ");
        for (int i = 0; i < nb_inputs; ++i)
        {
            printf("%f ", x_sol[i]);
        }
        printf(" ]\n");
        printf("f_feas = ");
        printf("%f \n", outputs_sol[2]);
        printf("Constraints = [ %f %f %f ]\n",
               outputs_sol[0], outputs_sol[1], outputs_sol[3]);
    }
    else
    {
        printf("Best infeasible solution found (least infeasible with lowest f): \n");
        printf("x_infeas = [ ");
        for (int i = 0; i < nb_inputs; ++i)
        {
            printf("%f ", x_sol[i]);
        }
        printf(" ]\n");
        printf("f_infeas = ");
        printf("%f \n", outputs_sol[2]);
        printf("Constraints = [ %f %f %f ]\n",
               outputs_sol[0], outputs_sol[1], outputs_sol[3]);
    }

    // NB: relaunch the problem will restart from the beginning
    printf("\n");
    printf("Relaunch second time solving !\n");
    printf("\n");

    addNomadBoolParam(nomad_pb, "DISPLAY_ALL_EVAL", false);
    // Set the dependant parameters to their default value
    // Needed because they have been set in the previous run and are not reset by default
    addNomadBoolParam(nomad_pb, "DISPLAY_INFEASIBLE", true);
    addNomadBoolParam(nomad_pb, "DISPLAY_UNSUCCESSFUL", false);

    run_flag = solveNomadProblem(nomad_result, nomad_pb, 1, x0, NULL);
    printf("Run status: %d\n", run_flag);

    nb_solutions = nbSolutionsNomadResult(nomad_result);
    exists_feas = feasibleSolutionsFoundNomadResult(nomad_result) &&
                  nb_solutions > 1;
    loadInputSolutionsNomadResult(x_sol, 1, nomad_result);
    loadOutputSolutionsNomadResult(outputs_sol, 1, nomad_result);

    // display found solutions
    if (exists_feas)
    {
        printf("Best feasible solution found: \n");
        printf("x_feas = [ ");
        for (int i = 0; i < nb_inputs; ++i)
        {
            printf("%f ", x_sol[i]);
        }
        printf(" ]\n");
        printf("f_feas = ");
        printf("%f \n", outputs_sol[2]);
        printf("c(x) = [ %f %f %f ]\n",
               outputs_sol[0], outputs_sol[1], outputs_sol[3]);
    }
    else
    {
        printf("Best infeasible solution found: \n");
        printf("x_infeas = [ ");
        for (int i = 0; i < nb_inputs; ++i)
        {
            printf("%f ", x_sol[i]);
        }
        printf(" ]\n");
        printf("f_infeas = ");
        printf("%f \n", outputs_sol[2]);
        printf("c(x) = [ %f %f %f ]\n",
               outputs_sol[0], outputs_sol[1], outputs_sol[3]);
    }

    freeNomadProblem(nomad_pb);
    freeNomadResult(nomad_result);

    return EXIT_SUCCESS;
}
