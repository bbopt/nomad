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
    double f = 1e+20, g1 = 1e+20, g2 = 1e+20;
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

    g1 = -prod2 + 0.75;
    g2 = sum2 - 7.5 * nb_inputs;

    f = 10 * g1 + 10 * g2;
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

    // for reproducibility
    addNomadValParam(nomad_pb, "NB_THREADS_OPENMP", 1);

    // and the number of blackbox allowed
    addNomadParam(nomad_pb, "MAX_BB_EVAL 1000");

    // run problem
    double x0[10] = {7.0, 7.0, 7.0, 7.0, 7.0,
                     7.0, 7.0, 7.0, 7.0, 7.0}; // starting point

    double x_feas_sol[10] = {0.0, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0, 0.0}; // feasible solution

    double x_inf_sol[10] = {0.0, 0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0, 0.0}; // infeasible solution

    double outputs_feas_sol[10] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                   0.0, 0.0, 0.0, 0.0, 0.0}; // feasible solution outputs

    double outputs_inf_sol[10] = {0.0, 0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 0.0, 0.0}; // infeasible solution outputs

    bool exists_feas, exists_infeas = false; // flag which indicates if the solution exists or not

    solveNomadProblem(nomad_pb, 1, x0,
                      &exists_feas, x_feas_sol, outputs_feas_sol,
                      &exists_infeas, x_inf_sol, outputs_inf_sol,
                      NULL);

    // display found solutions
    if (exists_feas)
    {
        printf("Best feasible solution found: \n");
        printf("x_feas = [ ");
        for (int i = 0; i < nb_inputs; ++i)
        {
            printf("%f ", x_feas_sol[i]);
        }
        printf(" ]\n");
        printf("f_feas = ");
        printf("%f \n", outputs_feas_sol[2]);
    }

    if (exists_infeas) // as a feasible solution has been found, no infeasible solution is given
    {
        printf("Best infeasible solution found (least infeasible with lowest f): \n");
        printf("x_infeas = [ ");
        for (int i = 0; i < nb_inputs; ++i)
        {
            printf("%f ", x_inf_sol[i]);
        }
        printf(" ]\n");
        printf("f_infeas = ");
        printf("%f \n", outputs_inf_sol[2]);
    }

    // NB: relaunch the problem will restart from the beginning
    printf("\n");
    printf("Relaunch second time solving !\n");
    printf("\n");

    addNomadBoolParam(nomad_pb, "DISPLAY_ALL_EVAL", false);

    exists_feas = false;
    exists_infeas = false;
    solveNomadProblem(nomad_pb, 1, x0,
                      &exists_feas, x_feas_sol, outputs_feas_sol,
                      &exists_infeas, x_inf_sol, outputs_inf_sol,
                      NULL);

    // display found solutions
    if (exists_feas)
    {
        printf("Best feasible solution found: \n");
        printf("x_feas = [ ");
        for (int i = 0; i < nb_inputs; ++i)
        {
            printf("%f ", x_feas_sol[i]);
        }
        printf(" ]\n");
        printf("f_feas = ");
        printf("%f \n", outputs_feas_sol[2]);
    }

    if (exists_infeas) // as a feasible solution has been found, no infeasible solution is given
    {
        printf("Best infeasible solution found: \n");
        printf("x_infeas = [ ");
        for (int i = 0; i < nb_inputs; ++i)
        {
            printf("%f ", x_inf_sol[i]);
        }
        printf(" ]\n");
        printf("f_infeas = ");
        printf("%f \n", outputs_inf_sol[2]);
    }

    freeNomadProblem(nomad_pb);

    return EXIT_SUCCESS;
}
