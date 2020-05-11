/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural            */
/*  Sciences and Engineering Research Council of Canada), InnovÉÉ (Innovation      */
/*  en Énergie Électrique) and IVADO (The Institute for Data Valorization)         */
/*                                                                                 */
/*  NOMAD v3 was created and developed by Charles Audet, Sebastien Le Digabel,     */
/*  Christophe Tribes and Viviane Rochon Montplaisir and was funded by AFOSR       */
/*  and Exxon Mobil.                                                               */
/*                                                                                 */
/*  NOMAD v1 and v2 were created and developed by Mark Abramson, Charles Audet,    */
/*  Gilles Couture, and John E. Dennis Jr., and were funded by AFOSR and           */
/*  Exxon Mobil.                                                                   */
/*                                                                                 */
/*  Contact information:                                                           */
/*    Polytechnique Montreal - GERAD                                               */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada              */
/*    e-mail: nomad@gerad.ca                                                       */
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
/*                                                                                 */
/*  This program is free software: you can redistribute it and/or modify it        */
/*  under the terms of the GNU Lesser General Public License as published by       */
/*  the Free Software Foundation, either version 3 of the License, or (at your     */
/*  option) any later version.                                                     */
/*                                                                                 */
/*  This program is distributed in the hope that it will be useful, but WITHOUT    */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or          */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License    */
/*  for more details.                                                              */
/*                                                                                 */
/*  You should have received a copy of the GNU Lesser General Public License       */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.           */
/*                                                                                 */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad           */
/*---------------------------------------------------------------------------------*/
/**
 \file   example1_c_api.c
 \brief  c api example for Nomad
 \author Ludovic Salomon
 \date   2020
 */
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "Interfaces/NomadStdCInterface.h"

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
bool example1_bb(int nb_inputs, double *x, int nb_outputs, double *bb_outputs, bool *count_eval, NomadUserDataPtr data)
{
    bool eval_ok = false;

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

    // count evaluation
    *count_eval = true;

    return eval_ok;
}

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main(int argc, char **argv)
{
    // fix essential parameters of the blackbox
    int nb_inputs = 10;
    int nb_outputs = 4;
    char type_bb_outputs[] = "PB PB OBJ EB";

    // the problem will terminate after 1000 evaluations
    int max_bb_eval = 1000;

    // create Nomad problem
    NomadProblem nomad_pb = createNomadProblem(example1_bb,
                                               nb_inputs,
                                               nb_outputs,
                                               NULL, NULL, // no lower neither upper bound
                                               NULL,       // no precision of type variable
                                               type_bb_outputs,
                                               max_bb_eval);

    // fix some external parameters
    double granularity_params[10] = {0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001,
                                     0.0000001, 0.0000001, 0.0000001, 0.0000001, 0.0000001};
    setNomadGranularityBBInputs(nomad_pb, granularity_params);

    setNomadDisplayDegree(nomad_pb, 2);
    setNomadDisplayAllEval(nomad_pb, false);
    setNomadDisplayUnsuccessful(nomad_pb, false);

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

    solveNomadProblem(nomad_pb, x0,
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
        printf("Best feasible solution found: \n");
        printf("x_infeas = [ ");
        for (int i = 0; i < nb_inputs; ++i)
        {
            printf("%f ", x_inf_sol[i]);
        }
        printf(" ]\n");
        printf("f_feas = ");
        printf("%f \n", outputs_inf_sol[2]);
    }

    // NB: relaunch the problem will restart from the beginning
    printf("\n");
    printf("Relaunch second time solving !\n");
    printf("\n");
    solveNomadProblem(nomad_pb, x0,
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

    return EXIT_SUCCESS;
}