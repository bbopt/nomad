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
 \file   example2_c_api.c
 \brief  c api example 2 for Nomad (test different blackboxes)
 \author Ludovic Salomon
 \date   2020
 */
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "Interfaces/NomadStdCInterface.h"

bool moustache_bb(int nb_inputs, double *x, int nb_outputs, double *bb_outputs, bool *count_eval, NomadUserDataPtr data)
{
    bool eval_ok = false;

    double f = -x[0];
    double g = -(fabs(cos(x[0])) + 0.1) * sin(x[0]) + 2;
    double epsilon = 0.05 + 0.05 * (1 - 1 / (1 + fabs(x[0] - 11)));

    // fix bb_outputs
    bb_outputs[0] = f;
    bb_outputs[1] = g - epsilon - x[1];
    bb_outputs[2] = x[1] - g - epsilon;

    // count evaluation
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
                                               nb_inputs,
                                               nb_outputs,
                                               lb, ub,
                                               NULL, // no precision of type variable
                                               type_bb_outputs,
                                               max_bb_eval);

    setNomadDisplayDegree(nomad_pb, 2);
    setNomadDisplayAllEval(nomad_pb, false);
    setNomadDisplayUnsuccessful(nomad_pb, false);

    // run problem
    double x0[2] = {0, 2.0}; // starting point

    double x_feas_sol[2] = {0.0, 0.0}; // feasible solution

    double x_inf_sol[2] = {0.0, 0.0}; // infeasible solution

    double outputs_feas_sol[3] = {0.0, 0.0, 0.0}; // feasible solution outputs

    double outputs_inf_sol[3] = {0.0, 0.0, 0.0}; // infeasible solution outputs

    bool exists_feas, exists_infeas = false; // flag which indicates if the solution exists or not

    solveNomadProblem(nomad_pb, x0,
                      &exists_feas, x_feas_sol, outputs_feas_sol,
                      &exists_infeas, x_inf_sol, outputs_inf_sol,
                      NULL);
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

    // count eval
    *count_eval = true;

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

    // the evaluation has succeeded
    return true;
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
                                               dim,
                                               nb_outputs,
                                               lb, ub,
                                               NULL, // no precision of type variable
                                               type_bb_outputs,
                                               max_bb_eval);

    setNomadDisplayDegree(nomad_pb, 2);
    setNomadDisplayAllEval(nomad_pb, false);
    setNomadDisplayUnsuccessful(nomad_pb, false);

    // set non opportunistic eval
    setNomadOpportunisticEval(nomad_pb, false);

    // run problem
    double x0[7] = {3.000000000000000e+00,
                    7.500000000000000e-01,
                    22,
                    8.000000000000000e+00,
                    8.000000000000000e+00,
                    3.400000000000000e+00,
                    5.250000000000000e+00};

    double x_feas_sol[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // feasible solution

    double x_inf_sol[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // infeasible solution

    double outputs_feas_sol[12] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // feasible solution outputs

    double outputs_inf_sol[12] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // infeasible solution outputs

    bool exists_feas, exists_infeas = false; // flag which indicates if the solution exists or not

    solveNomadProblem(nomad_pb, x0,
                      &exists_feas, x_feas_sol, outputs_feas_sol,
                      &exists_infeas, x_inf_sol, outputs_inf_sol,
                      NULL);

    return 0;
}

int main(int argc, char **argv)
{
    solve_moustache_pb();
    solve_speedreducer_pb();
    return EXIT_SUCCESS;
}