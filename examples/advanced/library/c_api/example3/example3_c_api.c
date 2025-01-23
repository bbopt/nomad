#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "NomadStdCInterface.h"

// A conceptual marine design problem
// See "An Easy-To-use Real-world Multi-objective Optimization Suite problem"
// by R. Tanabe and H. Ishibuchi
bool marine_design_bb(int nb_inputs, double *x, int nb_outputs, double *bb_outputs, bool *count_eval, NomadUserDataPtr data)
{
    bool eval_ok = true;

    // Variables
    const double L = x[0];
    const double B = x[1];
    const double D = x[2];
    const double T = x[3];
    const double Vk = x[4];
    const double Cb = x[5];

    // Parameters
    const double handlingRate = 8000;
    const double roundTripMiles = 5000;
    const double fuelPrice = 100;
    const double g = 9.8065;

    // Equations
    const double verticalCenterOfBuoyancy = 0.53 * T;
    const double metacentricRadius = (0.085 * Cb - 0.002) * B * B / (T * Cb);
    const double verticalCenterOfGravity = 1 + 0.52 * D;

    const double displacement = 1.025 * L * B * T * Cb;
    const double a = 4977.06 * Cb * Cb - 8105.61 * Cb + 4456.51;
    const double b = -10847.2 * Cb * Cb + 12817 * Cb - 6960.32;
    const double V = 0.5144 * Vk;
    const double froudenumber = V / pow(g * L, 0.5);
    const double P = pow(displacement, 2.0/3) * pow(Vk, 3) / (a + b * froudenumber);

    const double steelWeight = 0.034 * pow(L, 1.7) * pow(B, 0.7) *
        pow(D, 0.4) * pow(Cb, 0.5);
    const double outfitWeight = 1.0 * pow(L, 0.8) * pow(B, 0.6) *
        pow(D, 0.3) * pow(Cb, 0.1);
    const double machineryWeight = 0.17 * pow(P, 0.9);
    const double lightShipWeight = steelWeight + outfitWeight + machineryWeight;

    const double deadweight = displacement - lightShipWeight;

    const double dailyConsumption = (0.19 * P * 24) / 1000 + 0.2;
    const double seaDays = (roundTripMiles * Vk) / 24;
    const double fuelCarried = dailyConsumption * (seaDays + 5);
    const double portCost = 6.3 * pow(deadweight, 0.8);

    const double cargoDeadweight = deadweight - fuelCarried - 2 * pow(deadweight, 0.5);
    const double portDays = 2 * (cargoDeadweight / handlingRate + 0.5);
    const double roundTripsPerYear = 350.0 / (seaDays + portDays);

    const double fuelCost = 1.05 * dailyConsumption * seaDays * fuelPrice;
    const double runningCosts = 40000 * pow(deadweight, 0.8);
    const double voyageCosts = (fuelCost + portCost) * roundTripsPerYear;
    const double shipCosts = 1.3 * (2000 * pow(steelWeight, 0.85) +
            3500 * outfitWeight + 2400 * pow(P, 0.8));
    const double capitalCosts = 0.2 * shipCosts;
    const double annualCosts = capitalCosts + runningCosts + voyageCosts;

    const double annualCargo = - cargoDeadweight * roundTripsPerYear;

    // Objectives
    const double f1 = annualCosts / annualCargo;
    const double f2 = lightShipWeight;
    const double f3 = annualCargo;

    // Constraints
    const double g1 = L / B - 6.0; // >= 0
    const double g2 = 15.0 - L / D; // >= 0
    const double g3 = 19.0 - L / T; // >= 0
    const double g4 = 0.45 * pow(deadweight, 0.31) - T; // >= 0
    const double g5 = 0.7 * D + 0.7 - T; // >= 0
    const double g6 = deadweight - 3000; // >= 0
    const double g7 = 500000 - deadweight; // >= 0
    const double g8 = 0.32 - froudenumber; // >= 0
    const double g9 = verticalCenterOfBuoyancy + metacentricRadius -
        verticalCenterOfGravity - 0.07 * B; // >= 0

    // fix bb_outputs
    bb_outputs[0] = f1;
    bb_outputs[1] = f2;
    bb_outputs[2] = f3;
    bb_outputs[3] = -g1;
    bb_outputs[4] = -g2;
    bb_outputs[5] = -g3;
    bb_outputs[6] = -g4;
    bb_outputs[7] = -g5;
    bb_outputs[8] = -g6;
    bb_outputs[9] = -g7;
    bb_outputs[10] = -g8;
    bb_outputs[11] = -g9;

    *count_eval = true;

    return eval_ok;
}

int solve_marine_design_pb()
{
    // fix essential parameters of the blackbox
    int nb_inputs = 6;
    int nb_outputs = 12;
    char type_bb_outputs[] = "OBJ OBJ OBJ PB PB PB PB PB PB PB PB PB";

    // the problem will terminate after 1500 evaluations
    int max_bb_eval = 1500;

    // fix lower and upper bounds.
    double lb[] = {150, 20, 13, 10, 14, 0.63};
    double ub[] = {274.32, 32.31, 25, 11.71, 18, 0.75};

    // create Nomad problem
    NomadProblem nomad_pb = createNomadProblem(marine_design_bb,
                                               NULL,
                                               nb_inputs,
                                               nb_outputs);

    // Fix parameters without the NOMAD terminology

    // Main parameters
    addNomadArrayOfDoubleParam(nomad_pb, "LOWER_BOUND", lb);
    addNomadArrayOfDoubleParam(nomad_pb, "UPPER_BOUND", ub);
    addNomadStringParam(nomad_pb, "BB_OUTPUT_TYPE", type_bb_outputs);

    addNomadValParam(nomad_pb, "MAX_BB_EVAL", max_bb_eval);

    // Display options
    addNomadValParam(nomad_pb, "DISPLAY_DEGREE", 2);
    addNomadBoolParam(nomad_pb, "DISPLAY_ALL_EVAL", false);
    addNomadBoolParam(nomad_pb, "DISPLAY_UNSUCCESSFUL", false);

    // As we are in a multiobjective context, we need to explicitly
    // set the choice of the algorithm used
    addNomadBoolParam(nomad_pb, "DMULTIMADS_OPTIMIZATION", true);

    // For the multiobjective case, we cannot use ORTHO N+1 QUAD
    addNomadParam(nomad_pb, "DIRECTION_TYPE ORTHO N+1 NEG");

    // Here, deactivate QUAD_MODEL_SEARCH
    addNomadBoolParam(nomad_pb, "QUAD_MODEL_SEARCH", false);

    // Change options for NM strategy
    addNomadStringParam(nomad_pb, "DMULTIMADS_NM_STRATEGY", "MULTI");

    // A line initialization is practically more efficient than giving a single point for
    // multiobjective optimization.
    // The interesting reader can report to the following reference for more information.
    // Direct Multisearch for multiobjective optimization
    // by A.L. Custodio, J.F.A. Madeira, A.I.F. Vaz and L.N. Vicente, 2011.
    double x0s[6 * 6]; // starting points
    for (size_t j = 0; j < nb_inputs; ++j)
    {
        for (size_t i = 0; i < nb_inputs; ++i)
        {
            x0s[j * nb_inputs + i] = lb[i] + (double) j * (ub[i] - lb[i]) / (nb_inputs - 1);
        }
    }

    // Run problem
    NomadResult nomad_result = createNomadResult();
    int run_flag = solveNomadProblem(nomad_result, nomad_pb,
                                     nb_inputs, x0s, NULL);
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

    double* x_solutions = malloc(nb_inputs * nb_solutions * sizeof(double));
    double* output_solutions = malloc(nb_outputs * nb_solutions * sizeof(double));
    loadInputSolutionsNomadResult(x_solutions, nb_solutions, nomad_result);
    loadOutputSolutionsNomadResult(output_solutions, nb_solutions, nomad_result);
    printf("Solutions:\n");
    for (int index = 0; index < nb_solutions; ++index)
    {
        printf("sol %d: x = [", index+1);
        for (int i = 0; i < nb_inputs; ++i)
        {
            printf(" %f", x_solutions[index * nb_inputs + i]);
        }
        printf("]; f = [%f %f %f]\n",
               output_solutions[index * nb_outputs],
               output_solutions[index * nb_outputs + 1],
               output_solutions[index * nb_outputs + 2]);
    }

    free(x_solutions);
    free(output_solutions);
    freeNomadProblem(nomad_pb);
    freeNomadResult(nomad_result);
    return 0;
}

int main()
{
    solve_marine_design_pb();
    return EXIT_SUCCESS;
}
