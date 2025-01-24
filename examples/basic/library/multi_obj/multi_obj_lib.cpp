/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4 is owned by                                 */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD 4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada,             */
/*  NSERC (Natural Sciences and Engineering Research Council of Canada),           */
/*  InnovÉÉ (Innovation en Énergie Électrique) and IVADO (The Institute            */
/*  for Data Valorization)                                                         */
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
 \file   multi_obj_lib.cpp
 \brief  Library example for nomad (DMulti-MADS algorithm)
 \author Ludovic Salomon
 \date   2024
 */

#include "Nomad/nomad.hpp"
#include "Type/DMultiMadsSearchStrategyType.hpp"

// A conceptual marine design problem
// See "An Easy-To-use Real-world Multi-objective Optimization Suite problem"
// by R. Tanabe and H. Ishibuchi
class MarineDesignProblem : public NOMAD::Evaluator
{
public:
    explicit MarineDesignProblem(const std::shared_ptr<NOMAD::EvalParameters>& evalParams)
    : NOMAD::Evaluator(evalParams, NOMAD::EvalType::BB)
    {}

    ~MarineDesignProblem() override = default;

    bool eval_x(NOMAD::EvalPoint& x, const NOMAD::Double& hMax, bool& countEval) const override
    {
        bool eval_ok = false;

        NOMAD::Double f1 = 1e20, f2 = 1e20, f3 = 1e20;

        try
        {
            // Variables
            const double L = x[0].todouble();
            const double B = x[1].todouble();
            const double D = x[2].todouble();
            const double T = x[3].todouble();
            const double Vk = x[4].todouble();
            const double Cb = x[5].todouble();

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
            const double froudenumber = V / std::pow(g * L, 0.5);
            const double P = std::pow(displacement, 2.0/3) * std::pow(Vk, 3) / (a + b * froudenumber);

            const double steelWeight = 0.034 * std::pow(L, 1.7) * std::pow(B, 0.7) *
                                       std::pow(D, 0.4) * std::pow(Cb, 0.5);
            const double outfitWeight = 1.0 * std::pow(L, 0.8) * std::pow(B, 0.6) *
                                       std::pow(D, 0.3) * std::pow(Cb, 0.1);
            const double machineryWeight = 0.17 * std::pow(P, 0.9);
            const double lightShipWeight = steelWeight + outfitWeight + machineryWeight;

            const double deadweight = displacement - lightShipWeight;

            const double dailyConsumption = (0.19 * P * 24) / 1000 + 0.2;
            const double seaDays = (roundTripMiles * Vk) / 24;
            const double fuelCarried = dailyConsumption * (seaDays + 5);
            const double portCost = 6.3 * std::pow(deadweight, 0.8);

            const double cargoDeadweight = deadweight - fuelCarried - 2 * std::pow(deadweight, 0.5);
            const double portDays = 2 * (cargoDeadweight / handlingRate + 0.5);
            const double roundTripsPerYear = 350.0 / (seaDays + portDays);

            const double fuelCost = 1.05 * dailyConsumption * seaDays * fuelPrice;
            const double runningCosts = 40000 * std::pow(deadweight, 0.8);
            const double voyageCosts = (fuelCost + portCost) * roundTripsPerYear;
            const double shipCosts = 1.3 * (2000 * std::pow(steelWeight, 0.85) +
                                     3500 * outfitWeight + 2400 * std::pow(P, 0.8));
            const double capitalCosts = 0.2 * shipCosts;
            const double annualCosts = capitalCosts + runningCosts + voyageCosts;

            const double annualCargo = - cargoDeadweight * roundTripsPerYear;

            // Objectives
            f1 = annualCosts / annualCargo;
            f2 = lightShipWeight;
            f3 = annualCargo;

            // Constraints
            const NOMAD::Double g1 = L / B - 6.0; // >= 0
            const NOMAD::Double g2 = 15.0 - L / D; // >= 0
            const NOMAD::Double g3 = 19.0 - L / T; // >= 0
            const NOMAD::Double g4 = 0.45 * std::pow(deadweight, 0.31) - T; // >= 0
            const NOMAD::Double g5 = 0.7 * D + 0.7 - T; // >= 0
            const NOMAD::Double g6 = deadweight - 3000; // >= 0
            const NOMAD::Double g7 = 500000 - deadweight; // >= 0
            const NOMAD::Double g8 = 0.32 - froudenumber; // >= 0
            const NOMAD::Double g9 = verticalCenterOfBuoyancy + metacentricRadius -
                                     verticalCenterOfGravity - 0.07 * B; // >= 0

            std::string bbo = f1.tostring() + " " + f2.tostring() + " " + f3.tostring();
            bbo += " " + (-g1).tostring();
            bbo += " " + (-g2).tostring();
            bbo += " " + (-g3).tostring();
            bbo += " " + (-g4).tostring();
            bbo += " " + (-g5).tostring();
            bbo += " " + (-g6).tostring();
            bbo += " " + (-g7).tostring();
            bbo += " " + (-g8).tostring();
            bbo += " " + (-g9).tostring();

            x.setBBO(bbo);

            eval_ok = true;
        }
        catch (std::exception& e)
        {
            std::string err("Exception: ");
            err += e.what();
            throw std::logic_error(err);
        }
        countEval = true;
        return eval_ok;
    }
};

// Main function
int main()
{
    NOMAD::MainStep TheMainStep;

    try
    {
        // Parameters creation
        auto params = std::make_shared<NOMAD::AllParameters>();

        // Dimensions of the blackbox, inputs and outputs
        const size_t n = 6;
        params->setAttributeValue("DIMENSION", n);

        NOMAD::ArrayOfDouble lb(n, 0);
        lb[0] = 150;
        lb[1] = 20;
        lb[2] = 13;
        lb[3] = 10;
        lb[4] = 14;
        lb[5] = 0.63;
        params->setAttributeValue("LOWER_BOUND", lb);

        NOMAD::ArrayOfDouble ub(n, 0);
        ub[0] = 274.32;
        ub[1] = 32.31;
        ub[2] = 25;
        ub[3] = 11.71;
        ub[4] = 18;
        ub[5] = 0.75;
        params->setAttributeValue("UPPER_BOUND", ub);

        // Outputs, constraints and objectives
        NOMAD::BBOutputTypeList bbOutputTypes;
        for (size_t i = 0; i < 3; ++i)
        {
            bbOutputTypes.emplace_back(NOMAD::BBOutputType::OBJ);
        }
        for (size_t i = 0; i < 9; ++i)
        {
            bbOutputTypes.emplace_back(NOMAD::BBOutputType::EB); // We want a two-phase approach
        }
        params->setAttributeValue("BB_OUTPUT_TYPE", bbOutputTypes);

        // A line initialization is practically more efficient than giving a single point for
        // multiobjective optimization.
        // The interesting reader can report to the following reference for more information.
        // Direct Multisearch for multiobjective optimization
        // by A.L. Custodio, J.F.A. Madeira, A.I.F. Vaz and L.N. Vicente, 2011.
        NOMAD::ArrayOfPoint x0s;
        for (size_t j = 0; j < n; ++j)
        {
            NOMAD::Point x0(n, 0);
            for (size_t i = 0; i < n; ++i)
            {
                x0[i] = lb[i] + (double) j * (ub[i] - lb[i]) / (n - 1);
            }
            x0s.push_back(x0);
        }
        // Starting point
        params->setAttributeValue("X0", x0s);

        // Algorithm parameters
        // 1- Terminate after this number of maximum blackbox evaluations.
        params->setAttributeValue("MAX_BB_EVAL", 3000);

        // 2- Use n+1 directions
        params->setAttributeValue("DIRECTION_TYPE", NOMAD::DirectionType::ORTHO_NP1_NEG);

        // 3- For multiobjective optimization, these parameters are required.
        params->setAttributeValue("DMULTIMADS_OPTIMIZATION", true);
        // For multiobjective optimization, sort cannot use the default quad model info.
        params->setAttributeValue("EVAL_QUEUE_SORT", NOMAD::EvalSortType::DIR_LAST_SUCCESS);
        params->setAttributeValue("NM_SEARCH", false); // Deactivate Nelder-Mead search
        // NB: by default, QUAD_MODEL_SEARCH is activated; to deactivate it, uncomment
        // params->setAttributeValue("QUAD_MODEL_SEARCH", true); // Deactivate Quad Model search
        // Change Quad Model search strategy for DMultiMads
        params->setAttributeValue("DMULTIMADS_QUAD_MODEL_STRATEGY", NOMAD::DMultiMadsQuadSearchType::DMS);

        // Advanced attributes for DMultiMads
        params->setAttributeValue("DMULTIMADS_SELECT_INCUMBENT_THRESHOLD", 2);

        // 4- Other useful parameters
        params->setAttributeValue("DISPLAY_DEGREE", 2);
        params->setAttributeValue("DISPLAY_ALL_EVAL",true);
        params->setAttributeValue("SOLUTION_FILE", std::string("sol.txt")); // Save the Pareto front approximation
        // params->setAttributeValue("HISTORY_FILE", std::string("history.txt")); // Save history. To uncomment if you want it.

        // Validate
        params->checkAndComply();

        // Run the solver
        TheMainStep.setAllParameters(params);
        auto ev = std::make_unique<MarineDesignProblem>(params->getEvalParams());
        TheMainStep.addEvaluator(std::move(ev));

        TheMainStep.start();
        TheMainStep.run();
        TheMainStep.end();
    }

    catch(std::exception &e)
    {
        std::cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
    }
    return 0;
}
