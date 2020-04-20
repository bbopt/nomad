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

#include "Algos/MainStep.hpp"
#include "Eval/Evaluator.hpp"
#include "Param/AllParameters.hpp"

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
class My_Evaluator : public NOMAD::Evaluator
{
public:
    My_Evaluator(const std::shared_ptr<NOMAD::EvalParameters>& evalParams)
        : NOMAD::Evaluator(evalParams, NOMAD::EvalType::BB)
    {}

    ~My_Evaluator() {}

    bool eval_x(NOMAD::EvalPoint &x, const NOMAD::Double& hMax, bool &countEval) const override
    {
        size_t n = x.size();

        NOMAD::Double f = 0.0;  // Objective value
        NOMAD::Double c1 = 0.0; // Constraint 1
        NOMAD::Double c2 = 0.0; // Constraint 2

        try
        {
            for (size_t i = 0; i < n; i++)
            {
                NOMAD::Double xi = x[i];
                c1 += (xi-1).pow2();
                c2 += (xi+1).pow2();
            }
            c1 = c1-25;
            c2 = 25-c2;
            if (c1*c1 + c2*c2 > hMax)
            {
                // Does not count as an evaluation.
                // This can be useful in case the constraints take less time
                // to compute than the function itself.
                countEval = false;
                std::string bbo = f.tostring() + " " + c1.tostring() + " " + c2.tostring();
                x.setBBO(bbo, _evalParams->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE"), getEvalType());
            }
            else
            {
                f = x[n-1];
                std::string bbo = f.tostring() + " " + c1.tostring() + " " + c2.tostring();
                x.setBBO(bbo, _evalParams->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE"), getEvalType());
                countEval = true; // count a black-box evaluation
            }
        }
        catch (std::exception &e)
        {
            std::string err("Exception: ");
            err += e.what();
            throw std::logic_error(err);
        }

        return true;     // the evaluation succeeded
    }
};


void initParams(NOMAD::AllParameters &p)
{
    // parameters creation
    size_t n = 5;   // Number of variables
    p.getPbParams()->setAttributeValue("DIMENSION", n);
    p.getEvalParams()->setAttributeValue("BB_OUTPUT_TYPE", NOMAD::stringToBBOutputTypeList("OBJ PB PB"));

    p.getPbParams()->setAttributeValue("X0", NOMAD::Point(n,0.0));  // starting point (0.0 0.0 0.0 0.0 0.0)
    p.getPbParams()->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(n, -6.0)); // all var. >= -6
    NOMAD::ArrayOfDouble ub(n);     // x_4 and x_5 have no bounds
    ub[0] = 5.0;                    // x_1 <= 5
    ub[1] = 6.0;                    // x_2 <= 6
    ub[2] = 7.0;                    // x_3 <= 7
    p.getPbParams()->setAttributeValue("UPPER_BOUND", ub);
    NOMAD::ArrayOfDouble granularity(n, 10e-2);
    p.getPbParams()->setAttributeValue("GRANULARITY", granularity);

    // the algorithm terminates after 100 black-box evaluations,
    // or 100 total evaluations, including cache hits and evalutions for
    // which countEval was false.
    p.getEvaluatorControlParams()->setAttributeValue("MAX_BB_EVAL", 100);
    p.getEvaluatorControlParams()->setAttributeValue("MAX_EVAL", 200);

    p.getRunParams()->setAttributeValue("H_MAX_0", NOMAD::Double(10000000));

    p.getDispParams()->setAttributeValue("DISPLAY_DEGREE", 2);
    p.getDispParams()->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("EVAL ( SOL ) OBJ CONS_H H_MAX"));

    p.getRunParams()->setAttributeValue("HOT_RESTART_READ_FILES", false);
    p.getRunParams()->setAttributeValue("HOT_RESTART_WRITE_FILES", false);

    // parameters validation
    p.checkAndComply();
}

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main (int argc, char **argv)
{
    auto TheMainStep = std::make_unique<NOMAD::MainStep>();

    // Initialize all parameters
    auto params = std::make_shared<NOMAD::AllParameters>();
    initParams(*params);
    TheMainStep->setAllParameters(params);

    // Custom evaluator creation
    std::unique_ptr<My_Evaluator> ev(new My_Evaluator(params->getEvalParams()));
    TheMainStep->setEvaluator(std::move(ev));

    try
    {
        // Algorithm creation and execution
        TheMainStep->start();
        TheMainStep->run();
        TheMainStep->end();
    }

    catch(std::exception &e)
    {
        std::cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
    }

    return 0;
}
