/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created by                                          */
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

#include "Nomad/nomad.hpp"

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
class My_Evaluator : public NOMAD::Evaluator
{
public:
    My_Evaluator(const std::shared_ptr<NOMAD::EvalParameters>& evalParams, NOMAD::EvalType evalType)
        : NOMAD::Evaluator(evalParams, evalType)
    {}

    ~My_Evaluator() {}

    bool eval_x(NOMAD::EvalPoint &x, const NOMAD::Double& hMax, bool &countEval) const override
    {
        bool eval_ok = false;
        size_t n = x.size();

        NOMAD::Double f = 0.0;  // Objective value
        NOMAD::Double c1 = 0.0; // Constraint 1
        NOMAD::Double c2 = 0.0; // Constraint 2

        try
        {
            if (NOMAD::EvalType::BB == _evalType)
            {
                for (size_t i = 0; i < n; i++)
                {
                    NOMAD::Double xi = x[i];
                    c1 += (xi-1).pow2();
                    c2 += (xi+1).pow2();
                }
                c1 = c1-25;
                c2 = 25-c2;

            }
            else
            {
                c1 = x[0]*x[0]-5;
                c2 = 5-x[1]*x[1];
            }

            f = x[n-1];
            std::string bbo = f.tostring() + " " + c1.tostring() + " " + c2.tostring();
            x.setBBO(bbo); // Eval type and Bb output types are determined automatically
            eval_ok = true;
        }
        catch (std::exception &e)
        {
            std::string err("Exception: ");
            err += e.what();
            throw std::logic_error(err);
        }


        countEval = true;  // count a black-box evaluation

        return eval_ok;     // the evaluation succeeded
    }

    // Wrapper around eval_x
    // If eval_block is not defined here, eval_x is called sequentially
    // for each point in the block, and a warning is shown.
    // The user may redefine eval_block to optimize parallelism management.
    std::vector<bool> eval_block(std::vector<std::shared_ptr<NOMAD::EvalPoint>> &block,
                                 const NOMAD::Double& hMax,
                                 std::vector<bool> &countEval) const override
    {
        std::vector<bool> evalOk(block.size(), false);
        countEval.resize(block.size(), false);

        for (size_t index = 0; index < block.size(); index++)
        {
            bool countEval1 = false;
            evalOk[index] = eval_x(*block[index], hMax, countEval1);
            countEval[index] = countEval1;
        }

        return evalOk;
    }
};




void initParams(NOMAD::AllParameters &p)
{
    // parameters creation
    size_t n = 6;   // Number of variables
    p.getPbParams()->setAttributeValue("DIMENSION", n);
    p.getEvalParams()->setAttributeValue("BB_OUTPUT_TYPE", NOMAD::stringToBBOutputTypeList("OBJ EB EB"));

    
    NOMAD::Point X0(n, 0.0);
    // X0[n-1] = -4.0; // starting point (0.0 0.0 0.0 0.0 0.0 -4.0)
    p.getPbParams()->setAttributeValue("X0", X0);
    p.getPbParams()->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(n, -6.0)); // all var. >= -6
    NOMAD::ArrayOfDouble ub(n);     // x_4 and x_5 have no bounds
    ub[0] = 5.0;                    // x_1 <= 5
    ub[1] = 6.0;                    // x_2 <= 6
    ub[2] = 7.0;                    // x_3 <= 7
    ub[n-1] = 6.0;                  // x_6 <= 6
    p.getPbParams()->setAttributeValue("UPPER_BOUND", ub);

    // the algorithm terminates after MAX_BB_EVAL black-box evaluations.
    p.getEvaluatorControlGlobalParams()->setAttributeValue("MAX_BB_EVAL", 1000);
    p.getEvaluatorControlGlobalParams()->setAttributeValue("EVAL_SURROGATE_COST", 10);
    p.getEvaluatorControlParams()->setAttributeValue("EVAL_QUEUE_SORT",NOMAD::EvalSortType::SURROGATE);
    

    // Using surrogate for sort and the options below => no quad model is used.
    p.getRunParams()->setAttributeValue("QUAD_MODEL_SEARCH", false);
    p.getRunParams()->setAttributeValue("NM_SEARCH", false);
    p.getRunParams()->setAttributeValue("DIRECTION_TYPE", NOMAD::DirectionType::ORTHO_NP1_NEG);

    p.getDispParams()->setAttributeValue("DISPLAY_DEGREE", 2);


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

    // Custom BB evaluator creation
    auto evBB = std::make_unique<My_Evaluator>(params->getEvalParams(),NOMAD::EvalType::BB);
    TheMainStep->addEvaluator(std::move(evBB));

    // Custom SURROGATE evaluator creation
    auto evSurrogate = std::make_unique<My_Evaluator>(params->getEvalParams(),NOMAD::EvalType::SURROGATE);
    TheMainStep->addEvaluator(std::move(evSurrogate));
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
