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

/*--------------------------------------------------------------------------*/
/*  Example of a program that makes NOMAD do a local opportunistic stop     */
/*  of queued evaluations from a Poll step when a user criterion is met     */
/*--------------------------------------------------------------------------*/
#include "Nomad/nomad.hpp"
#include "Algos/EvcInterface.hpp"
#include "Algos/Mads/MadsMegaIteration.hpp"
#include "Algos/Mads/SearchMethodAlgo.hpp"
#include "Cache/CacheBase.hpp"
#include "Type/EvalSortType.hpp"
#include "Algos/AlgoStopReasons.hpp"
#include "Util/AllStopReasons.hpp"
#include "Output/OutputQueue.hpp"

// Global variable to store current best feasible point
// for custom opportunistic stop of step
NOMAD::Double currentBestFeasF;

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
class My_Evaluator : public NOMAD::Evaluator
{
private:

public:
    My_Evaluator(const std::shared_ptr<NOMAD::EvalParameters>& evalParams)
    : NOMAD::Evaluator(evalParams, NOMAD::EvalType::BB)
    {}

    ~My_Evaluator() {}

    bool eval_x(NOMAD::EvalPoint &x, const NOMAD::Double &hMax, bool &countEval) const override;
};


/*----------------------------------------*/
/*           user-defined eval_x          */
/*----------------------------------------*/
bool My_Evaluator::eval_x(NOMAD::EvalPoint &x,
                          const NOMAD::Double &hMax,
                          bool &countEval) const
{
    
    NOMAD::Double f = pow ( 10 * (x[1].todouble() - pow(x[0].todouble(),2) ) , 2 );
    f += pow ( 1 - x[0].todouble() , 2 );
    std::string bbo = f.tostring();
    x.setBBO(bbo);

    countEval = true; // count a black-box evaluation

    return true;       // the evaluation succeeded
}


void initAllParams( std::shared_ptr<NOMAD::AllParameters> allParams, const size_t n)
{

    // Parameters creation
    allParams->setAttributeValue("DIMENSION", n);
    // 100 black-box evaluations
    allParams->setAttributeValue("MAX_BB_EVAL", 400);
    // Starting point
    allParams->setAttributeValue("X0", NOMAD::Point(n, 0.0) );

    // Bounds
    allParams->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(n, -10.0));
    allParams->setAttributeValue("UPPER_BOUND", NOMAD::ArrayOfDouble(n, 10.0));

    // Constraints and objective
    NOMAD::BBOutputTypeList bbOutputTypes;
    bbOutputTypes.push_back(NOMAD::BBOutputType::OBJ);
    allParams->setAttributeValue("BB_OUTPUT_TYPE", bbOutputTypes );

    allParams->setAttributeValue("DISPLAY_DEGREE", 3);
    allParams->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("bbe ( sol ) obj"));
    allParams->setAttributeValue("DISPLAY_ALL_EVAL", true);

    //  Opportunistic eval must be activated (default is activated!). Let us sort with a poor sorting strategy (quadratic model is much better). Default criterion for opportunism is disabled and replaced by a custom opportunistic criterion provided by a user callback (see below).
    allParams->setAttributeValue("EVAL_OPPORTUNISTIC", true);
    allParams->setAttributeValue("EVAL_QUEUE_SORT",NOMAD::EvalSortType::LEXICOGRAPHICAL);

    allParams->setAttributeValue("DIRECTION_TYPE",NOMAD::DirectionType::ORTHO_2N);
    allParams->setAttributeValue("QUAD_MODEL_SEARCH", false);
    allParams->setAttributeValue("NM_SEARCH", false);
    allParams->setAttributeValue("SPECULATIVE_SEARCH", false);

    // Parameters validation
    allParams->checkAndComply();


}


/*-----------------------------------------------------*/
/* After each evaluation, verify if                    */
/* remaining evaluations in queue should be cancelled. */
/*-----------------------------------------------------*/
void customEvalCB(NOMAD::EvalQueuePointPtr & evalQueuePoint, bool &opportunisticEvalStop, bool &opportunisticIterStop)
{
    // Consider only BB evals
    if (NOMAD::EvalType::BB == evalQueuePoint->getEvalType() )
    {
        // Consider only feasible points
        if (evalQueuePoint->isFeasible(NOMAD::EvalType::BB))
    {
            // Update my current best feasible point
            if (!currentBestFeasF.isDefined())
            {
                currentBestFeasF = evalQueuePoint->getF(NOMAD::EvalType::BB);
                return;
            }

            auto mystep = evalQueuePoint->getGenStep();

            // Opportunism only if enough reduction is obtained (optim f is 0)
            auto FMinOpport = currentBestFeasF - 0.05*currentBestFeasF.abs();
            if (NOMAD::stepTypeToString(mystep).find("Poll") != string::npos &&
                evalQueuePoint->getF(NOMAD::EvalType::BB) < FMinOpport)
            {
                opportunisticEvalStop=true;
                std::cout<<"*****************************************************"<< std::endl;
                std::cout<<"Opportunistic stop in Poll on f sufficient decrease. "<< std::endl;
                std::cout<<"*****************************************************"<< std::endl;
                
                currentBestFeasF = evalQueuePoint->getF(NOMAD::EvalType::BB);
            }
        }
    }
}



/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv )
{
    // Dimension (Number of variables)
    // size_t n = 5;
    size_t n = 2;

    NOMAD::MainStep TheMainStep;

    // Set parameters
    auto params = std::make_shared<NOMAD::AllParameters>();
    initAllParams(params, n);
    TheMainStep.setAllParameters(params);

    // Custom Evaluator
    std::unique_ptr<My_Evaluator> ev(new My_Evaluator(params->getEvalParams()));
    TheMainStep.setEvaluator(std::move(ev));

    // Link callback function with user function defined locally
    NOMAD::EvalCallbackFunc<NOMAD::CallbackType::EVAL_OPPORTUNISTIC_CHECK> cbInter = customEvalCB;
    //  Add callback function run just after evaluation
    NOMAD::EvcInterface::getEvaluatorControl()->addEvalCallback<NOMAD::CallbackType::EVAL_OPPORTUNISTIC_CHECK>(cbInter);

    // The run
    TheMainStep.start();
    TheMainStep.run();
    TheMainStep.end();

    return 1;
}
