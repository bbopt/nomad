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
/*  example of a program that makes NOMAD do a local stop from user output  */
/*  The user criterion is based upon extra output                           */
/*  Feature is similar to NOMAD 3 STAT_SUM and STAT_SUM_TARGET              */
/*--------------------------------------------------------------------------*/
#include "Nomad/nomad.hpp"
#include "Algos/EvcInterface.hpp"
#include "Algos/Mads/MadsMegaIteration.hpp"
#include "Algos/Mads/SearchMethodAlgo.hpp"
#include "Cache/CacheBase.hpp"
#include "Type/EvalSortType.hpp"
#include "Algos/AlgoStopReasons.hpp"
#include "Util/AllStopReasons.hpp"


// User stats used for stoping optimization
/* This feature is implemented in Nomad 3 as STAT_SUM and STAT_SUM_TARGET
   Nomad 4 has more flexibility by using user evaluation callbacks
 */
NOMAD::Double user_stat_sum = 0;
NOMAD::Double user_stat_target = 10;

NOMAD::BBOutputTypeList bbOutputTypes;

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
    NOMAD::Double c1 = 0.0 , c2 = 0.0;
    for ( int i = 0 ; i < 5 ; i++ )
    {

        c1 += (x[i]-1).pow2();
        c2 += (x[i]+1).pow2();
    }
    NOMAD::Double constr1 = c1-25;
    NOMAD::Double constr2 = 25-c2;
    std::string bbo = x[4].tostring();
    bbo += " " + constr1.tostring();
    bbo += " " + constr2.tostring();
    
    NOMAD::Double bboExtraFeasible_1 =0, bboExtraFeasible_2 = 0;
    if (constr1 <= -5)
        bboExtraFeasible_1 = 1 ;
    if (constr2 <= -5)
        bboExtraFeasible_2 = 1;
    bbo += " " + bboExtraFeasible_1.tostring();
    bbo += " " + bboExtraFeasible_2.tostring();
    x.setBBO(bbo);

    countEval = true; // count a black-box evaluation

    return true;       // the evaluation succeeded
}

void initAllParams(std::shared_ptr<NOMAD::AllParameters> allParams)
{

    const size_t n = 5;
    
    // Parameters creation
    allParams->setAttributeValue("DIMENSION", n);
    
    // Starting point
    allParams->setAttributeValue("X0", NOMAD::Point(n, 0.0) );

    // Bounds
    allParams->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(n, -6.0 )); // all var. >= -6
    NOMAD::ArrayOfDouble ub(n);
    ub[0] = 5.0;    // x_1 <= 5
    ub[1] = 6.0;    // x_2 <= 6
    ub[2] = 7.0;    // x_3 <= 7
    allParams->setAttributeValue("UPPER_BOUND", ub);

    // Constraints and objective
    bbOutputTypes.push_back(NOMAD::BBOutputType::Type::OBJ);
    bbOutputTypes.push_back(NOMAD::BBOutputType::Type::EB);
    bbOutputTypes.push_back(NOMAD::BBOutputType::Type::EB);
    bbOutputTypes.push_back(NOMAD::BBOutputType::Type::BBO_UNDEFINED);
    bbOutputTypes.push_back(NOMAD::BBOutputType::Type::BBO_UNDEFINED);
    allParams->setAttributeValue("BB_OUTPUT_TYPE", bbOutputTypes );

    allParams->setAttributeValue("DISPLAY_DEGREE", 2);
    allParams->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("bbe ( sol ) bbo"));
    allParams->setAttributeValue("DISPLAY_ALL_EVAL", false);
    
    allParams->setAttributeValue("MAX_BB_EVAL",100);
    
    // Parameters validation requested to have access to their value.
    allParams->checkAndComply();

}


/*------------------------------------------------------------------------*/
/* After each evaluation compute user stats on extra outputs and check on */
/* user stat stopping criterion                                           */
/*------------------------------------------------------------------------*/
void customEvalStopCB( NOMAD::EvalQueuePointPtr & evalQueuePoint, bool & globalStop)
{
    globalStop = false;
    if (nullptr != evalQueuePoint)
    {
        auto eval = evalQueuePoint->getEval(NOMAD::EvalType::BB);
        if (nullptr != eval)
        {
            auto extra_bbo = eval->getBBOutput().getExtraOutputs(bbOutputTypes) ;
            // Dummy computation of user stat
            if (extra_bbo[0] == 1 && extra_bbo[1] == 1)
            {
                user_stat_sum++;
            }
            // Dummy criterion for user stat stopping
            if ( user_stat_sum == user_stat_target)
            {
                globalStop = true;
            }
        }
    }
}



/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv )
{

    
    NOMAD::MainStep TheMainStep;
        
    // Set parameters
    auto params = std::make_shared<NOMAD::AllParameters>();
    initAllParams(params);
    TheMainStep.setAllParameters(params);
    
    // Custom Evaluator creation
    auto ev = std::make_unique<My_Evaluator>(params->getEvalParams());
    TheMainStep.addEvaluator(std::move(ev));
    
    // Link callback function with user function defined locally
    NOMAD::EvalCallbackFunc<NOMAD::CallbackType::EVAL_STOP_CHECK> cbFailCheck = customEvalStopCB;
    
    // Add callback function for eval check and stop management.
    NOMAD::EvcInterface::getEvaluatorControl()->addEvalCallback<NOMAD::CallbackType::EVAL_STOP_CHECK>(cbFailCheck);

    
    // The run
    TheMainStep.start();
    TheMainStep.run();
    TheMainStep.end();
        
    return 1;
}
