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
/*  example of a program that makes NOMAD do a local step stop (search)     */
/*  after a user criterion                                                  */
/*--------------------------------------------------------------------------*/
#include "Nomad/nomad.hpp"
#include "Algos/EvcInterface.hpp"
#include "Algos/Mads/MadsMegaIteration.hpp"
#include "Algos/Mads/SearchMethodAlgo.hpp"
#include "Cache/CacheBase.hpp"
#include "Type/EvalSortType.hpp"
#include "Algos/AlgoStopReasons.hpp"
#include "Util/AllStopReasons.hpp"
#include "../../Output/OutputQueue.hpp"

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
    auto bbOutputType = _evalParams->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE");
    std::string bbo = x[4].tostring();
    bbo += " " + constr1.tostring();
    bbo += " " + constr2.tostring();
    x.setBBO(bbo, bbOutputType, getEvalType());

    countEval = true; // count a black-box evaluation

    return true;       // the evaluation succeeded
}


void initAllParams( std::shared_ptr<NOMAD::AllParameters> allParams, const size_t n)
{
    
    // Parameters creation
    allParams->setAttributeValue("DIMENSION", n);
    // 100 black-box evaluations
    allParams->setAttributeValue("MAX_BB_EVAL", 100);
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
    NOMAD::BBOutputTypeList bbOutputTypes;
    bbOutputTypes.push_back(NOMAD::BBOutputType::OBJ);
    bbOutputTypes.push_back(NOMAD::BBOutputType::PB);
    bbOutputTypes.push_back(NOMAD::BBOutputType::PB);
    allParams->setAttributeValue("BB_OUTPUT_TYPE", bbOutputTypes );

    allParams->setAttributeValue("DISPLAY_DEGREE", 3);
    allParams->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("bbe ( sol ) obj"));
    allParams->setAttributeValue("DISPLAY_ALL_EVAL", true);
    
    // Algo
    // allParams->setAttributeValue("QUAD_MODEL_SEARCH", false);
    // allParams->setAttributeValue("NM_SEARCH", true);
    //allParams->setAttributeValue("DIRECTION_TYPE", NOMAD::DirectionType::ORTHO_2N);
    // allParams->setAttributeValue("EVAL_QUEUE_SORT", NOMAD::EvalSortType::DIR_LAST_SUCCESS);

    // Parameters validation
    allParams->checkAndComply();
    

}


/*-----------------------------------------------------*/
/* After each evaluation, verify if                    */
/* remaining evaluations in queue should be cancelled. */
/*-----------------------------------------------------*/
void customEvalCB(NOMAD::EvalQueuePointPtr & evalQueuePoint, bool &opportunisticStop)
{
    if (NOMAD::EvalType::BB == evalQueuePoint->getEvalType() )
    {
        
        auto bbe = NOMAD::EvcInterface::getEvaluatorControl()->getBbEval();
        
        
        auto mystep = evalQueuePoint->getGenStep();
        
        if(bbe>10 && mystep==NOMAD::StepType::SEARCH_METHOD_SPECULATIVE)
        {
            opportunisticStop=true;
            std::cout<<"Custom opportunistic stop in speculative search"<<std::endl; // FIX : here we have a problem : nelder mead is still run after the speculative search
        }
        
        // if (bbe%3==0)
        // {
        //     opportunisticStop = true;
        // }
    }
    
}


void cbPostprocessing(const NOMAD::Step & step, bool &stop)
{
    stop = false;
    auto evc = NOMAD::EvcInterface::getEvaluatorControl();
   
    // This is postprocessing for BB only
    if (NOMAD::EvalType::BB != evc->getCurrentEvalType())
    {
        return;
    }
    
    auto evcStopReason = evc->getStopReason(-1);
    std::cout<<"cbPostprocessing D"<<std::endl;
    if ( evcStopReason.checkStopType(NOMAD::EvalMainThreadStopType::CUSTOM_OPPORTUNISTIC_STOP))
    {
        // NOTE: To output information, it is safer to use "Add" instead of "AddOutputInfo" in callback to avoid segmentation faults
        // when accessing step name as callbacks may be called deep in code, e.g. in EvaluatorControl just after evaluation 
        NOMAD::OutputQueue::Add("Custom opportunistic triggered. We don't want a global stop. Return stop = false. Reset evcStopReason to continue.", NOMAD::OutputLevel::LEVEL_INFO);
        NOMAD::OutputQueue::Flush();


        evc->setStopReason(-1, NOMAD::EvalMainThreadStopType::STARTED);


        // Is this a search ? Then stop its iteration.
        NOMAD::Search * searchStep = step.getParentOfType<NOMAD::Search*>(false);
        if(nullptr!= searchStep)
        {
            std::cout<<"Let's do a user stop of the search " << searchStep->getName() <<std::endl;

            // stop the current search and go the the end of main algorithm megaIteration
            searchStep->getAllStopReasons()->set(NOMAD::IterStopType::USER_ITER_STOP);
            
            // Is this an algo ? Look for an algorithm among the parents. It should be a search algorithm, that is, not the root Algorithm. Stop it completely if found.
            // If we simply search for an algorithm among the parents we may endup with the root Mads algorithm. 
            auto algoSM = step.getFirstAlgorithm();
            if (algoSM != step.getRootAlgorithm())
            {
                std::cout<<"Let's do a user stop of the search algo " << algoSM->getName() <<std::endl;
                algoSM->getAllStopReasons()->set(NOMAD::IterStopType::USER_ALGO_STOP);
            }

        }
    }
    
    std::cout<<"cbPostprocessing F"<<std::endl;
}



/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv )
{
    // Dimension (Number of variables)
    size_t n = 5;
    
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
    // Add callback function run just after evaluation
    NOMAD::EvcInterface::getEvaluatorControl()->addEvalCallback<NOMAD::CallbackType::EVAL_OPPORTUNISTIC_CHECK>(cbInter);
    
    // The run
    TheMainStep.start();
    TheMainStep.run();
    TheMainStep.end();
        
    return 1;
}
