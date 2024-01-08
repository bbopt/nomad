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

//
// IMPORTANT
//
// PB is Styrene run in hybrid library/batch mode.
// The styrene standalone executable is registered with the function call
// setAttributeValue("BB_EXE",xxxx)
// To run this optimization, the program must be executed in a path where styrene truth executable is available.
// Styrene sources are available at https://github.com/bbopt/styrene and must be compiled prior to run this optimization.

void initAllParams( std::shared_ptr<NOMAD::AllParameters> allParams)
{
    const int n = 8;
    
    // Parameters creation
    allParams->setAttributeValue("DIMENSION", n);
    // 100 black-box evaluations
    allParams->setAttributeValue("MAX_BB_EVAL", 100);
    // Starting point
    std::vector<double> x0 = { 54, 66, 86, 8, 29, 51, 32, 15};
    allParams->setAttributeValue("X0", NOMAD::Point(x0) );
    allParams->setAttributeValue("BB_EXE", std::string("./truth.exe"));
    

    // Bounds
    allParams->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(n, 0.0 ));
    allParams->setAttributeValue("UPPER_BOUND", NOMAD::ArrayOfDouble(n, 100.0 ));

    // Constraints and objective
    NOMAD::BBOutputTypeList bbOutputTypes= {NOMAD::BBOutputType::EB, NOMAD::BBOutputType::EB, NOMAD::BBOutputType::EB, NOMAD::BBOutputType::EB, NOMAD::BBOutputType::PB, NOMAD::BBOutputType::PB, NOMAD::BBOutputType::PB, NOMAD::BBOutputType::PB, NOMAD::BBOutputType::PB, NOMAD::BBOutputType::PB, NOMAD::BBOutputType::PB, NOMAD::BBOutputType::OBJ};
    allParams->setAttributeValue("BB_OUTPUT_TYPE", bbOutputTypes );

    allParams->setAttributeValue("DISPLAY_DEGREE", 2);
    allParams->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("bbe ( sol ) obj"));


    // Parameters validation
    allParams->checkAndComply();
    

}


/*------------------------------------------------------------------------*/
/* After failed evaluation, may apply a special treatment to eval point   */
/*------------------------------------------------------------------------*/
void customFailEvalCB( NOMAD::EvalQueuePointPtr & evalQueuePoint)
{
    if (nullptr != evalQueuePoint)
    {
    auto eval = evalQueuePoint->getEval(NOMAD::EvalType::BB);
    if ( nullptr != eval )
    {
        if (! eval->isBBOutputComplete())
        {
            eval->setBBO(std::string("1 1 1 1 1 1 1 1 1 1 1 1"), eval->getBBOutputTypeList(), true);
            std::cout<<"BBoutput was incomplete. BBO is reset to be complete."<<std::endl;
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
    
    // Link callback function with user function defined locally
    NOMAD::EvalCallbackFunc<NOMAD::CallbackType::EVAL_FAIL_CHECK> cbFailCheck = customFailEvalCB;
    
    // Add callback function for fail eval check and management.
    NOMAD::EvcInterface::getEvaluatorControl()->addEvalCallback<NOMAD::CallbackType::EVAL_FAIL_CHECK>(cbFailCheck);

    
    // The run
    TheMainStep.start();
    TheMainStep.run();
    TheMainStep.end();
        
    return 1;
}
