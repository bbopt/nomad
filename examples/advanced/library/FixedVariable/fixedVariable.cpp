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

#include "Algos/EvcInterface.hpp"
#include "Algos/MainStep.hpp"
#include "Cache/CacheBase.hpp"
#include "Eval/Evaluator.hpp"
#include "Param/AllParameters.hpp"



void initParams1(NOMAD::AllParameters &p)
{
    // parameters creation
    size_t n = 5;   // Number of variables
    p.getPbParams()->setAttributeValue("DIMENSION", n);
    p.getEvalParams()->setAttributeValue("BB_EXE", std::string("./u.exe"));

    NOMAD::Point x0(n, 6.0);
    x0[0] = 4;
    x0[4] = 5;
    p.getPbParams()->setAttributeValue("X0", x0);  // starting point ( 4 6 6 6 5 )
    NOMAD::ArrayOfDouble lowerBound(n);
    lowerBound[0] = -10.0;
    lowerBound[4] = -5.0;
    p.getPbParams()->setAttributeValue("LOWER_BOUND", lowerBound); // ( -10.0 - - - -5.0 )
    NOMAD::ArrayOfDouble upperBound(n, 31.0);
    p.getPbParams()->setAttributeValue("UPPER_BOUND", upperBound);  // * 31.0
    NOMAD::ArrayOfDouble granularity(n, 0.01);
    p.getPbParams()->setAttributeValue("GRANULARITY", granularity);

    p.getEvaluatorControlGlobalParams()->setAttributeValue("MAX_EVAL", size_t(1000));
    p.getEvaluatorControlGlobalParams()->setAttributeValue("MAX_BB_EVAL", size_t(1000));

    p.getDispParams()->setAttributeValue("DISPLAY_DEGREE", 2);
    p.getDispParams()->setAttributeValue("DISPLAY_ALL_EVAL", false);
    p.getDispParams()->setAttributeValue("DISPLAY_INFEASIBLE", true);
    p.getDispParams()->setAttributeValue("DISPLAY_UNSUCCESSFUL", false);
    p.getDispParams()->setAttributeValue("MAX_DISPLAY_STEP_LEVEL", 20);
    p.getDispParams()->setAttributeValue("DISPLAY_STATS", NOMAD::ArrayOfString("BBE THREAD_NUM ( SOL ) OBJ"));

    p.getEvaluatorControlGlobalParams()->setAttributeValue("TMP_DIR", std::string("/tmp"));
    p.getEvalParams()->setAttributeValue("BB_OUTPUT_TYPE", NOMAD::stringToBBOutputTypeList("PB OBJ"));

    p.getRunParams()->setAttributeValue("ADD_SEED_TO_FILE_NAMES", false);
    p.getEvaluatorControlParams()->setAttributeValue("OPPORTUNISTIC_EVAL", true);
    p.getRunParams()->setAttributeValue("H_MAX_0", NOMAD::Double(200000));

    // This does not work, to fix.
    //p.getPbParams()->setAttributeValue("FIXED_VARIABLE", "0-2");
    NOMAD::Point fixedVariable(n);
    fixedVariable[0] = x0[0];
    fixedVariable[1] = x0[1];
    fixedVariable[2] = x0[2];
    p.getPbParams()->setAttributeValue("FIXED_VARIABLE", fixedVariable);

    p.getRunParams()->setAttributeValue("HOT_RESTART_ON_USER_INTERRUPT", false);
    p.getRunParams()->setAttributeValue("HOT_RESTART_READ_FILES", false);
    p.getRunParams()->setAttributeValue("HOT_RESTART_WRITE_FILES", false);

    // This problem does not work well with NM_SEARCH.
    p.getRunParams()->setAttributeValue("NM_SEARCH", false);

    p.getCacheParams()->setAttributeValue("CACHE_FILE", std::string("cache.txt"));
    p.getRunParams()->setAttributeValue("REJECT_UNKNOWN_PARAMETERS", false);


    // parameters validation
    p.checkAndComply();
}

void initParams2(NOMAD::AllParameters &p, const NOMAD::Point& x0)
{
    size_t n = p.getPbParams()->getAttributeValue<size_t>("DIMENSION");

    p.getPbParams()->setAttributeValue("X0", x0);   // starting point provided
                                                    // Now, we can also simply unset X0.
                                                    // Best point from cache will be used.

    NOMAD::Point fixedVariable(n);
    fixedVariable[2] = x0[2];
    fixedVariable[3] = x0[3];
    p.getPbParams()->setAttributeValue("FIXED_VARIABLE", fixedVariable);

    //p.getEvaluatorControlGlobalParams()->setAttributeValue("MAX_EVAL", 2000);
    //p.getEvaluatorControlGlobalParams()->setAttributeValue("MAX_BB_EVAL", 2000);

    // parameters validation
    p.checkAndComply();
}

void initParams3(NOMAD::AllParameters &p, const NOMAD::Point& x0)
{
    size_t n = p.getPbParams()->getAttributeValue<size_t>("DIMENSION");
    p.getPbParams()->setAttributeValue("X0", x0);

    // Fixed variable values will be computed when X0 is found -
    // the call to checkAndComply() will take care of it.
    NOMAD::Point fixedVariable(n);
    fixedVariable[3].setToBeDefined();
    fixedVariable[4].setToBeDefined();
    p.getPbParams()->setAttributeValue("FIXED_VARIABLE", fixedVariable);

    // parameters validation
    p.checkAndComply();
}

void initParamsFinal(NOMAD::AllParameters &p, const NOMAD::Point& x0)
{
    size_t n = p.getPbParams()->getAttributeValue<size_t>("DIMENSION");
    p.getPbParams()->setAttributeValue("X0", x0);

    // Reset fixed variable.
    NOMAD::Point fixedVariable(n);
    p.getPbParams()->setAttributeValue("FIXED_VARIABLE", fixedVariable);

    // parameters validation
    p.checkAndComply();
}

/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main (int argc, char **argv)
{
    auto TheMainStep = std::make_unique<NOMAD::MainStep>();

    // Initialize parameters
    // Part 1: FIXED_VARIABLE 0-2
    auto params = std::make_shared<NOMAD::AllParameters>();
    initParams1(*params);
    TheMainStep->setAllParameters(params);

    try
    {
        // Algorithm creation and execution
        TheMainStep->start();
        // Here, clear the cache, ensure we start the test fresh.
        NOMAD::CacheBase::getInstance()->clear();
        TheMainStep->run();
        TheMainStep->end();
    }

    catch(std::exception &e)
    {
        std::cerr << "\nFirst run has been interrupted (" << e.what() << ")\n\n";
    }

    // Part 2: FIXED_VARIABLE 2-3
    // Use best point from cache to set X0.
    NOMAD::CacheBase::getInstance()->resetNbCacheHits();
    NOMAD::EvcInterface::getEvaluatorControl()->setNbEval(0);
    std::vector<NOMAD::EvalPoint> bestFeasList;
    NOMAD::CacheBase::getInstance()->findBestFeas(bestFeasList, NOMAD::Point(),
                                                  NOMAD::EvalType::BB, nullptr);
    // NB. Assuming the list is non-empty.
    NOMAD::Point x02 = *(bestFeasList[0].getX());
    initParams2(*params, x02);
    try
    {
        TheMainStep->start();
        TheMainStep->run();
        TheMainStep->end();
    }
    catch(std::exception &e)
    {
        std::cerr << "\nSecond run has been interrupted (" << e.what() << ")\n\n";
    }

    // Part 3: FIXED_VARIABLE 3-4
    NOMAD::CacheBase::getInstance()->resetNbCacheHits();
    NOMAD::EvcInterface::getEvaluatorControl()->setNbEval(0);
    NOMAD::CacheBase::getInstance()->findBestFeas(bestFeasList, NOMAD::Point(),
                                                  NOMAD::EvalType::BB, nullptr);
    NOMAD::Point x03 = *(bestFeasList[0].getX());
    initParams3(*params , x03);
    try
    {
        TheMainStep->start();
        TheMainStep->run();
        TheMainStep->end();
    }
    catch(std::exception &e)
    {
        std::cerr << "\nThird run has been interrupted (" << e.what() << ")\n\n";
    }

    // Final part: No fixed variable
    NOMAD::CacheBase::getInstance()->findBestFeas(bestFeasList, NOMAD::Point(),
                                                  NOMAD::EvalType::BB, nullptr);
    NOMAD::Point x0final = *(bestFeasList[0].getX());
    initParamsFinal(*params,x0final);
    try
    {
        TheMainStep->start();
        TheMainStep->run();
        TheMainStep->end();
    }
    catch(std::exception &e)
    {
        std::cerr << "\nFinal run has been interrupted (" << e.what() << ")\n\n";
    }


    return 0;
}
