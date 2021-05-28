/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0 has been created by                                        */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0 is owned by                               */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada,            */
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
/*--------------------------------------------------------*/
/*  how to use the NOMAD library to suggest trial points  */
/*--------------------------------------------------------*/
#include "Nomad/nomad.hpp"


bool eval_xs(const NOMAD::ArrayOfPoint &xs, std::vector<NOMAD::ArrayOfDouble>& fxs)
{
    bool eval_ok = false;

    size_t nPoints = xs.size();

    NOMAD::ArrayOfDouble output(1);

    try
    {
        std::cout << "Evaluation of points" << std::endl;
        for (size_t i = 0; i < nPoints ; i++)
        {

            // Rosenbrock 2D
            if (xs[i].size() != 2)
            {
                throw NOMAD::Exception(__FILE__,__LINE__,"Dimension should be 2");
            }
            output[0] = pow ( 10 * (xs[i][1].todouble() - pow(xs[i][0].todouble(),2) ) , 2 );
            output[0] += pow ( 1 - xs[i][0].todouble() , 2 );
            fxs.push_back(output);
            std::cout << xs[i] << " -> " << output << std::endl;
        }

        eval_ok = true;
    }
    catch (std::exception &e)
    {
        std::string err("Exception: ");
        err += e.what();
        throw std::logic_error(err);
    }

    return eval_ok;
}

void initParams(std::shared_ptr<NOMAD::AllParameters>& params, const std::string & cacheFileName, const std::string & initialFrameSizeAsString, const std::string & hmax0AsString, bool useCacheAndMegaSearchPoll )
{
    // Problem parameters
    params->setAttributeValue("DIMENSION", 2);             // number of variables
    params->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(2, -10.0));
    params->setAttributeValue("UPPER_BOUND", NOMAD::ArrayOfDouble(2, 10.0));

    NOMAD::BBOutputTypeList bbot;   // Definition of output types
    bbot.push_back(NOMAD::BBOutputType::OBJ);
    params->setAttributeValue("BB_OUTPUT_TYPE", bbot);

    if (! useCacheAndMegaSearchPoll)
    {
        // Use a cache for suggest (no lh).
        
        params->setAttributeValue("CACHE_FILE", cacheFileName);

        // Algorithm parameters
        params->setAttributeValue("MEGA_SEARCH_POLL",true); // Enable Mads MegaSearchPoll to perform suggest.
    }
    else
    {
        // For the first suggest, use lh_eval (may use a cache (if available))
        
        // Algorithm parameters
        params->setAttributeValue("LH_EVAL", 5);
    }
        
    params->readParamLine(initialFrameSizeAsString);
    params->readParamLine(hmax0AsString);


    // Display parameters
    //params->setAttributeValue("DISPLAY_DEGREE", 4);

    // parameters validation
    params->checkAndComply();
}




/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main(int argc, char ** argv)
{
    try
    {
        // Use separate MainStep for Suggest and Observe
        auto SuggestMainStep = std::make_unique<NOMAD::MainStep>();
        auto ObserveMainStep = std::make_unique<NOMAD::MainStep>();

        const size_t maxIteration = 10;

        std::string initialFrameSizeAsString = "INITIAL_FRAME_SIZE ( 0.1 0.1)";
        std::string hmax0AsString = "H_MAX_0 INF";
        size_t iterationCount;
        for (iterationCount = 0 ; iterationCount < maxIteration ; iterationCount ++)
        {
            std::string cacheFileName="cache"+std::to_string(iterationCount)+".txt";
            SuggestMainStep->setName("Suggest " + NOMAD::itos(iterationCount));
            ObserveMainStep->setName("Observe " + NOMAD::itos(iterationCount));

            // Use separate Parameters for Suggest and Observe
            auto paramsForSuggestPtr = std::make_shared<NOMAD::AllParameters>();

            if(0 == iterationCount)
            {
                initParams(paramsForSuggestPtr,cacheFileName,initialFrameSizeAsString,hmax0AsString,true /*first iteration do not use cache*/);
            }
            else
            {
                initParams(paramsForSuggestPtr,cacheFileName,initialFrameSizeAsString,hmax0AsString,false /* not first iteration, use cache*/);
            }

            SuggestMainStep->setAllParameters(paramsForSuggestPtr);

            // MainStep runs suggest
            auto xs = SuggestMainStep->suggest();
            if (0 == xs.size())
            {
                std::cout << "No more points to suggest at iteration " << iterationCount << "." << std::endl;
                // Could not suggest any more points. Break.
                break;
            }

            // Not sure what to do with this shared_ptr. The setAllParameters requires a shared_ptr!
            paramsForSuggestPtr.reset();

            std::vector<NOMAD::ArrayOfDouble> fxs;
            eval_xs(xs,fxs);

            //THIS IS IMPORTANT
            // Need to reset cache, evaluator control, subproblem manager, seed
            NOMAD::MainStep::resetComponentsBetweenOptimization();


            // Parameters creation (important to create a fresh one because Suggest modifies its params (X0 from cache))
            auto paramsForObservePtr = std::make_shared<NOMAD::AllParameters>();
            
            if(0 == iterationCount)
            {
                initParams(paramsForObservePtr,cacheFileName,initialFrameSizeAsString,hmax0AsString,true /*first iteration does not use cache*/);
            }
            else
            {
                initParams(paramsForObservePtr,cacheFileName,initialFrameSizeAsString,hmax0AsString,false /* not first iteration, use cache*/);
            }
            
            ObserveMainStep->setAllParameters(paramsForObservePtr);

            std::cout << "Parameters before observe:" << std::endl;
            std::string paramName = "H_MAX_0";
            std::cout << paramName << " " << paramsForObservePtr->getAttributeValue<NOMAD::Double>(paramName) << std::endl;
            paramName = "INITIAL_FRAME_SIZE";
            std::cout << paramName << " ( " << paramsForObservePtr->getAttributeValue<NOMAD::ArrayOfDouble>(paramName) << " )" << std::endl;

            // MainStep runs observe
            std::string updatedCacheFileName="cache"+std::to_string(iterationCount+1)+".txt";

            auto updatedParams = ObserveMainStep->observe(xs,fxs,updatedCacheFileName);

            // Not sure what to do with this shared_ptr. The setAllParameters requires a shared_ptr!
            paramsForObservePtr.reset();

            std::cout << "Updated parameters: " << std::endl;
            for (auto p : updatedParams)
            {
                std::cout << p << std::endl;
                NOMAD::ParameterEntry pe(p);

                if (pe.getName() == "INITIAL_FRAME_SIZE")
                    initialFrameSizeAsString = p ;
                else if (pe.getName() == "H_MAX_0")
                    hmax0AsString = p   ;
                else
                {
                    std::string err("Exception: The update parameter is unknown");
                    err += "\n" + p;
                    throw std::logic_error(err);
                }

            }

            std::cout << "==============================================" <<std::endl;

            //THIS IS IMPORTANT
            // Need to reset cache, evaluator control, subproblem manager, seed
            NOMAD::MainStep::resetComponentsBetweenOptimization();

        }
        
        // The word "termination" is used to detect normal execution during tests.
        std::cout << "Completed the suggest and observe loop. Normal termination. Number of iterations completed: " << iterationCount << std::endl;
    }

    catch (std::exception &e)
    {
        std::cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
    }

    return EXIT_SUCCESS;
}
