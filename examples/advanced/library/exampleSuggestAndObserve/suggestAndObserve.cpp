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


void initParams(std::shared_ptr<NOMAD::AllParameters>& params, const std::string & cacheFileName )
{
    // Problem parameters
    params->setAttributeValue("DIMENSION", 2);             // number of variables
    params->setAttributeValue("LOWER_BOUND", NOMAD::ArrayOfDouble(2, -10.0));
    params->setAttributeValue("UPPER_BOUND", NOMAD::ArrayOfDouble(2, 10.0));

    NOMAD::BBOutputTypeList bbot;   // Definition of output types
    bbot.push_back(NOMAD::BBOutputType::OBJ);
    params->setAttributeValue("BB_OUTPUT_TYPE", bbot);

    
    params->setAttributeValue("CACHE_FILE", cacheFileName);
    
        // Algorithm parameters
    params->setAttributeValue("INITIAL_FRAME_SIZE",NOMAD::ArrayOfDouble(2,0.1) );
    params->setAttributeValue("MEGA_SEARCH_POLL",true); // Enable Mads MegaSearchPoll to perform suggest.

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
        auto SuggestMainStep = std::make_unique<NOMAD::MainStep>();

        // Parameters creation
        auto paramsForSuggest = std::make_shared<NOMAD::AllParameters>();
        initParams(paramsForSuggest,"cache0.txt");
        SuggestMainStep->setAllParameters(paramsForSuggest);

        // MainStep runs suggest
        auto xs = SuggestMainStep->suggest();

        std::vector<NOMAD::ArrayOfDouble> fxs;
        eval_xs(xs,fxs);
        
        //THIS IS IMPORTANT (see comments in function)
        NOMAD::MainStep::resetCache();

        auto ObserveMainStep = std::make_unique<NOMAD::MainStep>();

        // Parameters creation (important to create a fresh one because Suggest modifies its params (X0 from cache))
        auto paramsForObserve = std::make_shared<NOMAD::AllParameters>();
        initParams(paramsForObserve,"cache0.txt");  // Source cache (cache0.txt) different than destination cache (see below when calling observe)
        ObserveMainStep->setAllParameters(paramsForObserve);

        // MainStep runs observe
        std::cout << "Parameters before observe:" << std::endl;
        std::string paramName = "H_MAX_0";
        std::cout << paramName << " " << paramsForObserve->getAttributeValue<NOMAD::Double>(paramName) << std::endl;
        paramName = "INITIAL_FRAME_SIZE";
        std::cout << paramName << " ( " << paramsForObserve->getAttributeValue<NOMAD::ArrayOfDouble>(paramName) << " )" << std::endl;
        auto updateParams = ObserveMainStep->observe(xs,fxs,"cache1.txt");
        std::cout << "Updated parameters: " << std::endl;
        for (auto p : updateParams)
        {
            std::cout << p << std::endl;
        }
        
        // The word "termination" is used to detect normal execution during tests.
        std::cout << "Completion of a single pass of suggest and observe. Normal termination." << std::endl;
    }

    catch (std::exception &e)
    {
        std::cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
    }

    return EXIT_SUCCESS;
}
