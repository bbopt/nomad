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
/*--------------------------------------------------------*/
/*  how to use the NOMAD library to suggest trial points  */
/*--------------------------------------------------------*/
#include "Nomad/nomad.hpp"
#include "Cache/CacheSet.hpp"
#include "utils.hpp"


/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main(int argc, char ** argv)
{
    try
    {
        if (argc < 5)
        {
            std::cerr << "Usage: " << argv[0] << " <xs file> <fxs file> <input cache file> <output cache file> [optional param file]" << std::endl;
            return 1;
        }

        auto main = std::make_unique<NOMAD::MainStep>();

        std::string cacheFileName = argv[3];

        auto params = std::make_shared<NOMAD::AllParameters>();
        std::vector<std::string> additionalParams;
        if (argc >= 4)
        {
            additionalParams = readAdditionalParams(argv[5]);
        }
        initParams(params, cacheFileName, additionalParams);
        main->setAllParameters(params);
        NOMAD::OutputQueue::getInstance()->initParameters(params->getDispParams());

        // Initialize cache
        NOMAD::CacheSet::setInstance(params->getCacheParams(),
                            params->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE"));

        // Read xs and fxs files
        NOMAD::ArrayOfPoint xs;
        readXFile(argv[1], xs);
        std::vector<NOMAD::ArrayOfDouble> fxs = readBBOutputFile(argv[2], params->getAttributeValue<NOMAD::BBOutputTypeList>("BB_OUTPUT_TYPE").size());

        // MainStep runs observe
        std::string updatedCacheFileName = argv[4];

        auto updatedParams = main->observe(xs, fxs, updatedCacheFileName);

        for (auto p : updatedParams)
        {
            std::cout << p << std::endl;
            NOMAD::ParameterEntry pe(p);
        }
    }

    catch (std::exception &e)
    {
        std::cerr << "\n" << argv[0] << " has been interrupted (" << e.what() << ")\n\n";
    }

    return EXIT_SUCCESS;
}
