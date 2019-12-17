/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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
#include "../../Algos/MainStep.hpp"
#include "../../Algos/Mads/SearchMethod.hpp"

#include "../../Algos/EvcInterface.hpp"

// NM specific
#include "../../Algos/NelderMead/NMInitialization.hpp"
#include "../../Algos/NelderMead/NMMegaIteration.hpp"
#include "../../Algos/NelderMead/NM.hpp"

void NOMAD::NM::init()
{
    _name = "NM";
    if ( _runParams->getAttributeValue<bool>("GENERATE_ALL_POINTS_BEFORE_EVAL") )
    {
        _name += " One Iteration";
    }

    // Instanciate NM initialization class
    _initialization = std::make_unique<NOMAD::NMInitialization>( this );
}

void NOMAD::NM::startImp()
{

    // Comment to appear at the end of stats lines
    NOMAD::MainStep::setAlgoComment("(NM)");

    // All stop reasons are reset.
    _stopReasons->setStarted();

    // Reset the lapBbEval counter for this sub-optimization
    NOMAD::EvcInterface::getEvaluatorControl()->resetLapBbEval();

    _initialization->start();
    _initialization->run();
    _initialization->end();

}

bool NOMAD::NM::runImp()
{
    bool successful = false;

    NOMAD::SuccessType bestSuccess = NOMAD::SuccessType::NOT_EVALUATED;

    if ( ! _stopReasons->checkTerminate() )
    {
        size_t k = 0;   // Iteration number

        // Barrier constructor automatically finds the best points in the cache.
        auto hMax = _runParams->getAttributeValue<NOMAD::Double>("H_MAX_0");
        auto barrier = std::make_shared<NOMAD::Barrier>(hMax, getSubFixedVariable(), getEvalType());
        NOMAD::SuccessType megaIterSuccess = NOMAD::SuccessType::NOT_EVALUATED;

        if (nullptr != _megaIteration)
        {
            // Case hot restart
            k       = _megaIteration->getK();
            barrier = _megaIteration->getBarrier();
            megaIterSuccess = _megaIteration->getSuccessType();
        }

        while (!_termination->terminate(k))
        {
            // Create an MegaIteration: manage multiple iterations.
            NOMAD::NMMegaIteration megaIteration(this, k, barrier, megaIterSuccess);
            megaIteration.start();
            bool currentMegaIterSuccess = megaIteration.run();
            megaIteration.end();

            successful = successful || currentMegaIterSuccess;

            // Remember these values to construct the next MegaIteration.
            k       = megaIteration.getK();
            barrier = megaIteration.getBarrier();
            megaIterSuccess = megaIteration.getSuccessType();

            if ( megaIterSuccess > bestSuccess )
            {
                bestSuccess = megaIterSuccess;
            }

            if (_userInterrupt)
            {
                hotRestartOnUserInterrupt();
            }
        }

        // TODO -->for hot restart make sure to save the simplex (maybe as X0s)
        // _megaIteration is used for hot restart (read
        // and write), as well as to keep values used in Mads::end(). Update it here.
        _megaIteration = std::make_shared<NOMAD::NMMegaIteration>(this, k, barrier, megaIterSuccess);

        _termination->start();
        _termination->run();
        _termination->end();

        // CT TODO Maybe move this (and the equivalent in Mads) into Algorithm::end(). Need to add _megaIterSuccess and _bestSuccess as private attributes of Algorithm. Do not forget initialization.
        // Update the SearchMethod success type with best success found.
        if ( successful )
        {
            // The parent can be a SearchMethod (NM-Mads Search) or not (NM standalone optimization)
            auto searchMethodConst = dynamic_cast<const NOMAD::SearchMethod*>(_parentStep);

            if ( searchMethodConst != nullptr )
            {
                auto searchMethod = const_cast<NOMAD::SearchMethod*>(searchMethodConst);
                searchMethod->setSuccessType(bestSuccess);
            }

        }

    }

    return successful;
}


void NOMAD::NM::readInformationForHotRestart()
{
    // Restart from where we were before.
    // For this, we need to read some files.
    // Note: Cache file is treated independently from hot restart file.

    if (_runParams->getAttributeValue<bool>("HOT_RESTART_READ_FILES"))
    {
        // Verify the files exist and are readable.
        std::string hotRestartFile = _runParams->getAttributeValue<std::string>("HOT_RESTART_FILE");
        if (NOMAD::checkReadFile(hotRestartFile))
        {
            std::cout << "Read hot restart file " << hotRestartFile << std::endl;

            // Create a GMesh and a MegaIteration with default values, to be filled
            // by istream is.
            // TODO Fix potential bug with Hot Restart
            // Note: Assuming the barrier read is in the same subspace as the current subspace.
            // This could be fixed if we write and read the barrier in full subspace.
            auto barrier = std::make_shared<NOMAD::Barrier>();
            int k = 0;
            NOMAD::SuccessType success = NOMAD::SuccessType::NOT_EVALUATED;


            _megaIteration = std::make_shared<NOMAD::NMMegaIteration>(this, k, barrier, success);

            // Here we use NM::operator>>
            NOMAD::read<NM>(*this, hotRestartFile);
        }
    }
}
