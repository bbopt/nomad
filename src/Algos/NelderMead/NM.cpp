
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/Mads/SearchMethodBase.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Util/fileutils.hpp"

// NM specific
#include "../../Algos/NelderMead/NM.hpp"
#include "../../Algos/NelderMead/NMInitialization.hpp"
#include "../../Algos/NelderMead/NMMegaIteration.hpp"

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
    setAlgoComment("(NM)");

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

        std::shared_ptr<NOMAD::Barrier> barrier = nullptr;

        if (_runParams->getAttributeValue<bool>("NM_OPTIMIZATION"))
        {
            // Barrier was computed by Initialization.
            barrier = _initialization->getBarrier();
        }
        else
        {
            // Get barrier from upper MadsMegaIteration, if available.
            auto madsMegaIter = getParentOfType<NOMAD::MadsMegaIteration*>(false);
            if (nullptr != madsMegaIter)
            {
                barrier = madsMegaIter->getBarrier();
            }
        }

        NOMAD::SuccessType megaIterSuccess = NOMAD::SuccessType::NOT_EVALUATED;

        // TODO fix this case
        /*
        if (nullptr != _megaIteration)
        {
            // Case hot restart
            k       = _megaIteration->getK();
            barrier = _megaIteration->getBarrier();
            megaIterSuccess = _megaIteration->getSuccessType();
        }
        */

        while (!_termination->terminate(k))
        {
            // Create a MegaIteration: manage multiple iterations.
            NOMAD::NMMegaIteration megaIteration(this, k, barrier, megaIterSuccess);
            megaIteration.start();
            bool currentMegaIterSuccess = megaIteration.run();
            megaIteration.end();

            successful = successful || currentMegaIterSuccess;

            // Remember these values to construct the next MegaIteration.
            k       = megaIteration.getNextK();
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
            // The parent can be a SearchMethod (NM-Mads Search) or not (that is NM is a standalone optimization)
            auto searchMethodConst = dynamic_cast<const NOMAD::SearchMethodBase*>(_parentStep);

            if (searchMethodConst != nullptr)
            {
                auto searchMethod = const_cast<NOMAD::SearchMethodBase*>(searchMethodConst);
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
