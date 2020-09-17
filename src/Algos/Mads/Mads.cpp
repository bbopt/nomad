
#include "../../Algos/Mads/GMesh.hpp"
#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/Mads/MadsInitialization.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Output/OutputQueue.hpp"
#include "../../Util/fileutils.hpp"

void NOMAD::Mads::init()
{
    _name = "MADS";

    // Instanciate Mads initialization class
    _initialization = std::make_unique<NOMAD::MadsInitialization>( this );

}


bool NOMAD::Mads::runImp()
{
    size_t k = 0;   // Iteration number
    NOMAD::SuccessType megaIterSuccess = NOMAD::SuccessType::NOT_EVALUATED;

    bool runOk = true;

    if (!_termination->terminate(k))
    {
        std::shared_ptr<NOMAD::MeshBase> mesh;

        // TODO Review case hot restart.
        /*
        if (nullptr != _megaIteration)
        {
            // Case hot restart
            k       = _megaIteration->getK();
            barrier = _megaIteration->getBarrier();

            // Downcast from MegaIteration to MadsMegaIteration
            mesh    = (std::dynamic_pointer_cast<NOMAD::MadsMegaIteration> (_megaIteration ))->getMesh();
            megaIterSuccess = _megaIteration->getSuccessType();
        }
        */

        mesh = dynamic_cast<NOMAD::MadsInitialization*>(_initialization.get())->getMesh();
        auto barrier = _initialization->getBarrier();

        // Mads member _megaIteration is used for hot restart (read and write),
        // as well as to keep values used in Mads::end(), and may be used for _termination.
        // Update it here.
        _megaIteration = std::make_shared<NOMAD::MadsMegaIteration>(this, k, barrier, mesh, megaIterSuccess);


        while (!_termination->terminate(k))
        {
            // Create an MegaIteration: manage multiple iterations on different
            // meshes and with different frame centers at the same time.
            NOMAD::MadsMegaIteration megaIteration(this, k, barrier, mesh, megaIterSuccess);
            megaIteration.start();
            megaIteration.run();
            megaIteration.end();

            // Remember these values to construct the next MegaIteration.
            k       = megaIteration.getNextK();
            barrier = megaIteration.getBarrier();
            mesh    = megaIteration.getMesh();
            megaIterSuccess = megaIteration.getSuccessType();

            if (_userInterrupt)
            {
                hotRestartOnUserInterrupt();
            }
        }
    }

    else
    {
        runOk = false;
    }

    _termination->start();
    _termination->run();
    _termination->end();

    return runOk;
}


void NOMAD::Mads::readInformationForHotRestart()
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
            std::string s = "Read hot restart file " + hotRestartFile;
            NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_NORMAL);

            // Create a GMesh and an MadsMegaIteration with default values, to be filled
            // by istream is.
            // NOTE: Working in full dimension
            auto barrier = std::make_shared<NOMAD::Barrier>(NOMAD::INF, NOMAD::Point(), NOMAD::EvalType::BB);
            std::shared_ptr<NOMAD::MeshBase> mesh = std::make_shared<NOMAD::GMesh>(_pbParams);

            _megaIteration = std::make_shared<NOMAD::MadsMegaIteration>(this, 0, barrier, mesh, NOMAD::SuccessType::NOT_EVALUATED);

            // Here we use Algorithm::operator>>
            NOMAD::read<NOMAD::Mads>(*this, hotRestartFile);
        }
    }
}
