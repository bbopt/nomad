
#include "../../Algos/Mads/MadsInitialization.hpp"
#include "../../Algos/SSDMads/SSDMads.hpp"
#include "../../Algos/SSDMads/SSDMadsMegaIteration.hpp"


void NOMAD::SSDMads::init()
{
    _name = "SSD-MADS";
    verifyParentNotNull();

    // Instanciate Mads initialization class
    _initialization = std::make_unique<NOMAD::MadsInitialization>( this );
}


void NOMAD::SSDMads::readInformationForHotRestart()
{
}


bool NOMAD::SSDMads::runImp()
{
    size_t k = 0;   // Iteration number
    NOMAD::SuccessType megaIterSuccess = NOMAD::SuccessType::NOT_EVALUATED;

    bool runOk = true;

    // Note: _initialization is run in Algorithm::startImp().

    if (!_termination->terminate(k))
    {
        std::shared_ptr<NOMAD::MeshBase> mesh = dynamic_cast<NOMAD::MadsInitialization*>(_initialization.get())->getMesh();


        auto barrier = _initialization->getBarrier();

        // Mads member _megaIteration is used for hot restart (read and write),
        // as well as to keep values used in Mads::end(), and may be used for _termination.
        // Update it here.
        _megaIteration = std::make_shared<NOMAD::SSDMadsMegaIteration>(this, k, barrier, mesh, megaIterSuccess);

        while (!_termination->terminate(k))
        {

            // Create a MegaIteration to manage the pollster worker and the regular workers.
            NOMAD::SSDMadsMegaIteration ssdMegaIteration(this, k, barrier, mesh, megaIterSuccess);
            ssdMegaIteration.start();
            ssdMegaIteration.run();
            ssdMegaIteration.end();

            // Remember these values to construct the next MegaIteration.
            k       = ssdMegaIteration.getNextK();
            barrier = ssdMegaIteration.getBarrier();
            mesh    = ssdMegaIteration.getMesh();
            megaIterSuccess = ssdMegaIteration.getSuccessType();

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
