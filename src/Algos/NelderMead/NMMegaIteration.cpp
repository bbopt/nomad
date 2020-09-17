
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/NelderMead/NMMegaIteration.hpp"
#include "../../Output/OutputQueue.hpp"


void NOMAD::NMMegaIteration::init()
{
    _name = getAlgoName() + NOMAD::MegaIteration::getName();

    // Get barrier from upper MadsMegaIteration, if available.
    auto madsMegaIter = getParentOfType<NOMAD::MadsMegaIteration*>(false);
    if (nullptr != madsMegaIter)
    {
        _barrier = madsMegaIter->getBarrier();
    }
}


void NOMAD::NMMegaIteration::startImp()
{
    // Create a Nelder Mead iteration for a frame center.
    // Use xFeas or xInf if XFeas is not available.
    // During NM we use a single iteration object with several start, run, end for the various iterations of the algorithm.

    if ( ! _stopReasons->checkTerminate() )
    {
        // MegaIteration's barrier member is already in sub dimension.
        auto bestXFeas = _barrier->getFirstXFeas();
        auto bestXInf  = _barrier->getFirstXInf();

        // Note: getParentOfType with argument "false" gets over the "Algorithm" parents.
        // Here, we are looking for a MadsMegaIteration which would be ancestor of
        // the NM (Algorithm) parent.
        auto madsMegaIter = getParentOfType<NOMAD::MadsMegaIteration*>(false);
        std::shared_ptr<NOMAD::MeshBase> mesh = nullptr;

        if ( madsMegaIter != nullptr )
        {
            mesh = madsMegaIter->getMesh();
        }

        if (nullptr != bestXFeas)
        {
            _nmIteration = std::make_shared<NOMAD::NMIteration>(this,
                                    std::make_shared<NOMAD::EvalPoint>(*bestXFeas),
                                    _k,
                                    mesh);
            _k++;
        }
        else if (nullptr != bestXInf)
        {
            _nmIteration = std::make_shared<NOMAD::NMIteration>(this,
                                    std::make_shared<NOMAD::EvalPoint>(*bestXInf),
                                    _k,
                                    mesh);
            _k++;
        }

        OUTPUT_DEBUG_START
        auto frameCenter = _nmIteration->getFrameCenter();
        AddOutputDebug("Frame center: " + frameCenter->display());
        auto previousFrameCenter = frameCenter->getPointFrom();
        AddOutputDebug("Previous frame center: " + (previousFrameCenter ? previousFrameCenter->display() : "NULL"));
        OUTPUT_DEBUG_END
    }
}


bool NOMAD::NMMegaIteration::runImp()
{
    bool successful = false;
    std::string s;

    if ( _stopReasons->checkTerminate() )
    {
        OUTPUT_DEBUG_START
        s = _name + ": stopReason = " + _stopReasons->getStopReasonAsString() ;
        AddOutputDebug(s);
        OUTPUT_DEBUG_END
        return false;
    }

    if ( _nmIteration == nullptr )
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "No iteration to run");
    }

    const size_t maxIter = _runParams->getAttributeValue<size_t>("MAX_ITERATION_PER_MEGAITERATION");
    size_t nbMegaIter = 0;
    while ( ! _stopReasons->checkTerminate() && nbMegaIter < maxIter )
    {
        _nmIteration->start();

        bool iterSuccessful = _nmIteration->run();          // Is this iteration successful
        successful = iterSuccessful || successful;  // Is the whole MegaIteration successful

        _nmIteration->end();

        if (iterSuccessful)
        {
            OUTPUT_DEBUG_START
            s = _name + ": new success " + NOMAD::enumStr(getSuccessType());
            AddOutputDebug(s);
            OUTPUT_DEBUG_END
        }

        if (_userInterrupt)
        {
            hotRestartOnUserInterrupt();
        }

        nbMegaIter++;
    }
    OUTPUT_DEBUG_START
    // Display MegaIteration's stop reason
    AddOutputDebug(_name + " stop reason set to: " + _stopReasons->getStopReasonAsString());
    OUTPUT_DEBUG_END

    // MegaIteration is a success if either a better xFeas or
    // a dominating or partial success for xInf was found.
    // See Algorithm 12.2 from DFBO.

    // return true if we have a partial or full success.
    return successful;
}


void NOMAD::NMMegaIteration::display( std::ostream& os ) const
{
// TODO display simplex
//    os << "MAIN_MESH " << std::endl;
//    os << *_mainMesh ;
    NOMAD::MegaIteration::display(os);
}


void NOMAD::NMMegaIteration::read(  std::istream& is )
{
    // TODO read simplex
    // Set up structures to gather member info
    size_t k=0;
    // Read line by line
    std::string name;
    while (is >> name && is.good() && !is.eof())
    {
//        if ("MAIN_MESH" == name)
//        {
//            if (nullptr != _mainMesh)
//            {
//                is >> *_mainMesh;
//            }
//            else
//            {
//                std::string err = "Error: Reading a mesh onto a NULL pointer";
//                std::cerr << err;
//            }
//        }
//        else
        if ("ITERATION_COUNT" == name)
        {
            is >> k;
        }
        else if ("BARRIER" == name)
        {
            if (nullptr != _barrier)
            {
                is >> *_barrier;
            }
            else
            {
                std::string err = "Error: Reading a Barrier onto a NULL pointer";
                std::cerr << err;
            }
        }
        else
        {
            for (size_t i = 0; i < name.size(); i++)
            {
                is.unget();
            }
            break;
        }
    }

    setK(k);
}


std::ostream& NOMAD::operator<<(std::ostream& os, const NOMAD::NMMegaIteration& megaIteration)
{
    megaIteration.display ( os );
    return os;
}


std::istream& NOMAD::operator>>(std::istream& is, NOMAD::NMMegaIteration& megaIteration)
{

    megaIteration.read( is );
    return is;

}
