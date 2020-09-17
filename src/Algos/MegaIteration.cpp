
#include "../Algos/MegaIteration.hpp"


// Constructor
NOMAD::MegaIteration::MegaIteration(const Step* parentStep,
                              size_t k,
                              std::shared_ptr<Barrier> barrier,
                              SuccessType success)
  : Step(parentStep),
    _iterList(),
    _barrier(barrier),
    _k(k),
    _nbIterRun(0),
    _megaIterationSuccess(success)
{
    if (nullptr == _barrier)
    {
        throw Exception(__FILE__, __LINE__, "MegaIteration constructor: barrier must not be NULL.");
    }

    init();
}


void NOMAD::MegaIteration::init()
{
    _name = getAlgoName() + "MegaIteration " + std::to_string(_k);
    verifyParentNotNull();
}


size_t NOMAD::MegaIteration::getNextK() const
{
    size_t nextK;
    if (_nbIterRun > 0)
    {
        // Next MegaIteration number will start where the Iterations ended.
        nextK = _k + _nbIterRun;
    }
    else
    {
        // Must increment MegaIteration number, cannot have
        // twice the same number.
        nextK = _k + 1;
    }

    return nextK;
}


void NOMAD::MegaIteration::endImp()
{
    if (_runParams->getAttributeValue<bool>("USER_CALLS_ENABLED"))
    {
        bool stop = false;
        runCallback(NOMAD::CallbackType::MEGA_ITERATION_END, *this, stop);
        if (!_stopReasons->checkTerminate() && stop)
        {
            _stopReasons->set(NOMAD::BaseStopType::USER_STOPPED);
        }
    }

    _iterList.clear();
}


void NOMAD::MegaIteration::computeMaxXFeasXInf(size_t &maxXFeas, size_t &maxXInf)
{
    const size_t maxIter = _runParams->getAttributeValue<size_t>("MAX_ITERATION_PER_MEGAITERATION");
    const size_t maxXFeas0 = maxXFeas;
    const size_t maxXInf0 = maxXInf;

    // If maxXFeas + maxXInf does not exceed maxIter, do nothing.
    if (maxXFeas + maxXInf > maxIter)
    {
        if (maxXFeas <= maxIter / 2)
        {
            // Use all xFeas, and the remaining in xInf.
            maxXInf = maxIter - maxXFeas;
        }
        else if (maxXInf < maxIter / 2)
        {
            // Use all xInf, and the remaining in xFeas.
            maxXFeas = maxIter - maxXInf;
        }
        else
        {
            // Both the number of xFeas and xInf is over.
            // Use half of maxIter for xFeas and half for xInf.
            maxXInf = maxIter / 2;
            maxXFeas = maxIter - maxXInf;
        }
        if (maxXFeas + maxXInf > maxIter)
        {
            // This case should not happen and should be debugged.
            std::cerr << "Warning: Bad computation in computeMaxXFeasXInf. maxIter = " << maxIter << " maxXFeas = " << maxXFeas << " (was " << maxXFeas0 << ") maxXInf = " << maxXInf << " (was " << maxXInf0 << ")" << std::endl;
        }
    }
}


void NOMAD::MegaIteration::read(std::istream& is)
{
    // Set up structures to gather member info
    size_t k = 0;

    // Read line by line
    std::string name;
    while (is >> name && is.good() && !is.eof())
    {
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
                std::cerr << err << std::endl;
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


void NOMAD::MegaIteration::display(  std::ostream& os ) const
{
    os << "ITERATION_COUNT " << _k << std::endl;
    os << "BARRIER " << std::endl;
    os << *_barrier;
}


std::ostream& NOMAD::operator<<(std::ostream& os, const NOMAD::MegaIteration& megaIteration)
{

    megaIteration.display( os );
    return os;
}


std::istream& NOMAD::operator>>(std::istream& is, NOMAD::MegaIteration& megaIteration)
{
    megaIteration.read(is);
    return is;
}
