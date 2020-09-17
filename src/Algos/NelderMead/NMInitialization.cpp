
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/NelderMead/NMInitialization.hpp"
#include "../../Algos/SubproblemManager.hpp"


void NOMAD::NMInitialization::init()
{
    _name = getAlgoName() + "Initialization";

    _nmStopReason = NOMAD::AlgoStopReasons<NOMAD::NMStopType>::get( _stopReasons );
}


bool NOMAD::NMInitialization::runImp()
{
    bool doContinue = ! _stopReasons->checkTerminate();

    if (doContinue)
    {
        // For a standalone NM, evaluate the trial points generated during start (simplex is created later)
        // Otherwise, there are no trial points available
        evalTrialPoints(this);
        doContinue = ! _stopReasons->checkTerminate();
        if ( ! doContinue )
            _nmStopReason->set(NOMAD::NMStopType::INITIAL_FAILED);

    }
    return doContinue;
}

void NOMAD::NMInitialization::startImp()
{

    if ( ! _stopReasons->checkTerminate() )
    {
        // If needed, generate trial points and put them in cache to form simplex
        // For a standalone optimization (NM_OPTIMIZATION true), initial trial points must be generated to form a valid simplex around x0. Otherwise, the cache will be used to construct the simplex.
        auto nm_opt = _runParams->getAttributeValue<bool>("NM_OPTIMIZATION");
        if ( nm_opt && ! checkCacheCanFormSimplex() )
        {
            generateTrialPoints();
        }
    }

}


void NOMAD::NMInitialization::endImp()
{
    // Construct _barrier member with evaluated _trialPoints for future use
    // _trialPoints are already updated with Evals.
    if (_trialPoints.size() > 0)
    {
        std::vector<NOMAD::EvalPoint> evalPointList;
        std::copy(_trialPoints.begin(), _trialPoints.end(),
                          std::back_inserter(evalPointList));
        auto hMax = _runParams->getAttributeValue<NOMAD::Double>("H_MAX_0");
        _barrier = std::make_shared<NOMAD::Barrier>(hMax, NOMAD::SubproblemManager::getSubFixedVariable(this), NOMAD::EvcInterface::getEvaluatorControl()->getEvalType(), evalPointList);
    }
}


bool NOMAD::NMInitialization::checkCacheCanFormSimplex()
{
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    if ( NOMAD::CacheBase::getInstance()->size() < n+1 )
        return false;
    // TODO
    return false;

}

// Generate trial points to form a simplex
void NOMAD::NMInitialization::generateTrialPoints()
{
    NOMAD::Point x0 = _pbParams->getAttributeValue<NOMAD::Point>("X0");
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");

    if (!x0.isComplete() || x0.size() != n)
    {
        std::string err = "Initialization: evalY0: Invalid X0 " + x0.display();
        size_t cacheSize = NOMAD::CacheBase::getInstance()->size();
        if (cacheSize > 0)
        {
            err += ". Hint: Try not setting X0 so that the cache is used (";
            err += std::to_string(cacheSize) + " points)";
        }
        else
        {
            err += ". Cache is empty.";
        }
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    NOMAD::EvalPoint evalPoint_x0(x0);
    insertTrialPoint(evalPoint_x0);
    OUTPUT_INFO_START
    AddOutputInfo("Using X0: " + evalPoint_x0.display());
    OUTPUT_INFO_END

    // Method to generate simplex points using X0 adapted from fminsearch (matlab)
    const NOMAD::Double usualDelta = 0.05;    //  x0 + 5 percent
    const NOMAD::Double zeroDelta = 0.00025;  //
    for ( size_t j = 0 ; j < n ; j++ )
    {
        NOMAD::EvalPoint trialPoint(x0);
        if ( trialPoint[j] != 0 )
            trialPoint[j] *= (1 + usualDelta );
        else
            trialPoint[j] = zeroDelta;

        insertTrialPoint(trialPoint);
    }

    OUTPUT_INFO_START
    NOMAD::OutputQueue::Flush();
    OUTPUT_INFO_END
}
