
#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/Mads/NMSearchMethod.hpp"
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/NelderMead/NM.hpp"
#include "../../Algos/NelderMead/NMAllReflective.hpp"

void NOMAD::NMSearchMethod::init()
{
    if ( _runParams->getAttributeValue<bool>("GENERATE_ALL_POINTS_BEFORE_EVAL") )
    {
        _name = "Search (Nelder Mead single pass)";
    }
    else
    {
        _name = "Search (Nelder Mead optimization)";
    }
    //setComment("(NMSearch)");

    auto nmSearch = _runParams->getAttributeValue<bool>("NM_SEARCH");
    setEnabled(nmSearch);

    if (nmSearch)
    {
        // Set the lap counter
        auto nmFactor = _runParams->getAttributeValue<size_t>("NM_SEARCH_MAX_TRIAL_PTS_NFACTOR");
        auto dim = _pbParams->getAttributeValue<size_t>("DIMENSION");
        if (nmFactor < NOMAD::INF_SIZE_T)
        {
            NOMAD::EvcInterface::getEvaluatorControl()->setLapMaxBbEval( dim*nmFactor );
        }
    }
}


bool NOMAD::NMSearchMethod::runImp()
{
    // NM is an algorithm with its own stop reasons.
    auto nmStopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::NMStopType>>();

    // Create the NM algorithm with its own stop reason
    auto nm = std::make_shared<NOMAD::NM>(this,
                                          nmStopReasons ,
                                          _runParams,
                                          _pbParams);
    nm->setEndDisplay(false);

    nm->start();
    bool foundBetter = nm->run();
    nm->end();

    return foundBetter;
}


void NOMAD::NMSearchMethod::generateTrialPointsImp()
{
    // The trial points of one iteration of NM reflective steps are generated (not evaluated).
    // The trial points are Reflect, Expansion, Inside and Outside Contraction NM points

    auto madsIteration = getParentOfType<MadsIteration*>();

    NOMAD::NMAllReflective allReflective(this, madsIteration->getFrameCenter(), madsIteration->getMesh());
    allReflective.start();
    allReflective.end();

    // Pass the generated trial pts to this
    auto trialPtsNM = allReflective.getTrialPoints();
    for (auto point : trialPtsNM)
    {
        insertTrialPoint(point);
    }

}
