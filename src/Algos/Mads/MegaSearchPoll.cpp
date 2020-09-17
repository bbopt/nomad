
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/Mads/MegaSearchPoll.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Mads/Search.hpp"
#include "../../Algos/Mads/Poll.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::MegaSearchPoll::init()
{
    _name = "MegaSearchPoll";
    verifyParentNotNull();

    auto megaIter = dynamic_cast<const NOMAD::MadsMegaIteration*>( _megaIterAncestor );
    if (nullptr == megaIter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"An instance of class MegaSearch must have a MadsMegaIteration among its ancestors");
    }

}


const std::shared_ptr<NOMAD::MadsIteration> NOMAD::MegaSearchPoll::getIterForPoint(const NOMAD::EvalPoint& point) const
{
    return _iterForPoint[point];
}


void NOMAD::MegaSearchPoll::startImp()
{

    // Generate trial points using poll and search and merge them
    generateTrialPoints();

}


bool NOMAD::MegaSearchPoll::runImp()
{
    bool foundBetter = false;

    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);
    }

    // Update MegaIteration success type with best success found.
    _megaIterAncestor->setSuccessType(_success);

    return foundBetter;
}


void NOMAD::MegaSearchPoll::endImp()
{
    postProcessing(NOMAD::EvcInterface::getEvaluatorControl()->getEvalType());
}


// Generate new points to evaluate from Poll and Search
void NOMAD::MegaSearchPoll::generateTrialPoints()
{
    verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, true);
    OUTPUT_INFO_START
    AddOutputInfo("Generate points for " + _name, true, false);
    OUTPUT_INFO_END

    NOMAD::EvalPointSet trialPoints;

    // Generate points for all frame centers, all meshes, all Search and Poll strategies.
    for (size_t i = 0; i < _megaIterAncestor->getNbIterations() ; i++)
    {
        // downcast from Iteration to MadsIteration
        const std::shared_ptr<NOMAD::MadsIteration> & iter = std::dynamic_pointer_cast<NOMAD::MadsIteration> ( _megaIterAncestor->getIter(i));

        if ( iter == nullptr )
            throw NOMAD::Exception(__FILE__, __LINE__, "Cannot convert to MadsIteration shared pointer");

        // Generate trial points for Search (all enabled search methods) and Poll.
        // Note: Search and Poll generateTrialPoints() methods both
        // take care of verifying that the generated are on mesh, and also
        // update the "PointFrom" with the Iteration frame center.
        NOMAD::Search search(iter.get() );
        search.generateTrialPoints();
        auto trialPointsSearch = search.getTrialPoints() ;

        NOMAD::Poll poll(iter.get() );
        poll.generateTrialPoints();
        auto trialPointsPoll = poll.getTrialPoints();

        // Merge two sets and remove duplicates
        // Naive implementation. Easier to understand - I could not make std::merge,
        // std::unique or std::set_union work fine.
        // Caveat: Multiple EvalPoints copy.
        for (auto point : trialPointsSearch)
        {
            insertTrialPoint( point );
            // Remember which iteration generated these points
            auto pointIterPair = std::pair<NOMAD::EvalPoint, std::shared_ptr<NOMAD::MadsIteration>>(point, std::make_shared<NOMAD::MadsIteration>(*iter));
            _iterForPoint.insert(pointIterPair);
        }
        for (auto point : trialPointsPoll)
        {
            insertTrialPoint( point );

            // Remember which iteration generated these points
            auto pointIterPair = std::pair<NOMAD::EvalPoint, std::shared_ptr<NOMAD::MadsIteration>>(point, std::make_shared<NOMAD::MadsIteration>(*iter));
            _iterForPoint.insert(pointIterPair);
        }

    }

    OUTPUT_INFO_START
    AddOutputInfo("Generated " + NOMAD::itos(getTrialPointsCount()) + " points");
    AddOutputInfo("Generate points for " + _name, false, true);
    OUTPUT_INFO_END

}
