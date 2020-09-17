
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/Mads/SearchMethodBase.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::SearchMethodBase::init()
{
    // A search method must have a parent
    verifyParentNotNull();

}


void NOMAD::SearchMethodBase::endImp()
{
    // Compute hMax and update Barrier.
    postProcessing(NOMAD::EvcInterface::getEvaluatorControl()->getEvalType());

    // Need to reimplement end() to set a stop reason for Mads based on the search method stop reason
}


void NOMAD::SearchMethodBase::generateTrialPoints()
{

    OUTPUT_INFO_START
    AddOutputInfo("Generate points for " + _name, true, false);
    OUTPUT_INFO_END

    generateTrialPointsImp();

    // Snap the points to bounds and mesh
    auto searchMethodPoints = getTrialPoints();
    auto lowerBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");
    auto upperBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");

    std::list<NOMAD::EvalPoint> snappedTrialPoints;
    for (auto point : searchMethodPoints)
    {
        if (snapPointToBoundsAndProjectOnMesh(point,lowerBound,upperBound))
        {
            snappedTrialPoints.push_back(NOMAD::EvalPoint(point));
            OUTPUT_INFO_START
            std::string s = "Snap point " + point.display();
            AddOutputInfo(s);
            OUTPUT_INFO_END
        }
    }

    // Re-insert snapped trial points
    clearTrialPoints();
    for (auto point : snappedTrialPoints)
    {
        insertTrialPoint(point);
    }

    OUTPUT_INFO_START
    AddOutputInfo("Generated " + std::to_string(getTrialPointsCount()) + " points");
    AddOutputInfo("Generate points for " + _name, false, true);
    OUTPUT_INFO_END


    // The trial points must know what frame center originated them.
    updatePointsWithFrameCenter();

}
