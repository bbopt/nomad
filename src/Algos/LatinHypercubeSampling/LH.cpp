
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/LatinHypercubeSampling/LH.hpp"
#include "../../Math/LHS.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::LH::init()
{
    _name = "Latin Hypercube Sampling";
    verifyParentNotNull();

}

void NOMAD::LH::startImp()
{

    // Comment to appear at the end of stats lines
    setAlgoComment("(LH)");

    generateTrialPoints();

}

void NOMAD::LH::generateTrialPoints()
{
    OUTPUT_INFO_START
    AddOutputInfo("Generate points for " + _name, true, false);
    OUTPUT_INFO_END

    auto lhEvals = _runParams->getAttributeValue<size_t>("LH_EVAL");
    if (NOMAD::INF_SIZE_T == lhEvals)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "The number of evaluations for LH cannot be infinite.");
    }

    auto n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    auto lowerBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");

    if (!lowerBound.isComplete())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,_name + " requires a complete lower bound vector");
    }

    auto upperBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");
    if (!upperBound.isComplete())
    {
        throw NOMAD::Exception(__FILE__,__LINE__,_name + " requires a complete upper bound vector");
    }

    // Apply Latin Hypercube algorithm
    NOMAD::LHS lhs(n, lhEvals, lowerBound, upperBound);
    auto pointVector = lhs.Sample();

    for (auto point : pointVector)
    {
        // Make an EvalPoint from the Point.
        // We do not need the Eval part of EvalPoint right now,
        // but it will be used soon. Could be refactored, but
        // not high priority. Note that an EvalPointSet compares
        // the Point part of the EvalPoints only.

        // Projection without scale
        NOMAD::EvalPoint evalPoint(point);

        // Test if the point is inserted correctly
        bool inserted = insertTrialPoint(evalPoint);
        OUTPUT_INFO_START
        std::string s = "Generated point";
        s += (inserted) ? ": " : " not inserted: ";
        s += evalPoint.display();
        AddOutputInfo(s);
        OUTPUT_INFO_END
    }

    OUTPUT_INFO_START
    AddOutputInfo("Generated " + std::to_string(getTrialPointsCount()) + " points");
    AddOutputInfo("Generate points for " + _name, false, true);
    OUTPUT_INFO_END

}

bool NOMAD::LH::runImp()
{
    bool foundBetter = false;

    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);
    }
    auto LHStopReasons = NOMAD::AlgoStopReasons<NOMAD::LHStopType>::get( _stopReasons );
    if (  _stopReasons->testIf( NOMAD::EvalStopType::ALL_POINTS_EVALUATED ) )
    {
        LHStopReasons->set( NOMAD::LHStopType::ALL_POINTS_EVALUATED );
    }

    return foundBetter;
}

void NOMAD::LH::endImp()
{
    // Remove any remaining points from eval queue.
    EvcInterface::getEvaluatorControl()->clearQueue();

    // reset to the previous stats comment
    resetPreviousAlgoComment();

}
