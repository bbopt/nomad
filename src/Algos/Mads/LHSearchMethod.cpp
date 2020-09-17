
#include "../../Algos/Mads/LHSearchMethod.hpp"
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Math/LHS.hpp"
#include "../../Type/LHSearchType.hpp"

void NOMAD::LHSearchMethod::init()
{
    _name = "Latin Hypercube Search Method";
    //setComment("(LHSearch)");

    auto lhSearch = _runParams->getAttributeValue<NOMAD::LHSearchType>("LH_SEARCH");
    setEnabled(lhSearch.isEnabled());
}


void NOMAD::LHSearchMethod::generateTrialPointsImp()
{

    if (nullptr == _iterAncestor)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"LHSearchMethod: must have an iteration ancestor");
    }
    auto mesh = _iterAncestor->getMesh();
    if (nullptr == mesh)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"LHSearchMethod: must have a mesh");
    }
    auto frameCenter = _iterAncestor->getFrameCenter();
    if (nullptr == frameCenter)
    {
        throw NOMAD::Exception(__FILE__,__LINE__,"LHSearchMethod: must have a frameCenter");
    }

    auto lhSearch = _runParams->getAttributeValue<NOMAD::LHSearchType>("LH_SEARCH");
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    size_t p = (0 == _iterAncestor->getK()) ? lhSearch.getNbInitial() : lhSearch.getNbIteration();
    auto lowerBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");
    auto upperBound = _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");


    // Update undefined values of lower and upper bounds to use values based
    // on DeltaFrameSize.
    // Based on the code in NOMAD 3, but slightly different.
    // If we used INF values instead of these, we get huge values for the
    // generated points. It is not elegant.
    NOMAD::ArrayOfDouble deltaFrameSize = mesh->getDeltaFrameSize();
    NOMAD::Double scaleFactor = sqrt(-log(NOMAD::DEFAULT_EPSILON));

    for (size_t i = 0; i < n; i++)
    {
        if (!lowerBound[i].isDefined())
        {
            lowerBound[i] = (*frameCenter)[i] - 10.0 * deltaFrameSize[i] * scaleFactor;
        }
        if (!upperBound[i].isDefined())
        {
            upperBound[i] = (*frameCenter)[i] + 10.0 * deltaFrameSize[i] * scaleFactor;
        }
    }

    // Apply Latin Hypercube algorithm
    NOMAD::LHS lhs(n, p, lowerBound, upperBound);
    auto pointVector = lhs.Sample();

    // Insert the point. Projection on mesh and snap to bounds is done in SearchMethod
    for (auto point : pointVector)
    {
        // Insert point (if possible)
        insertTrialPoint(NOMAD::EvalPoint(point));

    }
}
