
#include "../Algos/CacheInterface.hpp"
#include "../Algos/SubproblemManager.hpp"
#include "../Cache/CacheBase.hpp"


void NOMAD::CacheInterface::init()
{
    _fixedVariable = NOMAD::SubproblemManager::getSubFixedVariable(_step);

}


bool NOMAD::CacheInterface::smartInsert(const NOMAD::EvalPoint &evalPoint,
                                        const short maxNumberEval,
                                        const EvalType& evalType)
{
    // Always insert full dimension points.
    NOMAD::EvalPoint evalPointFull = evalPoint.makeFullSpacePointFromFixed(_fixedVariable);
    return NOMAD::CacheBase::getInstance()->smartInsert(evalPointFull, maxNumberEval, evalType);
}


size_t NOMAD::CacheInterface::find(const NOMAD::Point x,
                                   NOMAD::EvalPoint &evalPoint)
{
    // Look for full dimension points.
    NOMAD::Point xFull = x.makeFullSpacePointFromFixed(_fixedVariable);
    size_t nbFound = NOMAD::CacheBase::getInstance()->find(xFull, evalPoint);
    evalPoint = evalPoint.makeSubSpacePointFromFixed(_fixedVariable);
    return nbFound;
}


size_t NOMAD::CacheInterface::findBestFeas(std::vector<NOMAD::EvalPoint> &evalPointList,
                                           const EvalType& evalType,
                                           const Eval* refeval) const
{
    // Cache holds the full dimension points.
    // Return a list of sub dimension points.
    NOMAD::CacheBase::getInstance()->findBestFeas(evalPointList, _fixedVariable,
                                                  evalType, refeval);

    NOMAD::convertPointListToSub(evalPointList, _fixedVariable);

    return evalPointList.size();
}


size_t NOMAD::CacheInterface::findBestInf(std::vector<NOMAD::EvalPoint>& evalPointList,
                                          const NOMAD::Double& hMax,
                                          const EvalType& evalType,
                                          const Eval* refeval) const
{
    // Cache holds the full dimension points.
    // Return a list of sub dimension points.
    NOMAD::CacheBase::getInstance()->findBestInf(evalPointList, hMax,
                                                 _fixedVariable, evalType,
                                                 refeval);

    NOMAD::convertPointListToSub(evalPointList, _fixedVariable);

    return evalPointList.size();
}


size_t NOMAD::CacheInterface::find(std::function<bool(const NOMAD::EvalPoint&)> crit1,
                                   std::vector<NOMAD::EvalPoint> &evalPointList,
                                   bool findInSubspace ) const
{
    if ( findInSubspace )
    {
        // Lambda function to test if an eval point is in the current subspace (its  variables are consistent with the fixed values in this interface)
        auto critSubSpace1 = [&](const EvalPoint& evalPoint){return evalPoint.hasFixed(_fixedVariable);};

        // Make sure to convert an eval point coming from cache into subspace before calling crit1 function.
        auto critSubSpace2 = [&](const EvalPoint& evalPoint){ const NOMAD::EvalPoint & xSub = evalPoint.makeSubSpacePointFromFixed(_fixedVariable); return crit1(xSub);};

        NOMAD::CacheBase::getInstance()->find(critSubSpace1, critSubSpace2, evalPointList);

    }
    else
    {
        NOMAD::CacheBase::getInstance()->find(crit1, evalPointList);
    }

    NOMAD::convertPointListToSub(evalPointList, _fixedVariable);

    return evalPointList.size();
}

size_t NOMAD::CacheInterface::getAllPoints(std::vector<NOMAD::EvalPoint> &evalPointList) const
{
    NOMAD::CacheBase::getInstance()->find(
        [this] (const NOMAD::EvalPoint& evalPoint) {
            return evalPoint.hasFixed(_fixedVariable);
        },
        evalPointList);


    NOMAD::convertPointListToSub(evalPointList, _fixedVariable);

    return evalPointList.size();
}
