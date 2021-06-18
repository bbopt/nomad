/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created by                                          */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4 is owned by                                 */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD 4 has been funded by Rio Tinto, Hydro-Québec, Huawei-Canada,             */
/*  NSERC (Natural Sciences and Engineering Research Council of Canada),           */
/*  InnovÉÉ (Innovation en Énergie Électrique) and IVADO (The Institute            */
/*  for Data Valorization)                                                         */
/*                                                                                 */
/*  NOMAD v3 was created and developed by Charles Audet, Sebastien Le Digabel,     */
/*  Christophe Tribes and Viviane Rochon Montplaisir and was funded by AFOSR       */
/*  and Exxon Mobil.                                                               */
/*                                                                                 */
/*  NOMAD v1 and v2 were created and developed by Mark Abramson, Charles Audet,    */
/*  Gilles Couture, and John E. Dennis Jr., and were funded by AFOSR and           */
/*  Exxon Mobil.                                                                   */
/*                                                                                 */
/*  Contact information:                                                           */
/*    Polytechnique Montreal - GERAD                                               */
/*    C.P. 6079, Succ. Centre-ville, Montreal (Quebec) H3C 3A7 Canada              */
/*    e-mail: nomad@gerad.ca                                                       */
/*                                                                                 */
/*  This program is free software: you can redistribute it and/or modify it        */
/*  under the terms of the GNU Lesser General Public License as published by       */
/*  the Free Software Foundation, either version 3 of the License, or (at your     */
/*  option) any later version.                                                     */
/*                                                                                 */
/*  This program is distributed in the hope that it will be useful, but WITHOUT    */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or          */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License    */
/*  for more details.                                                              */
/*                                                                                 */
/*  You should have received a copy of the GNU Lesser General Public License       */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.           */
/*                                                                                 */
/*  You can find information on the NOMAD software at www.gerad.ca/nomad           */
/*---------------------------------------------------------------------------------*/

#include "../Algos/CacheInterface.hpp"
#include "../Algos/SubproblemManager.hpp"
#include "../Cache/CacheBase.hpp"


void NOMAD::CacheInterface::init()
{
    _fixedVariable = NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(_step);

}


bool NOMAD::CacheInterface::smartInsert(const NOMAD::EvalPoint &evalPoint,
                                        const short maxNumberEval,
                                        const NOMAD::EvalType& evalType)
{
    // Always insert full dimension points.
    NOMAD::EvalPoint evalPointFull = evalPoint.makeFullSpacePointFromFixed(_fixedVariable);
    return NOMAD::CacheBase::getInstance()->smartInsert(evalPointFull, maxNumberEval, evalType);
}


size_t NOMAD::CacheInterface::find(const NOMAD::Point& x,
                                   NOMAD::EvalPoint &evalPoint,
                                   const NOMAD::EvalType& evalType)
{
    // Look for full dimension points.
    NOMAD::Point xFull = x.makeFullSpacePointFromFixed(_fixedVariable);
    size_t nbFound = NOMAD::CacheBase::getInstance()->find(xFull, evalPoint, evalType);
    if (nbFound > 0)
    {
        evalPoint = evalPoint.makeSubSpacePointFromFixed(_fixedVariable);
    }
    return nbFound;
}


size_t NOMAD::CacheInterface::findBestFeas(std::vector<NOMAD::EvalPoint> &evalPointList,
                                           const NOMAD::EvalType& evalType,
                                           const NOMAD::ComputeType& computeType,
                                           const NOMAD::Eval* refeval) const
{
    // Cache holds the full dimension points.
    // Return a list of sub dimension points.
    NOMAD::CacheBase::getInstance()->findBestFeas(evalPointList, _fixedVariable,
                                                  evalType, computeType, refeval);

    NOMAD::convertPointListToSub(evalPointList, _fixedVariable);

    return evalPointList.size();
}


size_t NOMAD::CacheInterface::findBestInf(std::vector<NOMAD::EvalPoint>& evalPointList,
                                          const NOMAD::Double& hMax,
                                          const NOMAD::EvalType& evalType,
                                          const NOMAD::ComputeType& computeType,
                                          const NOMAD::Eval* refeval) const
{
    // Cache holds the full dimension points.
    // Return a list of sub dimension points.
    NOMAD::CacheBase::getInstance()->findBestInf(evalPointList, hMax,
                                                 _fixedVariable, evalType,
                                                 computeType, refeval);

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
        auto critSubSpace1 = [&](const NOMAD::EvalPoint& evalPoint){return evalPoint.hasFixed(_fixedVariable);};

        // Make sure to convert an eval point coming from cache into subspace before calling crit1 function.
        auto critSubSpace2 = [&](const NOMAD::EvalPoint& evalPoint){ const NOMAD::EvalPoint & xSub = evalPoint.makeSubSpacePointFromFixed(_fixedVariable); return crit1(xSub);};

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
