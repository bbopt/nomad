/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural Science    */
/*  and Engineering Research Council of Canada), INOVEE (Innovation en Energie     */
/*  Electrique and IVADO (The Institute for Data Valorization)                     */
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
/*    phone : 1-514-340-6053 #6928                                                 */
/*    fax   : 1-514-340-5665                                                       */
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

#include <algorithm>    // For std::merge and std::unique
#include <sstream>

#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/Mads/Search.hpp"
#include "../../Algos/Mads/Poll.hpp"

void NOMAD::MadsIteration::init()
{
    _name = getAlgoName() + NOMAD::Iteration::getName();
}


bool NOMAD::MadsIteration::runImp()
{
    verifyGenerateAllPointsBeforeEval(__PRETTY_FUNCTION__, false);

    bool iterationSuccess = false;
    NOMAD::SuccessType bestSuccessYet = NOMAD::SuccessType::NOT_EVALUATED;

    // Parameter Update is handled at the upper level - MegaIteration.

    // 1. Search
    if ( ! _stopReasons->checkTerminate() )
    {
        NOMAD::Search search(this );
        search.start();
        iterationSuccess = search.run();

        NOMAD::SuccessType success = search.getSuccessType();
        if (success > bestSuccessYet)
        {
            bestSuccessYet = success;
        }
        search.end();

    }

    if ( ! _stopReasons->checkTerminate() )
    {
        if (iterationSuccess)
        {
            AddOutputInfo("Search Successful. Enlarge Delta frame size.");
        }
        else
        {
            // 2. Poll
            NOMAD::Poll poll( this );
            poll.start();
            // Iteration is a success if either a better xFeas or
            // a better xInf (partial success or dominating) xInf was found.
            // See Algorithm 12.2 from DFBO.
            iterationSuccess = poll.run();

            NOMAD::SuccessType success = poll.getSuccessType();
            if (success > bestSuccessYet)
            {
                bestSuccessYet = success;
            }
            poll.end();
        }
    }

    setSuccessType(bestSuccessYet);

    // End of the iteration: iterationSuccess is true if we have a success (partial or full).
    return iterationSuccess;
}


bool NOMAD::MadsIteration::isMainIteration() const
{
    // This MadsIteration is the main iteration if it has the same mesh and k
    // as its parent MadsMegaIteration, and if the poll center is the first point of the MadsMegaIteration's barrier.
    bool ret = false;

    auto megaIter = dynamic_cast<const NOMAD::MadsMegaIteration*>(getParentOfType<NOMAD::MadsMegaIteration*>());
    if (nullptr != megaIter)
    {
        ret = (megaIter->getMesh() == _mesh && megaIter->getK() == _k);

        if (ret)
        {
            auto firstBarrierPoint = megaIter->getBarrier()->getFirstXFeas();
            if (nullptr == firstBarrierPoint)
            {
                firstBarrierPoint = megaIter->getBarrier()->getFirstXInf();
            }
            ret = (*_frameCenter == *firstBarrierPoint);
        }
    }

    return ret;
}

