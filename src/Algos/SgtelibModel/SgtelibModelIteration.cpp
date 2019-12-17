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

#include <sstream>

#include "../../Algos/SgtelibModel/SgtelibModel.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelIteration.hpp"
#include "../../Algos/SgtelibModel/SgtelibModelUpdate.hpp"

void NOMAD::SgtelibModelIteration::init()
{
    _name = getAlgoName() + NOMAD::Iteration::getName();

    // Initialize optimize member - model optimizer on sgte
    auto modelAlgo = dynamic_cast<const NOMAD::SgtelibModel*>(getParentOfType<NOMAD::SgtelibModel*>());
    _optimize = std::make_shared<NOMAD::SgtelibModelOptimize>(modelAlgo,
                                                    _runParams, _pbParams);
}


void NOMAD::SgtelibModelIteration::startImp()
{
    // Model update
    // Use the cache to determine a sgtelib model
    // Update has a side effect of setting the _ready member.
    NOMAD::SgtelibModelUpdate update(this);
    update.start();
    update.run();
    update.end();
}


bool NOMAD::SgtelibModelIteration::runImp()
{
    //verifyGenerateAllPointsBeforeEval("Iteration::run()", false);

    bool optimizationOk = false;

    // Model Update is handled in start().

    auto modelAlgo = dynamic_cast<const NOMAD::SgtelibModel*>(getParentOfType<NOMAD::SgtelibModel*>());
    if (!_stopReasons->checkTerminate() && modelAlgo->isReady())
    {
        // Use the optimizer to find oracle points on this model
        // If mesh available, add mesh and frame size arguments.
        NOMAD::ArrayOfDouble initialMeshSize;
        NOMAD::ArrayOfDouble initialFrameSize;

        auto mesh = modelAlgo->getMesh();
        if (nullptr != mesh)
        {
            initialMeshSize = mesh->getdeltaMeshSize();
            initialFrameSize = mesh->getDeltaFrameSize();
        }

        // Setup Pb parameters just before optimization.
        // This way, we get the best X0s.
        _optimize->setupPbParameters(modelAlgo->getExtendedLowerBound(),
                                     modelAlgo->getExtendedUpperBound(),
                                     initialMeshSize,
                                     initialFrameSize);

        _optimize->start();
        optimizationOk = _optimize->run();
        _optimize->end();
    }

    // End of the iteration: return value of optimizationOk
    return optimizationOk;
}


// Oracle points are the best points found in sub optimization on sgte model.
const NOMAD::EvalPointSet& NOMAD::SgtelibModelIteration::getOraclePoints() const
{
    return _optimize->getOraclePoints();
}

