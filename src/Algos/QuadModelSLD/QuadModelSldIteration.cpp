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

#include "../../Algos/QuadModelSLD/QuadModelSldAlgo.hpp"
#include "../../Algos/QuadModelSLD/QuadModelSldIteration.hpp"
#include "../../Algos/QuadModelSLD/QuadModelSldOptimize.hpp"
#include "../../Algos/QuadModelSLD/QuadModelSldUpdate.hpp"


void NOMAD::QuadModelSldIteration::init()
{

    // Count the number of constraints
    const auto bbot = NOMAD::QuadModelSldAlgo::getBBOutputType();

    // Init the TrainingSet
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");

    // The quadratic model
    _model = std::shared_ptr<NOMAD::QuadModelSld>(new NOMAD::QuadModelSld(bbot ,n));
    
    if (_trialPoints.size() > 0)
    {
        _useForSortingTrialPoints = true;
        setStepType(NOMAD::StepType::QUAD_MODEL_SORT);
    }

}

std::string NOMAD::QuadModelSldIteration::getName() const
{
    if (_useForSortingTrialPoints)
    {
        return NOMAD::stepTypeToString(_stepType) + " #" + std::to_string(_k);
    }
    else
    {
        return NOMAD::Iteration::getName();
    }
}


void NOMAD::QuadModelSldIteration::startImp()
{

    // Select the sample points to construct the model. Use a center pt and the cache
    
    NOMAD::QuadModelSldUpdate update(this, _trialPoints);
    update.start();
    bool updateSuccess = update.run();
    update.end();

    if ( ! updateSuccess && ! _useForSortingTrialPoints)
    {
        auto qmsStopReason = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get ( getAllStopReasons() );

        // The initial update is not a success. If the stop reason is not set to terminate we set a default stop reason for initialization.
        if ( !_stopReasons->checkTerminate() )
            qmsStopReason->set( NOMAD::ModelStopType::INITIAL_FAIL);
    }
}


bool NOMAD::QuadModelSldIteration::runImp()
{

    bool iterationSuccess = false;

    // Initialize optimize member on model
    NOMAD::QuadModelSldOptimize optimize (this, _pbParams);

    // Model Update is handled in start().
    if (!_stopReasons->checkTerminate() && _model->check() )
    {
        // Optimize to find oracle points on this model
        optimize.start();
        iterationSuccess = optimize.run();
        optimize.end();
    }

    // Update MegaIteration success type
    NOMAD::SuccessType success = optimize.getSuccessType();
    auto megaIter = getParentOfType<NOMAD::MegaIteration*>();
    megaIter->setSuccessType(success);

    // End of the iteration: iterationSuccess is true if we have a success.
    return iterationSuccess;

}
