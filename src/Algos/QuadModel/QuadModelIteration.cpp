/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4.0.0 has been created by                                      */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  The copyright of NOMAD - version 4.0.0 is owned by                             */
/*                 Charles Audet               - Polytechnique Montreal            */
/*                 Sebastien Le Digabel        - Polytechnique Montreal            */
/*                 Viviane Rochon Montplaisir  - Polytechnique Montreal            */
/*                 Christophe Tribes           - Polytechnique Montreal            */
/*                                                                                 */
/*  NOMAD v4 has been funded by Rio Tinto, Hydro-Québec, NSERC (Natural            */
/*  Sciences and Engineering Research Council of Canada), InnovÉÉ (Innovation      */
/*  en Énergie Électrique) and IVADO (The Institute for Data Valorization)         */
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

#include "../../Algos/QuadModel/QuadModelAlgo.hpp"
#include "../../Algos/QuadModel/QuadModelIteration.hpp"
#include "../../Algos/QuadModel/QuadModelOptimize.hpp"
#include "../../Algos/QuadModel/QuadModelUpdate.hpp"
#include "../../../ext/sgtelib/src/Surrogate_Factory.hpp"

void NOMAD::QuadModelIteration::reset()
{
    if (nullptr != _model)
    {
        _model.reset();
    }

    if (nullptr != _trainingSet)
    {
        _trainingSet.reset();
    }
}


void NOMAD::QuadModelIteration::init()
{
    _name = getAlgoName() + NOMAD::Iteration::getName();

    // Count the number of constraints
    const auto bbot = NOMAD::QuadModelAlgo::getBBOutputType();
    size_t nbConstraints = NOMAD::getNbConstraints(bbot);

    // Init the TrainingSet
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    SGTELIB::Matrix empty_X("empty_X", 0, static_cast<int>(n));
    SGTELIB::Matrix empty_Z("empty_Z", 0, static_cast<int>(nbConstraints+1));
    _trainingSet = std::make_shared<SGTELIB::TrainingSet>(empty_X, empty_Z);

    // The quadratic model uses Sgtelib
    _model = std::shared_ptr<SGTELIB::Surrogate>(SGTELIB::Surrogate_Factory(*_trainingSet, "TYPE PRS"));

}


void NOMAD::QuadModelIteration::startImp()
{
    incK();

    // Select the sample points to construct the model. Use a center pt and the cache
    NOMAD::QuadModelUpdate update(this);
    update.start();
    bool updateSuccess = update.run();
    update.end();

    if ( ! updateSuccess )
    {
        auto qmsStopReason = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get ( getAllStopReasons() );

        // The initial update is not a success. If the stop reason is not set to terminate we set a default stop reason for initialization.
        if ( !_stopReasons->checkTerminate() )
            qmsStopReason->set( NOMAD::ModelStopType::INITIAL_FAIL);
    }
}


bool NOMAD::QuadModelIteration::runImp()
{

    bool iterationSuccess = false;

    // Initialize optimize member - model optimizer on sgte
    NOMAD::QuadModelOptimize optimize (this, _pbParams);

    // Model Update is handled in start().
    if (!_stopReasons->checkTerminate() && _model->is_ready() )
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
