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

#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Mads/SgtelibSearchMethod.hpp"
#ifdef USE_SGTELIB
#include "../../Algos/SgtelibModel/SgtelibModel.hpp"
#endif
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"
//
// Reference: File Sgtelib_Model_Search.cpp in NOMAD 3.9.1
// Author: Bastien Talgorn

void NOMAD::SgtelibSearchMethod::init()
{
    setStepType(NOMAD::StepType::SEARCH_METHOD_SGTELIB_MODEL);
    verifyParentNotNull();

    const auto parentSearch = getParentStep()->getParentOfType<NOMAD::SgtelibSearchMethod*>(false);
    setEnabled((nullptr == parentSearch) && _runParams->getAttributeValue<bool>("SGTELIB_MODEL_SEARCH"));
#ifndef USE_SGTELIB
    if (isEnabled())
    {
        OUTPUT_INFO_START
        AddOutputInfo(getName() + " cannot be performed because NOMAD is compiled without sgtelib library");
        OUTPUT_INFO_END
        setEnabled(false);
    }
#endif

#ifdef USE_SGTELIB
    // Check that there is exactly one objective
    if (isEnabled())
    {
        const auto bbot = NOMAD::SgtelibModel::getBBOutputType();
        auto nbObj = NOMAD::getNbObj(bbot);
        if (0 == nbObj)
        {
            OUTPUT_INFO_START
            AddOutputInfo(getName() + " not performed when there is no objective function");
            OUTPUT_INFO_END
            setEnabled(false);
        }
        else if (nbObj > 1)
        {
            OUTPUT_INFO_START
            AddOutputInfo(getName() + " not performed on multi-objective function");
            OUTPUT_INFO_END
            setEnabled(false);
        }

        auto modelDisplay = _runParams->getAttributeValue<std::string>("SGTELIB_MODEL_DISPLAY");
        _displayLevel = modelDisplay.empty()
                            ? NOMAD::OutputLevel::LEVEL_DEBUGDEBUG
                            : NOMAD::OutputLevel::LEVEL_INFO;

        // Create the SgtelibModel algorithm with its own stop reasons
        auto stopReasons = std::make_shared<NOMAD::AlgoStopReasons<NOMAD::ModelStopType>>();
        auto barrier = getMegaIterationBarrier();
        const NOMAD::MadsIteration* iteration = getParentOfType<NOMAD::MadsIteration*>();
        auto mesh = iteration->getMesh();
        _modelAlgo = std::make_shared<NOMAD::SgtelibModel>(this, stopReasons,
                                                           barrier, _runParams,
                                                           _pbParams, mesh);
        _modelAlgo->setEndDisplay(false);
    }
#endif
}


bool NOMAD::SgtelibSearchMethod::runImp()
{
    // SgtelibModel algorithm _modelAlgo was created in init().
    // Use generateTrialPoints(). It will call createOraclePoints().
    generateTrialPoints();

    bool foundBetter = evalTrialPoints(this);

    return foundBetter;
}


void NOMAD::SgtelibSearchMethod::generateTrialPointsImp()
{
#ifdef USE_SGTELIB
    std::string s;
    NOMAD::EvalPointSet oraclePoints;

    // Get Iteration
    const NOMAD::MadsIteration* iteration = getParentOfType<NOMAD::MadsIteration*>();

    if (!_stopReasons->checkTerminate())
    {
        // Initial displays
        OUTPUT_INFO_START
        s = "Number of cache points: " + std::to_string(NOMAD::CacheBase::getInstance()->size());
        AddOutputInfo(s, _displayLevel);
        s = "Mesh size parameter: " + iteration->getMesh()->getdeltaMeshSize().display();
        AddOutputInfo(s, _displayLevel);
        NOMAD::OutputQueue::Flush();
        OUTPUT_INFO_END

        // Here, NOMAD 3 uses parameter SGTELIB_MODEL_TRIALS: Max number of
        // sgtelib model search failures before going to the poll step.
        // Not used.
        //const size_t kkmax = _runParams->getAttributeValue<size_t>("SGTELIB_MODEL_SEARCH_TRIALS");
        /*----------------*/
        /*  oracle points */
        /*----------------*/
        oraclePoints = _modelAlgo->createOraclePoints();

        if (0 == oraclePoints.size())
        {
            OUTPUT_INFO_START
            s = "Failed generating points. Stop " + getName();
            AddOutputInfo(s, _displayLevel);
            OUTPUT_INFO_END

            auto sgtelibModelStopReasons = NOMAD::AlgoStopReasons<NOMAD::ModelStopType>::get(_modelAlgo->getAllStopReasons());
            if (nullptr == sgtelibModelStopReasons)
            {
                throw NOMAD::Exception(__FILE__, __LINE__, "SgtelibModel Algorithm must have a Sgtelib stop reason");
            }
            sgtelibModelStopReasons->set(NOMAD::ModelStopType::ORACLE_FAIL);
        }
        else
        {
            // Add oracle point to _trialPoints.
            // SearchMethodBase will take care of projecting trial points to mesh.
            _trialPoints = oraclePoints;
        }
    }

#endif
}   // end generateTrialPoints


void NOMAD::SgtelibSearchMethod::getBestProjection(const NOMAD::Point& incumbent,
                                    const NOMAD::ArrayOfDouble& deltaMeshSize,
                                    std::shared_ptr<NOMAD::Point> x)
{
    // Issue #383: Use Projection class
}














