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
/**
 \file   Subproblem.cpp
 \brief  Subproblem of lesser dimension than the original problem
 \author Viviane Rochon Montplaisir
 \date   February 2019
 */
#include "Subproblem.hpp"
#include "../Output/OutputQueue.hpp"


/*------------*/
/* Subproblem */
/*------------*/
void NOMAD::Subproblem::init()
{
    if (nullptr == _refPbParams)
    {
        throw NOMAD::Exception(__FILE__, __LINE__,
                               "A valid PbParameters must be provided to the Subproblem constructor.");
    }

    // Compute new dimension
    auto nbFixed = _fixedVariable.nbDefined();
    const size_t refDimension = _refPbParams->getAttributeValue<size_t>("DIMENSION");
    _dimension = refDimension - nbFixed;

    std::string s = "FIXED_VARIABLE set to " + _fixedVariable.display();
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_INFO);

    setupProblemParameters();
}


NOMAD::Subproblem::~Subproblem()
{
}


void NOMAD::Subproblem::setupProblemParameters()
{
    // NOTE: If a new parameter with dimension is added to the class PbParameters,
    // this method will break.
    // It could be generalized by going through each parameter, and adjust it only
    // if it is an ArrayOfDouble, Point, or Dimension. Backlog task.
    const size_t refDimension = _refPbParams->getAttributeValue<size_t>("DIMENSION");
    size_t n = _dimension;

    _subPbParams = std::make_shared<NOMAD::PbParameters>(*_refPbParams);
    _subPbParams->setAttributeValue("DIMENSION", n);

    // Get reference values for parameters
    const auto refX0s           = _refPbParams->getAttributeValue<NOMAD::ArrayOfPoint>("X0");
    const auto refLowerBound    = _refPbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND");
    const auto refUpperBound    = _refPbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND");
    const auto refBBInputType   = _refPbParams->getAttributeValue<NOMAD::BBInputTypeList>("BB_INPUT_TYPE");
    const auto refInitMeshSize  = _refPbParams->getAttributeValue<NOMAD::ArrayOfDouble>("INITIAL_MESH_SIZE");
    const auto refInitFrameSize  = _refPbParams->getAttributeValue<NOMAD::ArrayOfDouble>("INITIAL_FRAME_SIZE");
    const auto refMinMeshSize   = _refPbParams->getAttributeValue<NOMAD::ArrayOfDouble>("MIN_MESH_SIZE");
    const auto refMinFrameSize   = _refPbParams->getAttributeValue<NOMAD::ArrayOfDouble>("MIN_FRAME_SIZE");
    const auto refGranularity   = _refPbParams->getAttributeValue<NOMAD::ArrayOfDouble>("GRANULARITY");

    // Initialize new arrays
    NOMAD::ArrayOfPoint x0s;
    for (size_t x0index = 0; x0index < refX0s.size(); x0index++)
    {
        NOMAD::Point x0(n);
        x0s.push_back(x0);
    }
    NOMAD::Point fixedVariable(n);
    NOMAD::ArrayOfDouble lb(n), ub(n);
    NOMAD::BBInputTypeList bbInputType;
    NOMAD::ArrayOfDouble initialMeshSize(n), initialFrameSize(n), minMeshSize(n), minFrameSize(n);
    NOMAD::ArrayOfDouble granularity(n);

    // Compute new values, simply using the values on the positions of non-fixed variables.
    size_t i = 0;
    for (size_t refIndex = 0; refIndex < refDimension; refIndex++)
    {
        if (_fixedVariable[refIndex].isDefined())
        {
            continue;
        }

        for (size_t x0index = 0; x0index < refX0s.size(); x0index++)
        {
            auto refX0 = refX0s[x0index];
            x0s[x0index][i] = refX0[refIndex];
        }
        lb[i] = refLowerBound[refIndex];
        ub[i] = refUpperBound[refIndex];
        bbInputType.push_back(refBBInputType[refIndex]);
        initialMeshSize[i] = refInitMeshSize[refIndex];
        initialFrameSize[i] = refInitFrameSize[refIndex];
        minMeshSize[i] = refMinMeshSize[refIndex];
        minFrameSize[i] = refMinFrameSize[refIndex];
        granularity[i] = refGranularity[refIndex];

        i++;
    }

    // Set new values to _subPbParams
    _subPbParams->setAttributeValue("X0", x0s);
    _subPbParams->setAttributeValue("FIXED_VARIABLE", fixedVariable);
    _subPbParams->setAttributeValue("LOWER_BOUND", lb);
    _subPbParams->setAttributeValue("UPPER_BOUND", ub);
    _subPbParams->setAttributeValue("BB_INPUT_TYPE", bbInputType);
    _subPbParams->setAttributeValue("INITIAL_MESH_SIZE", initialMeshSize);
    _subPbParams->setAttributeValue("INITIAL_FRAME_SIZE", initialFrameSize);
    _subPbParams->setAttributeValue("MIN_MESH_SIZE", minMeshSize);
    _subPbParams->setAttributeValue("MIN_FRAME_SIZE", minFrameSize);
    _subPbParams->setAttributeValue("GRANULARITY", granularity);

    
    _subPbParams->checkAndComply();

}
