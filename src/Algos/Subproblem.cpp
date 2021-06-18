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
/**
 \file   Subproblem.cpp
 \brief  Subproblem of lesser dimension than the original problem
 \author Viviane Rochon Montplaisir
 \date   February 2019
 */
#include "../Algos/Subproblem.hpp"
#include "../Output/OutputQueue.hpp"
#include "../Type/BBInputType.hpp"


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

    if (_fixedVariable.isEmpty())
    {
        std::string s = "Error: Fixed variable of dimension 0";
        throw NOMAD::Exception(__FILE__,__LINE__,s);
    }

    // Compute new dimension
    NOMAD::Point subFixedVariable = _refPbParams->getAttributeValue<NOMAD::Point>("FIXED_VARIABLE");
    _dimension = subFixedVariable.size() - subFixedVariable.nbDefined();


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
    const auto refListVariableGroup = _refPbParams->getAttributeValue<NOMAD::ListOfVariableGroup>("VARIABLE_GROUP");

    // Initialize new arrays
    NOMAD::ArrayOfPoint x0s;
    for (size_t x0index = 0; x0index < refX0s.size(); x0index++)
    {
        NOMAD::Point x0(n);
        x0s.push_back(x0);
    }
    NOMAD::Point subFixedVariable(n);
    NOMAD::ArrayOfDouble lb(n), ub(n);
    NOMAD::BBInputTypeList bbInputType;
    NOMAD::ArrayOfDouble initialMeshSize(n), initialFrameSize(n), minMeshSize(n), minFrameSize(n);
    NOMAD::ArrayOfDouble granularity(n);
    NOMAD::ListOfVariableGroup listVariableGroup = refListVariableGroup ;

    // Compute new fixed variable.
    // Current value of _fixedVariable contains only values from parent. Merge in values from _refPbParams.
    NOMAD::Point refFixedVariable = _refPbParams->getAttributeValue<NOMAD::Point>("FIXED_VARIABLE");

    // Compute new values, simply using the values on the positions of non-fixed variables.
    size_t i = 0;
    for (size_t refIndex = 0; refIndex < refDimension; refIndex++)
    {
        if (refFixedVariable[refIndex].isDefined())
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
    resetVariableGroupsAgainstFixedVariables(listVariableGroup, refFixedVariable);


    // Set new values to _subPbParams
    _subPbParams->setAttributeValue("X0", x0s);
    _subPbParams->setAttributeValue("FIXED_VARIABLE", subFixedVariable);    // No fixed variable defined in the subproblem
    _subPbParams->setAttributeValue("LOWER_BOUND", lb);
    _subPbParams->setAttributeValue("UPPER_BOUND", ub);
    _subPbParams->setAttributeValue("BB_INPUT_TYPE", bbInputType);
    _subPbParams->setAttributeValue("INITIAL_MESH_SIZE", initialMeshSize);
    _subPbParams->setAttributeValue("INITIAL_FRAME_SIZE", initialFrameSize);
    _subPbParams->setAttributeValue("MIN_MESH_SIZE", minMeshSize);
    _subPbParams->setAttributeValue("MIN_FRAME_SIZE", minFrameSize);
    _subPbParams->setAttributeValue("GRANULARITY", granularity);
    _subPbParams->setAttributeValue("VARIABLE_GROUP", listVariableGroup);

    _subPbParams->doNotShowWarnings();

    _subPbParams->checkAndComply();

    // Only now, conbine refFixedVariable into _fixedVariable.
    // Verify that refFixedVariable is of the dimension of the subproblem defined by _fixedVariable.
    if (refFixedVariable.size() == _fixedVariable.size())
    {
        refFixedVariable = refFixedVariable.makeSubSpacePointFromFixed(_fixedVariable);
    }
    const size_t subdim1 = _fixedVariable.size() - _fixedVariable.nbDefined();
    if (refFixedVariable.size() != subdim1)
    {
        std::string s = "Expecting FIXED_VARIABLE to be of size "  + std::to_string(subdim1);
        s += ". Current FIXED_VARIABLE is of size " + std::to_string(refFixedVariable.size());
        s += ": " + refFixedVariable.display();
        throw NOMAD::Exception(__FILE__,__LINE__, s);
    }

    for (size_t refIndex = 0, newIndex = 0; refIndex < _fixedVariable.size(); refIndex++)
    {
        if (!_fixedVariable[refIndex].isDefined())
        {
            _fixedVariable[refIndex] = refFixedVariable[newIndex];
            newIndex++;
        }
    }

    OUTPUT_INFO_START
    std::string s = "FIXED_VARIABLE set to " + _fixedVariable.display();
    NOMAD::OutputQueue::Add(s, NOMAD::OutputLevel::LEVEL_INFO);
    OUTPUT_INFO_END

}


// If a variable is fixed, its index must be removed from the group of variables
// All the indices above the fixed variable index must be decreased by one
void NOMAD::Subproblem::resetVariableGroupsAgainstFixedVariables(NOMAD::ListOfVariableGroup & lvg, const NOMAD::Point & fixedVar) const
{
    if (lvg.empty() || !fixedVar.isDefined())
    {
        return;
    }

    // Put the indices of fixed variables in a single set.
    const size_t n = fixedVar.size();
    std::set<size_t> indicesToRemove;
    for (size_t i = 0 ; i < n ; i++)
    {
        if (fixedVar[i].isDefined())
        {
            indicesToRemove.insert(i);
        }
    }

    NOMAD::ListOfVariableGroup updatedLvg;
    while (!indicesToRemove.empty())
    {
        updatedLvg.clear();
        auto itIndexBegin = indicesToRemove.begin();

        // Remove an index in a variable group. Decrement by one all indices (in all variable groups) above the index to remove.
        for (auto vg: lvg)
        {
            NOMAD::VariableGroup updatedVariableGroup;
            for (auto index : vg)
            {
                // Do not include an index equal to the index to remove.
                // Include and decrement indices above the index to remove.
                // Include indices below the index to remove.
                if (index > *itIndexBegin)
                {
                    updatedVariableGroup.insert(index-1);
                }
                else if (index < *itIndexBegin)
                {
                    updatedVariableGroup.insert(index);
                }

            }
            // A variable group can be empty -> do not include.
            if (!updatedVariableGroup.empty())
            {
                updatedLvg.push_back(updatedVariableGroup);
            }
        }

        // Remove index from the set of indices. Decrement remaining indices that are smaller than the index to remove.
        std::set<size_t> updatedIndicesToRemove;
        for (std::set<size_t>::iterator itIndex = ++itIndexBegin; itIndex != indicesToRemove.end() ; ++itIndex)
        {
            updatedIndicesToRemove.insert((*itIndex)-1);
        }
        indicesToRemove = updatedIndicesToRemove;
        lvg = updatedLvg;
    }

}

