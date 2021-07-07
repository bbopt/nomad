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

#include "../../Algos/Mads/PollMethodBase.hpp"
#include "../../Algos/SubproblemManager.hpp"
#include "../../Output/OutputQueue.hpp"


void NOMAD::PollMethodBase::init()
{
    // A poll method must have a parent
    verifyParentNotNull();
}

void NOMAD::PollMethodBase::generateTrialPoints()
{
    generateTrialPointsInternal(false);
}

void NOMAD::PollMethodBase::generate2NDirections(std::list<NOMAD::Direction> &directions, size_t n) const
{
    NOMAD::Direction dirUnit(n, 0.0);
    NOMAD::Direction::computeDirOnUnitSphere(dirUnit);

    // Householder Matrix
    NOMAD::Direction** H = new NOMAD::Direction*[2*n];

    // Ordering D_k alternates Hk and -Hk instead of [H_k -H_k]
    for (size_t i = 0; i < n; ++i)
    {
        directions.push_back(NOMAD::Direction(n, 0.0));
        H[i]   = &(directions.back());
        directions.push_back(NOMAD::Direction(n, 0.0));
        H[i+n] = &(directions.back());
    }
    // Householder transformations on the 2n directions on a unit n-sphere
    NOMAD::Direction::householder(dirUnit, true, H);
    delete [] H;
}


void NOMAD::PollMethodBase::generateTrialPointsNPlus1(const NOMAD::EvalPointSet& inputTrialPoints)
{
    // Clear trial points generated in the first pass, where 2N points were
    // generated and N points were evaluated. The evaluations are not part of
    // this data.
    _trialPoints.clear();
    // Then add back the N evaluated points, with evaluation this time.
    for (auto trialPoint : inputTrialPoints)
    {
        insertTrialPoint(trialPoint);
    }

    generateTrialPointsInternal(true);
}


// If isNPlus1 is true, _trialPoints containes N evaluated trial points,
// and the N+1th point needs to be generated.
void NOMAD::PollMethodBase::generateTrialPointsInternal(const bool isNPlus1)
{

    std::list<NOMAD::Direction> directionsFullSpace = generateFullSpaceScaledDirections(isNPlus1);
    
    OUTPUT_INFO_START
    std::string s = "Generate ";
    s+= (isNPlus1) ? "n+1th point" : "n points";
    s += " for " + getName();
    AddOutputInfo(s, true, false);
    OUTPUT_INFO_END

    OUTPUT_DEBUG_START
    for (auto dir : directionsFullSpace)
    {
        AddOutputDebug("Scaled and mesh projected poll direction: " + dir.display());
    }
    OUTPUT_DEBUG_END

    auto n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    // We need a frame center to start with.
    if (!_frameCenter.ArrayOfDouble::isDefined() || _frameCenter.size() != n)
    {
        std::string err = "Invalid frame center: " + _frameCenter.display();
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    OUTPUT_DEBUG_START
    AddOutputDebug("Frame center: " + _frameCenter.display());
    OUTPUT_DEBUG_END

    if (isNPlus1)
    {
        // Clear trial points so that they are not re-evaluated (or tried to be).
        // If not N+1 pass, we keep all points from all PollMethods, for instance,
        // Ortho 2N for primary poll center and Double for secondary poll center.
        clearTrialPoints();
    }
    for (std::list<NOMAD::Direction>::iterator it = directionsFullSpace.begin(); it != directionsFullSpace.end() ; ++it)
    {
        NOMAD::Point pt(n);

        // pt = frame center + direction
        for (size_t i = 0 ; i < n ; ++i )
        {
            pt[i] = _frameCenter[i] + (*it)[i];
        }

        auto evalPoint = NOMAD::EvalPoint(pt);
        evalPoint.setPointFrom(std::make_shared<NOMAD::EvalPoint>(_frameCenter), NOMAD::SubproblemManager::getInstance()->getSubFixedVariable(this));
        // Snap the points and the corresponding direction to the bounds
        if (snapPointToBoundsAndProjectOnMesh(evalPoint,
                                              _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND"),
                                              _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND")))
        {
            if (*evalPoint.getX() != *_frameCenter.getX())
            {
                // New EvalPoint to be evaluated.
                // Add it to the list.
                evalPoint.addGenStep(getStepType());
                bool inserted = insertTrialPoint(evalPoint);

                OUTPUT_INFO_START
                std::string s = "Generated point";
                s += (inserted) ? ": " : " not inserted: ";
                s += evalPoint.display();
                AddOutputInfo(s);
                OUTPUT_INFO_END
            }
            else
            {
                OUTPUT_INFO_START
                std::string s = "Generated point not inserted (equal to frame center): ";
                s += evalPoint.display();
                AddOutputInfo(s);
                OUTPUT_INFO_END
            }
        }
    }

    OUTPUT_INFO_START
    AddOutputInfo("Generated " + NOMAD::itos(getTrialPointsCount()) + " points");
    std::string s = "Generate ";
    s+= (isNPlus1) ? "n+1th point" : "n points";
    s += " for " + getName();
    AddOutputInfo(s, false, true);
    OUTPUT_INFO_END
}


std::list<NOMAD::Direction> NOMAD::PollMethodBase::generateFullSpaceScaledDirections(bool isNPlus1, std::shared_ptr<NOMAD::MeshBase> mesh)
{
    // Groups of variables.
    auto varGroups = _pbParams->getAttributeValue<NOMAD::ListOfVariableGroup>("VARIABLE_GROUP");;
    
    auto n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    
    std::list<NOMAD::Direction> directionsSubSpace, directionsFullSpace;
    
    if (varGroups.size() == 0)
    {
        if (!isNPlus1)
        {
            // Creation of the poll directions in the full space
            generateUnitPollDirections(directionsFullSpace,n);
        }
        else
        {
            // Creation of the N+1th direction in the full space
            generateNPlus1Direction(directionsFullSpace);

            // The trial points were needed to generate the n+1th direction.
            // We do not want to re-evaluate them, and the n+1th point will
            // be added to _trialPoints when it is generated. So here, we clear the current trial points.
            _trialPoints.clear();
        }
    }
    else
    {
        for (auto vg : varGroups)
        {
            size_t nVG = vg.size();

            if (!isNPlus1)
            {
                // Creation of the poll directions in the sub space of the variable group
                generateUnitPollDirections(directionsSubSpace,nVG);
            }
            else
            {
                // Creation of the N+1th poll direction in the sub space of the variable group
                generateNPlus1Direction(directionsSubSpace);
                // We do not want to re-evaluate points generated in the first pass.
                _trialPoints.clear();
            }

            // Convert sub space (in a group of variable) directions to full space directions (all variables)
            if (varGroups.size() > 1)
            {
                size_t vgIndex = 0; // For Output debug only
                for (std::list<NOMAD::Direction>::iterator it = directionsSubSpace.begin(); it != directionsSubSpace.end() ; ++it)
                {
                    // In full space, the direction for an index outside the group of variables is null
                    NOMAD::Direction fullSpaceDirection(n,0.0);

                    // Copy the the sub space direction elements to full space
                    size_t i = 0;
                    for (auto index: vg)
                    {
                        fullSpaceDirection[index] = (*it)[i++];
                    }
                    directionsFullSpace.push_back(fullSpaceDirection);
                    OUTPUT_DEBUG_START
                    AddOutputDebug("Unit poll direction for Variable Group " + std::to_string(vgIndex) + ": "+ fullSpaceDirection.display());
                    OUTPUT_DEBUG_END
                    vgIndex++;
                }
            }
            else
            {
                directionsFullSpace = directionsSubSpace;

                OUTPUT_DEBUG_START
                for (auto dir : directionsFullSpace)
                {
                    AddOutputDebug("Unit poll direction: " + dir.display());
                }
                OUTPUT_DEBUG_END

            }
        }
    }
    
    // Scale and project directions on the mesh
    scaleAndProjectOnMesh(directionsFullSpace, mesh);
    
    return directionsFullSpace;
    
}


void NOMAD::PollMethodBase::scaleAndProjectOnMesh(std::list<Direction> & dirs, std::shared_ptr<NOMAD::MeshBase> mesh)
{
    // Scale the directions and project on the mesh
    if (nullptr == mesh)
    {
        mesh = getIterationMesh();
    }
    if (nullptr == mesh)
    {
        std::string err("Iteration or Mesh not found.");
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    std::list<NOMAD::Direction>::iterator itDir;
    auto n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    for (itDir = dirs.begin(); itDir != dirs.end(); ++itDir)
    {
        Direction scaledDir(n,0.0);

        // Compute infinite norm for direction pointed by itDir.
        NOMAD::Double infiniteNorm = (*itDir).infiniteNorm();
        if (0 == infiniteNorm)
        {
            std::string err("Cannot handle an infinite norm of zero");
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }

        for (size_t i = 0; i < n; ++i)
        {
            // Scaling and projection on the mesh
            scaledDir[i] = mesh->scaleAndProjectOnMesh(i, (*itDir)[i] / infiniteNorm);
        }

        *itDir = scaledDir;
    }
}
