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

#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/Mads/Poll.hpp"
#include "../../Math/RNG.hpp"

void NOMAD::Poll::init()
{
    _name = "Poll";
    verifyParentNotNull();
}


void NOMAD::Poll::startImp()
{

    // Create EvalPoints and send them to EvaluatorControl
    generateTrialPoints();

    if ( 0 == getTrialPointsCount() )
    {
        auto madsStopReasons = NOMAD::AlgoStopReasons<NOMAD::MadsStopType>::get ( _stopReasons );
        madsStopReasons->set( NOMAD::MadsStopType::MESH_PREC_REACHED );
    }
}


bool NOMAD::Poll::runImp()
{
    bool foundBetter = false;

    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);
    }


    return foundBetter;
}


void NOMAD::Poll::endImp()
{
    postProcessing(getEvalType());
}


// Generate poll directions
void NOMAD::Poll::setPollDirections(std::list<NOMAD::Direction> &directions) const
{
    directions.clear();
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");

    NOMAD::Direction dirUnit(n, 0.0);
    bool dirComputed = NOMAD::Direction::computeDirOnUnitSphere(dirUnit);
    if (!dirComputed)
    {
        std::string err("Poll: setPollDirections: Cannot compute a random direction on unit sphere");
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    // Ortho MADS 2n
    // Householder Matrix
    NOMAD::Direction** H = new NOMAD::Direction*[2*n];
    std::list<NOMAD::Direction> dirsUnit;

    // Ordering D_k alternates Hk and -Hk instead of [H_k -H_k]
    // TODO rewrite this in a clearer way, and follow the book's notation
    for (size_t i = 0; i < n; ++i)
    {
        dirsUnit.push_back(NOMAD::Direction(n, 0.0));
        H[i]   = &(*(--dirsUnit.end()));
        dirsUnit.push_back(NOMAD::Direction(n, 0.0));
        H[i+n] = &(*(--dirsUnit.end()));
    }

    // Householder transformations on the 2n directions on a unit n-sphere
    NOMAD::Direction::householder(dirUnit, true, H);
    delete [] H;

    // Scale the directions and project on the mesh
    std::list<NOMAD::Direction>::const_iterator itDir;
    int k = 1;
    std::shared_ptr<NOMAD::MeshBase> mesh = getIterationMesh();
    if (nullptr == mesh)
    {
        std::string err("Poll: setPollDirections: Iteration or Mesh not found.");
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    // TODO rewrite this in a clearer way
    for (itDir = dirsUnit.begin(); itDir != dirsUnit.end(); ++itDir, ++k)
    {
        directions.push_back(NOMAD::Direction(n, 0.0));
        NOMAD::Direction* pd = &(*(--directions.end()));

        // Compute infinite norm for direction pointed by itDir.
        NOMAD::Double infiniteNorm = (*itDir).infiniteNorm();
        if (0 == infiniteNorm)
        {
            std::string err("Poll: setPollDirections: Cannot handle an infinite norm of zero");
            throw NOMAD::Exception(__FILE__, __LINE__, err);
        }

        for (size_t i = 0; i < n; ++i)
        {
            // Scaling and projection on the mesh
            (*pd)[i] = mesh->scaleAndProjectOnMesh(i, (*itDir)[i] / infiniteNorm);
        }
    }

    std::list<NOMAD::Direction>::const_iterator it;
    for (it = directions.begin(); it != directions.end(); ++it)
    {
        AddOutputDebug("Poll direction: " + it->display());
    }

    dirsUnit.clear();
}


// Generate new points to evaluate
void NOMAD::Poll::generateTrialPoints()
{
    AddOutputInfo("Generate points for " + _name, true, false);
    
    // Creation of the poll directions
    std::list<NOMAD::Direction> directions;
    setPollDirections(directions);
    
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    
    // We need a frame center to start with.
    auto frameCenter = getIterationFrameCenter();
    if (nullptr == frameCenter || !frameCenter->ArrayOfDouble::isDefined() || frameCenter->size() != n)
    {
        std::string err("NOMAD::Poll::generateTrialPoints: invalid frame center: ");
        if (nullptr != frameCenter)
        {
            err += frameCenter->display();
        }
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }
    
    AddOutputDebug("Frame center: " + frameCenter->display());
    
    for (std::list<NOMAD::Direction>::iterator it = directions.begin(); it != directions.end() ; ++it)
    {
        NOMAD::Point pt(n);
        
        // pt = frame center + direction
        for (size_t i = 0 ; i < n ; ++i )
        {
            pt[i] = (*frameCenter)[i] + (*it)[i];
        }
        
        // Snap the points and the corresponding direction to the bounds
        if (snapPointToBoundsAndProjectOnMesh(pt, _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND"),
                        _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND"),
                        frameCenter,
                        getIterationMesh()))
        {
            if (pt != *frameCenter->getX())
            {
                // New EvalPoint to be evaluated.
                // Add it to the list.
                bool inserted = insertTrialPoint(NOMAD::EvalPoint(pt));
                
                std::string s = "Generated point";
                s += (inserted) ? ": " : " not inserted: ";
                s += pt.display();
                AddOutputInfo(s);
            }
        }
    }
    
    verifyPointsAreOnMesh(getName());
    updatePointsWithFrameCenter();
    
    AddOutputInfo("Generated " + NOMAD::itos(getTrialPointsCount()) + " points");
    AddOutputInfo("Generate points for " + _name, false, true);
    
}


