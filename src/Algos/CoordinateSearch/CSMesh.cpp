/*---------------------------------------------------------------------------------*/
/*  NOMAD - Nonlinear Optimization by Mesh Adaptive Direct Search -                */
/*                                                                                 */
/*  NOMAD - Version 4 has been created and developed by                            */
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
#include "../../Algos/CoordinateSearch/CSMesh.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::CSMesh::init()
{
    // Compute and initialize values for _frameSize.
    initFrameSizeGranular(_initialFrameSize);

    _initFrameSize.reset(_n);
    _initFrameSize = _frameSize;

    // Expecting _minMeshSize to be fully defined.
    if (!_minMeshSize.isComplete())
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Expecting mesh minimum size to be fully defined.");
    }
    
    for (size_t i = 0; i < _n; i++)
    {
        if ( _initialMeshSize[i] < _granularity[i])
            {
                throw NOMAD::Exception(__FILE__, __LINE__, "CSMesh: MeshSize below granularity ");
            }
    }
}


/*----------------------------------------------------------*/
/*  check the stopping condition on the mesh parameters     */
/*----------------------------------------------------------*/

void NOMAD::CSMesh::checkMeshForStopping( std::shared_ptr<NOMAD::AllStopReasons> stopReasons ) const
{
    bool stop = true;

    // CSMesh is the mesh for CoordinateSearch
    auto csStopReasons = NOMAD::AlgoStopReasons<NOMAD::CSStopType>::get( stopReasons );


    // General case, when min mesh size is reached, stop reason
    // MIN_MESH_SIZE_REACHED is set.
    // Special case: if all variables are granular, MAX_EVAL is automatically
    // set. Always continue looking, even if the min mesh size is reached,
    // but only for a finite number of iterations.
    
    bool allGranular = true;

    for (size_t i = 0; i < _n; i++)
    {
        if (0.0 == _granularity[i])
        {
            allGranular = false;
            break;
        }
    }
    if (allGranular)
    {
        stop = false;
    }
    else
    {
        for (size_t i = 0; i < _n; i++)
        {
            if (getdeltaMeshSize(i).todouble() >= _minMeshSize[i].todouble() && _granularity[i] ==0) // Do not stop if a non-granular variable has its mesh coarser than the minMeshSize
            {
                stop = false;
                break;
            }
        }
    }
    if (stop)
    {
        csStopReasons->set ( NOMAD::CSStopType::MIN_MESH_SIZE_REACHED );
    }

    if (!stop && !allGranular && _minFrameSize.isDefined())
    {
        // Reset stop
        stop = true;
        for (size_t i = 0; i < _n; i++)
        {
            if (_minFrameSize[i].isDefined() && getDeltaFrameSize(i).todouble() >= _minFrameSize[i].todouble())
            {
                stop = false;
                break;
            }
        }
        if (stop)
        {
            csStopReasons->set ( NOMAD::CSStopType::MIN_FRAME_SIZE_REACHED );
        }
    }
    for (size_t i = 0; i < _n; i++)
    {
        if (getDeltaFrameSize(i) < _granularity[i])
        {
            csStopReasons->set (NOMAD::CSStopType::GRANULARITY_REACHED);
        }
    }
}


// Update mesh size.
void NOMAD::CSMesh::updatedeltaMeshSize()
{
    // In CSMesh, delta is already updated, as a side effect of the
    // other mesh parameters being updated.
}

// In the CS algorithm, the mesh is never enlarged.
// Return value: always false because the mesh never changes.
bool NOMAD::CSMesh::enlargeDeltaFrameSize(const NOMAD::Direction& direction)
{
    return false;
}


// Update frame size after unsuccessful CS Poll steps.
// In CSMesh, big Delta (frame size) and small delta (mesh size) are
// updated simultaneously as a result of the implementation.
// All variables are refined (no anisotropy created)
void NOMAD::CSMesh::refineDeltaFrameSize()
{
    for (size_t i = 0 ; i < _n ; i++)
    {
        // Compute the new value frameSize first.
        // We will do some verifications before setting them.
        NOMAD::Double frameSize = _frameSize[i];
        refineDeltaFrameSize(frameSize, _granularity[i]);

        // Verify delta mesh size does not go too small if we use the new values.
        NOMAD::Double olddeltaMeshSize = getdeltaMeshSize(_frameSize[i], _granularity[i]);
        if (_minMeshSize[i] <= olddeltaMeshSize)
        {
            // We can go lower
            _frameSize[i] = frameSize;
        }
    }
}


void NOMAD::CSMesh::refineDeltaFrameSize(NOMAD::Double &frameSize,
                                        const NOMAD::Double& /*granularity*/) const
{
    frameSize *= 0.5;
}


// Initialize frame size  according to initial
// frame size and granularity.
void NOMAD::CSMesh::initFrameSizeGranular(const NOMAD::ArrayOfDouble &initialFrameSize)
{
    if (!initialFrameSize.isDefined() || initialFrameSize.size() != _n)
    {
        std::ostringstream oss;
        oss << "CSMesh: initFrameSizeGranular: inconsistent dimension of the frame size.";
        oss << " initial frame size defined: " << initialFrameSize.isDefined();
        oss << " size: " << initialFrameSize.size();
        oss << " n: " << _n;
        throw NOMAD::Exception(__FILE__, __LINE__, oss.str());
    }

    _frameSize.reset(_n);


    NOMAD::Double dMin;

    for (size_t i = 0 ; i < _n ; i++)
    {
        if (_granularity[i].todouble() > 0)
        {
            dMin = _granularity[i];
        }
        else
        {
            dMin = 1;
        }

        NOMAD::Double div = initialFrameSize[i] / dMin;
        double exp = std::log10(div.abs().todouble());
        _frameSize[i] =  pow(div.todouble() * pow(10 , -exp),exp);
    }

}

/*--------------------------------------------------------------*/
/*  get rho (ratio frame/mesh size)                             */
/*  The ratio is constant for CS Mesh                           */
/*--------------------------------------------------------------*/
NOMAD::Double NOMAD::CSMesh::getRho(size_t i) const // kept the function because declared as pure virtual in declaration meshbase
{
    NOMAD::Double rho = 2 ;

    return rho;
}


/*--------------------------------------------------------------*/
/*  get delta (mesh size parameter )                            */
/* delta mesh size is equal to 2* frame size parameter with CS  */
/*--------------------------------------------------------------*/
NOMAD::ArrayOfDouble NOMAD::CSMesh::getdeltaMeshSize() const
{
    return MeshBase::getdeltaMeshSize();
}


NOMAD::Double NOMAD::CSMesh::getdeltaMeshSize(size_t i) const
{
    NOMAD::Double deltai = getdeltaMeshSize(_frameSize[i], _granularity[i]);

    return deltai;
}


NOMAD::Double NOMAD::CSMesh::getdeltaMeshSize(const NOMAD::Double& frameSize,
                                              const NOMAD::Double& /*granularity*/) const
{
    NOMAD::Double delta = frameSize * 0.5;

    return delta;
}


/*--------------------------------------------------------------*/
/*  get Delta_i  (frame size parameter)                         */
/*       Delta^k = gran * frameSize                             */
/*--------------------------------------------------------------*/
NOMAD::ArrayOfDouble NOMAD::CSMesh::getDeltaFrameSize() const
{
    return MeshBase::getDeltaFrameSize();
}


NOMAD::Double NOMAD::CSMesh::getDeltaFrameSize(const size_t i) const
{
    return getDeltaFrameSize(_granularity[i], _frameSize[i]);
}


NOMAD::Double NOMAD::CSMesh::getDeltaFrameSize(const NOMAD::Double& granularity, const NOMAD::Double& frameSize) const
{
    NOMAD::Double dMinGran = 1.0;

    if (granularity > 0)
    {
        dMinGran = granularity;
    }

    NOMAD::Double Delta = dMinGran * frameSize;

    return Delta;
}


NOMAD::ArrayOfDouble NOMAD::CSMesh::getDeltaFrameSizeCoarser() const
{
    return MeshBase::getDeltaFrameSizeCoarser();
}


NOMAD::Double NOMAD::CSMesh::getDeltaFrameSizeCoarser(const size_t i) const
{
    NOMAD::Double frameSize = _frameSize[i];
    
    frameSize *=2;

    return getDeltaFrameSize(_granularity[i], frameSize);
}


// This method is used by the input operator>>
void NOMAD::CSMesh::setDeltas(const size_t i,
                              const NOMAD::Double &deltaMeshSize,
                              const NOMAD::Double &deltaFrameSize)
{
    
    // Value to use for granularity (division so default = 1.0)
    
    NOMAD::Double gran = 1.0;
    if (0.0 < _granularity[i])
    {
        _granularity[i] = gran;
    }

    _frameSize[i] = deltaFrameSize;
}


void NOMAD::CSMesh::setDeltas(const NOMAD::ArrayOfDouble &deltaMeshSize,
                             const NOMAD::ArrayOfDouble &deltaFrameSize)
{
    MeshBase::setDeltas(deltaMeshSize, deltaFrameSize);
}





/*-----------------------------------------------------------*/
/*     scale and project on the mesh /not changed yet        */
/*-----------------------------------------------------------*/
NOMAD::Double NOMAD::CSMesh::scaleAndProjectOnMesh(size_t i, const NOMAD::Double &l) const
{
    const NOMAD::Double delta = getdeltaMeshSize(i);

    if (i < _n && _frameSize.isDefined()  && delta.isDefined())
    {
        NOMAD::Double d = getRho(i) * l;
        return d.roundd() * delta;
    }
    else
    {
        std::ostringstream oss;
        oss << "CSMesh: scaleAndProjectOnMesh cannot be performed.";
        oss << " i = ";
        oss << i;
        oss << " frame size defined : " << _frameSize.isDefined();
        oss << " delta mesh size defined: " << delta.isDefined();
        throw NOMAD::Exception(__FILE__, __LINE__, oss.str());
    }
}


// Scale and project, using infinite norm.
NOMAD::ArrayOfDouble NOMAD::CSMesh::scaleAndProjectOnMesh(
    const NOMAD::Direction &dir) const
{
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    NOMAD::ArrayOfDouble proj(n);

    NOMAD::Double infiniteNorm = dir.infiniteNorm();

    if (0 == infiniteNorm)
    {
        std::string err("CSMesh: scaleAndProjectOnMesh: Cannot handle an infinite norm of zero");
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    for (size_t i = 0; i < n; ++i)
    {
        // Scaling and projection on the mesh
        proj[i] = this->scaleAndProjectOnMesh(i, dir[i] / infiniteNorm);
    }

    return proj;
}


NOMAD::Point NOMAD::CSMesh::projectOnMesh(const NOMAD::Point& point,
                                         const NOMAD::Point& frameCenter) const
{
    // Projection on the mesh
    const NOMAD::Point& proj = point;
    auto delta = getdeltaMeshSize();
    // To avoid running around in circles
    const size_t maxNbTry = 10;

    for (size_t i = 0; i < point.size(); ++i)
    {
        const NOMAD::Double& deltaI = delta[i];
        bool frameCenterIsOnMesh = (frameCenter[i].isMultipleOf(deltaI));

        NOMAD::Double diffProjFrameCenter = proj[i] - frameCenter[i];
        // Value which will be used in verifyPointIsOnMesh
        NOMAD::Double verifValueI = (frameCenterIsOnMesh) ? proj[i]
                                           : diffProjFrameCenter;

        // Force verifValueI to be a multiple of deltaI.
        // nbTry = 0 means point is already on mesh.
        // nbTry = 1 means the projection worked.
        // nbTry > 1 means the process went hacky by forcing the value to work
        // for verifyPointIsOnMesh.
        size_t nbTry = 0;   // Limit on the number of tries
        while (!verifValueI.isMultipleOf(deltaI) && nbTry <= maxNbTry)
        {
            NOMAD::Double newVerifValueI;

            if (0 == nbTry)
            {
                // Use closest projection
                NOMAD::Double vHigh = verifValueI.nextMult(deltaI);
                NOMAD::Double vLow = - (-verifValueI).nextMult(deltaI);
                NOMAD::Double diffHigh = vHigh - verifValueI;
                NOMAD::Double diffLow = verifValueI - vLow;
                verifValueI = (diffLow < diffHigh) ? vLow
                                                   : (diffHigh < diffLow) ? vHigh
                                                   : (proj[i] < 0) ? vLow : vHigh;
            }
            else
            {
                // Go hacky
                verifValueI = (diffProjFrameCenter >= 0) ? verifValueI.nextMult(deltaI)
                                                         : - (-verifValueI).nextMult(deltaI);
            }

            proj[i] = (frameCenterIsOnMesh) ? verifValueI
                                            : verifValueI + frameCenter[i];

            // Recompute verifValue for more precision
            newVerifValueI = (frameCenterIsOnMesh) ? proj[i]
                                                   : proj[i] - frameCenter[i];

            nbTry++;

            // Special cases
            while (newVerifValueI != verifValueI && nbTry <= maxNbTry)
            {
                if (verifValueI >= 0)
                {
                    verifValueI = NOMAD::max(verifValueI, newVerifValueI);
                    verifValueI += NOMAD::DEFAULT_EPSILON;
                    verifValueI = verifValueI.nextMult(deltaI);
                }
                else
                {
                    verifValueI = NOMAD::min(verifValueI, newVerifValueI);
                    verifValueI -= NOMAD::DEFAULT_EPSILON;
                    verifValueI = - (-verifValueI).nextMult(deltaI);
                }

                proj[i] = (frameCenterIsOnMesh) ? verifValueI
                                               : verifValueI + frameCenter[i];

                 // Recompute verifValue for more precision
                newVerifValueI = (frameCenterIsOnMesh) ? proj[i]
                                                       : proj[i] - frameCenter[i];

                nbTry++;
            }

            verifValueI = newVerifValueI;
        }

        if (nbTry >= maxNbTry && !verifValueI.isMultipleOf(deltaI))
        {
            // Some values are just ill-conditioned.
            std::string s = "Warning: Could not project point (index " + std::to_string(i) + ") ";
            s += point.display() + " on mesh " + delta.display();
            s += " with frame center " + frameCenter.display();
            NOMAD::OutputInfo outputInfo("Mesh", s, NOMAD::OutputLevel::LEVEL_INFO);
            NOMAD::OutputQueue::Add(std::move(outputInfo));

            // Revert proj to its original value,
            // since the hack did not work.
            proj[i] = point[i];
        }
    }
    return proj;
}



