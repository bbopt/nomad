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
#include "../../Algos/Mads/GMesh.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::GMesh::init()
{
    // Compute and initialize values for _frameSizeMant and _frameSizeExp.
    initFrameSizeGranular(_initialFrameSize);

    _initFrameSizeExp.reset(_n);
    _initFrameSizeExp = _frameSizeExp;

    // Expecting _minMeshSize to be fully defined.
    if (!_minMeshSize.isComplete())
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Expecting mesh minimum size to be fully defined.");
    }

    // Sanity checks
    if (_enforceSanityChecks)
    {
        for (size_t i = 0 ; i < _n ; i++)
        {
            checkFrameSizeIntegrity(_frameSizeExp[i], _frameSizeMant[i]);
            checkDeltasGranularity(i, getdeltaMeshSize(i), getDeltaFrameSize(i));
        }
    }
}


/*----------------------------------------------------------*/
/*  check the stopping condition on the mesh parameters     */
/*----------------------------------------------------------*/

void NOMAD::GMesh::checkMeshForStopping( std::shared_ptr<NOMAD::AllStopReasons> stopReasons ) const
{
    bool stop = true;

    // GMesh is the mesh for Mads
    // stopReasons must be AlgoStopReasons<MadsStopType>
    auto madsStopReasons = NOMAD::AlgoStopReasons<NOMAD::MadsStopType>::get( stopReasons );


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
        madsStopReasons->set ( NOMAD::MadsStopType::MIN_MESH_SIZE_REACHED );
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
            madsStopReasons->set ( NOMAD::MadsStopType::MIN_FRAME_SIZE_REACHED );
        }
    }

}


// Update mesh size.
void NOMAD::GMesh::updatedeltaMeshSize()
{
    // In GMesh, delta is already updated, as a side effect of the
    // other mesh parameters being updated.
}

// Update frame size after a successful Search or Frame step.
// In GMesh, big Delta and small delta are updated simultaneously as a
// result of the implementation.
// Return value: true if the mesh changed, false otherwise.
bool NOMAD::GMesh::enlargeDeltaFrameSize(const NOMAD::Direction& direction,
                                         const NOMAD::Double& anisotropyFactor,
                                         bool anisotropicMesh)
{
    bool oneFrameSizeChanged = false;

    NOMAD::Double minRho = NOMAD::INF;
    for (size_t i = 0; i < _n ; i++)
    {
        if (0 == _granularity[i])
        {
            minRho = NOMAD::min(minRho, getRho(i));
        }
    }

    for (size_t i = 0 ; i < _n ; i++)
    {

        bool frameSizeIChanged = false;
        // Test for producing anisotropic mesh
        if ( ! anisotropicMesh
            || ( direction[i].abs() / getdeltaMeshSize(i) / getRho(i) > anisotropyFactor )
            || ( _granularity[i] == 0  && _frameSizeExp[i] < _initFrameSizeExp[i] && getRho(i) > minRho*minRho )
            )
        {
            getLargerMantExp(_frameSizeMant[i], _frameSizeExp[i]);
            frameSizeIChanged = true;
            oneFrameSizeChanged = true;
        }

        // Sanity checks
        if (_enforceSanityChecks && frameSizeIChanged)
        {
            checkFrameSizeIntegrity(_frameSizeExp[i], _frameSizeMant[i]);
            checkDeltasGranularity(i, getdeltaMeshSize(i), getDeltaFrameSize(i));
        }
    }

    return oneFrameSizeChanged;
}


// Update frame size after unsuccessful Search and Poll steps.
// In GMesh, big Delta (frame size) and small delta (mesh size) are
// updated simultaneously as a result of the implementation.
void NOMAD::GMesh::refineDeltaFrameSize()
{
    for (size_t i = 0 ; i < _n ; i++)
    {
        // Compute the new values frameSizeMant and frameSizeExp first.
        // We will do some verifications before setting them.
        NOMAD::Double frameSizeMant = _frameSizeMant[i];
        NOMAD::Double frameSizeExp  = _frameSizeExp[i];
        refineDeltaFrameSize(frameSizeMant, frameSizeExp, _granularity[i]);

        // Verify delta mesh size does not go too small if we use the new values.
        NOMAD::Double olddeltaMeshSize = getdeltaMeshSize(_frameSizeExp[i], _initFrameSizeExp[i], _granularity[i]);
        if (_minMeshSize[i] <= olddeltaMeshSize)
        {
            // We can go lower
            _frameSizeMant[i] = frameSizeMant;
            _frameSizeExp[i] = frameSizeExp;
        }

        // Sanity checks
        if (_enforceSanityChecks)
        {
            checkFrameSizeIntegrity(_frameSizeExp[i], _frameSizeMant[i]);
            checkDeltasGranularity(i, getdeltaMeshSize(i), getDeltaFrameSize(i));
        }
    }
}


void NOMAD::GMesh::refineDeltaFrameSize(NOMAD::Double &frameSizeMant,
                                        NOMAD::Double &frameSizeExp,
                                        const NOMAD::Double &granularity) const
{
    if (frameSizeMant == 1)
    {
        frameSizeMant = 5;
        --frameSizeExp;
    }
    else if (frameSizeMant == 2)
    {
        frameSizeMant = 1;
    }
    else
    {
        frameSizeMant = 2;
    }

    if (granularity > 0 && frameSizeExp <= 0 && frameSizeMant != 1)
    {
        frameSizeExp = 0;
        frameSizeMant = 1;
    }
}


// Initialize frame size mantisse and exponent according to initial
// frame size and granularity.
void NOMAD::GMesh::initFrameSizeGranular(const NOMAD::ArrayOfDouble &initialFrameSize)
{
    if (!initialFrameSize.isDefined() || initialFrameSize.size() != _n)
    {
        std::ostringstream oss;
        oss << "GMesh: initFrameSizeGranular: inconsistent dimension of the frame size.";
        oss << " initial frame size defined: " << initialFrameSize.isDefined();
        oss << " size: " << initialFrameSize.size();
        oss << " n: " << _n;
        throw NOMAD::Exception(__FILE__, __LINE__, oss.str());
    }

    _frameSizeExp.reset(_n);
    _frameSizeMant.reset(_n);

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
        int exp = roundFrameSizeExp(std::log10(div.abs().todouble()));
        _frameSizeExp[i]  = exp;
        _frameSizeMant[i] = roundFrameSizeMant(div.todouble() * pow(10 , -exp));
    }

}


int NOMAD::GMesh::roundFrameSizeExp(const NOMAD::Double &exp)
{
    int frameSizeExp = (int)(exp.todouble());

    return frameSizeExp;
}


int NOMAD::GMesh::roundFrameSizeMant(const NOMAD::Double &mant)
{
    // Round input mant to 1, 2 or 5 for _frameSizeMant element.
    int frameSizeMant = 0;

    if (mant < 1.5)
    {
        frameSizeMant = 1;
    }
    else if ( mant >= 1.5 && mant < 3.5 )
    {
        frameSizeMant = 2;
    }
    else
    {
        frameSizeMant = 5;
    }

    return frameSizeMant;
}


/*--------------------------------------------------------------*/
/*  get rho (ratio frame/mesh size)                              */
/*       rho^k = 10^(b^k-|b^k-b_0^k|)                           */
/*--------------------------------------------------------------*/
NOMAD::Double NOMAD::GMesh::getRho(size_t i) const
{
    NOMAD::Double rho;
    NOMAD::Double diff = _frameSizeExp[i] - _initFrameSizeExp[i];
    NOMAD::Double powDiff = pow(10.0, diff.abs().todouble());
    if (_granularity[i] > 0)
    {
        rho = _frameSizeMant[i] * NOMAD::min(pow(10.0, _frameSizeExp[i].todouble()), powDiff);
    }
    else
    {
        rho = _frameSizeMant[i] * powDiff;
    }

    return rho;
}


/*--------------------------------------------------------------*/
/*  get delta (mesh size parameter)                             */
/*       delta^k = 10^(b^k-|b^k-b_0^k|)                         */
/*       if (granularity > 0)                                   */
/*           delta^k = granularity * max (1.0, delta^k )        */
/*--------------------------------------------------------------*/
NOMAD::ArrayOfDouble NOMAD::GMesh::getdeltaMeshSize() const
{
    return MeshBase::getdeltaMeshSize();
}


NOMAD::Double NOMAD::GMesh::getdeltaMeshSize(size_t i) const
{
    NOMAD::Double deltai = getdeltaMeshSize(_frameSizeExp[i], _initFrameSizeExp[i], _granularity[i]);

    return deltai;
}


NOMAD::Double NOMAD::GMesh::getdeltaMeshSize(const NOMAD::Double& frameSizeExp,
                                             const NOMAD::Double& initFrameSizeExp,
                                             const NOMAD::Double& granularity) const
{
    NOMAD::Double diff  = frameSizeExp - initFrameSizeExp;
    NOMAD::Double exp   = frameSizeExp - diff.abs();
    NOMAD::Double delta = pow(10.0, exp.todouble());

    if (0.0 < granularity)
    {
        // Ensure delta is a multiple of granularity
        delta = granularity * max(1.0 , delta);
    }

    return delta;
}


/*--------------------------------------------------------------*/
/*  get Delta_i  (frame size parameter)                          */
/*       Delta^k = a^k *10^{b^k}                                */
/*--------------------------------------------------------------*/
NOMAD::ArrayOfDouble NOMAD::GMesh::getDeltaFrameSize() const
{
    return MeshBase::getDeltaFrameSize();
}


NOMAD::Double NOMAD::GMesh::getDeltaFrameSize(const size_t i) const
{
    return getDeltaFrameSize(_granularity[i], _frameSizeMant[i], _frameSizeExp[i]);
}


NOMAD::Double NOMAD::GMesh::getDeltaFrameSize(const NOMAD::Double& granularity, const NOMAD::Double& frameSizeMant, const NOMAD::Double& frameSizeExp) const
{
    NOMAD::Double dMinGran = 1.0;

    if (granularity > 0)
    {
        dMinGran = granularity;
    }

    NOMAD::Double Delta = dMinGran * frameSizeMant * pow(10, frameSizeExp.todouble());

    return Delta;
}


NOMAD::ArrayOfDouble NOMAD::GMesh::getDeltaFrameSizeCoarser() const
{
    return MeshBase::getDeltaFrameSizeCoarser();
}


NOMAD::Double NOMAD::GMesh::getDeltaFrameSizeCoarser(const size_t i) const
{
    NOMAD::Double frameSizeMant = _frameSizeMant[i];
    NOMAD::Double frameSizeExp  = _frameSizeExp[i];
    getLargerMantExp(frameSizeMant, frameSizeExp);

    return getDeltaFrameSize(_granularity[i], frameSizeMant, frameSizeExp);
}


// This method is used by the input operator>>
void NOMAD::GMesh::setDeltas(const size_t i,
                             const NOMAD::Double &deltaMeshSize,
                             const NOMAD::Double &deltaFrameSize)
{
    // Input checks
    checkDeltasGranularity(i, deltaMeshSize, deltaFrameSize);

    // Value to use for granularity (division so default = 1.0)
    NOMAD::Double gran = 1.0;
    if (0.0 < _granularity[i])
    {
        gran = _granularity[i];
    }

    NOMAD::Double mant;
    NOMAD::Double exp;

    // Compute mantisse first
    // There are only 3 cases: 1, 2, 5, so compute all
    // 3 possibilities and then assign the values that work.
    NOMAD::Double mant1 = deltaFrameSize / (1.0 * gran);
    NOMAD::Double mant2 = deltaFrameSize / (2.0 * gran);
    NOMAD::Double mant5 = deltaFrameSize / (5.0 * gran);
    NOMAD::Double exp1 = std::log10(mant1.todouble());
    NOMAD::Double exp2 = std::log10(mant2.todouble());
    NOMAD::Double exp5 = std::log10(mant5.todouble());

    // deltaFrameSize = gran * mant * 10^exp  (where gran is 1.0 if granularity is not defined)
    // => exp = log10(deltaFrameSize / (mant * gran))
    // exp must be an integer so verify which one of the 3 values exp1, exp2, exp5
    // is an integer and use that value for exp, and the corresponding value
    // 1, 2 or 5 for mant.
    if (exp1.isInteger())
    {
        mant = 1;
        exp = exp1;
    }
    else if (exp2.isInteger())
    {
        mant = 2;
        exp = exp2;
    }
    else
    {
        mant = 5;
        exp = exp5;
    }
    _frameSizeExp[i] = roundFrameSizeExp(exp);
    _frameSizeMant[i] = mant;

    // Sanity checks
    if (_enforceSanityChecks)
    {
        checkFrameSizeIntegrity(_frameSizeExp[i], _frameSizeMant[i]);
        checkSetDeltas(i, deltaMeshSize, deltaFrameSize);
        checkDeltasGranularity(i, getdeltaMeshSize(i), getDeltaFrameSize(i));
    }
}


void NOMAD::GMesh::setDeltas(const NOMAD::ArrayOfDouble &deltaMeshSize,
                             const NOMAD::ArrayOfDouble &deltaFrameSize)
{
    MeshBase::setDeltas(deltaMeshSize, deltaFrameSize);
}


void NOMAD::GMesh::checkDeltasGranularity(const size_t i,
                                          const NOMAD::Double &deltaMeshSize,
                                          const NOMAD::Double &deltaFrameSize) const
{
    // Sanity checks
    if (_granularity[i] > 0.0)
    {
        bool hasError = false;
        std::string err = "Error: setDeltas: ";
        if (!deltaMeshSize.isMultipleOf(_granularity[i]))
        {
            hasError = true;
            err += "deltaMeshSize at index " + std::to_string(i);
            err += " is not a multiple of granularity " + _granularity[i].tostring();
        }
        else if (!deltaFrameSize.isMultipleOf(_granularity[i]))
        {
            hasError = true;
            err += "deltaFrameSize at index " + std::to_string(i);
            err += " is not a multiple of granularity " + _granularity[i].tostring();
        }
        if (hasError)
        {
            std::cerr << err;
            throw NOMAD::Exception(__FILE__,__LINE__,err);
        }
    }
}


void NOMAD::GMesh::checkFrameSizeIntegrity(const NOMAD::Double &frameSizeExp,
                                          const NOMAD::Double &frameSizeMant) const
{
    // frameSizeExp must be an integer.
    // frameSizeMant must be 1, 2 or 5.

    bool hasError = false;
    std::string err = "Error: Integrity check";
    if (!frameSizeExp.isInteger())
    {
        hasError = true;
        err += " of frameSizeExp (" + frameSizeExp.tostring() + "): Should be integer.";
    }

    else if (frameSizeMant != 1.0 && frameSizeMant != 2.0 && frameSizeMant != 5.0)
    {
        hasError = true;
        err += " of frameSizeMant (" + frameSizeMant.tostring() + "): Should be 1, 2 or 5.";
    }

    if (hasError)
    {
        std::cerr << err;
        throw NOMAD::Exception(__FILE__,__LINE__,err);
    }
}


// Verify that setDeltas() gave the correct values.
void NOMAD::GMesh::checkSetDeltas(const size_t i,
                                  const NOMAD::Double &deltaMeshSize,
                                  const NOMAD::Double &deltaFrameSize) const
{
    bool hasError = false;
    std::string err = "Warning: setDeltas did not give good value";

    // Something might be wrong with setDeltas(), so double check.
    if (getdeltaMeshSize(i) != deltaMeshSize)
    {
        hasError = true;
        err += " for deltaMeshSize at index " + std::to_string(i);
        err += " Expected: " + deltaMeshSize.tostring();
        err += " computed: " + getdeltaMeshSize(i).tostring();
    }

    else if (getDeltaFrameSize(i) != deltaFrameSize)
    {
        hasError = true;
        err += " for deltaFrameSize at index " + std::to_string(i) + ".";
        err += " Expected: " + deltaFrameSize.tostring();
        err += " computed: " + getDeltaFrameSize(i).tostring();
    }

    if (hasError)
    {
        std::cerr << err << std::endl;
        throw NOMAD::Exception(__FILE__,__LINE__,err);
    }

}


/*-----------------------------------------------------------*/
/*              scale and project on the mesh                */
/*-----------------------------------------------------------*/
NOMAD::Double NOMAD::GMesh::scaleAndProjectOnMesh(size_t i, const NOMAD::Double &l) const
{
    const NOMAD::Double delta = getdeltaMeshSize(i);

    if (i < _n && _frameSizeMant.isDefined() && _frameSizeExp.isDefined() && delta.isDefined())
    {
        NOMAD::Double d = getRho(i) * l;
        return d.roundd() * delta;
    }
    else
    {
        std::ostringstream oss;
        oss << "GMesh: scaleAndProjectOnMesh cannot be performed.";
        oss << " i = ";
        oss << i;
        oss << " mantissa defined: " << _frameSizeMant.isDefined();
        oss << " exp defined: " << _frameSizeExp.isDefined();
        oss << " delta mesh size defined: " << delta.isDefined();
        throw NOMAD::Exception(__FILE__, __LINE__, oss.str());
    }
}


// Scale and project, using infinite norm.
NOMAD::ArrayOfDouble NOMAD::GMesh::scaleAndProjectOnMesh(
    const NOMAD::Direction &dir) const
{
    size_t n = _pbParams->getAttributeValue<size_t>("DIMENSION");
    NOMAD::ArrayOfDouble proj(n);

    NOMAD::Double infiniteNorm = dir.infiniteNorm();

    if (0 == infiniteNorm)
    {
        std::string err("GMesh: scaleAndProjectOnMesh: Cannot handle an infinite norm of zero");
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    for (size_t i = 0; i < n; ++i)
    {
        // Scaling and projection on the mesh
        proj[i] = this->scaleAndProjectOnMesh(i, dir[i] / infiniteNorm);
    }

    return proj;
}


NOMAD::Point NOMAD::GMesh::projectOnMesh(const NOMAD::Point& point,
                                         const NOMAD::Point& frameCenter) const
{
    // Projection on the mesh
    NOMAD::Point proj = point;
    auto delta = getdeltaMeshSize();
    // To avoid running around in circles
    const size_t maxNbTry = 10;

//    for (size_t i = 0; i < point.size(); ++i)
//    {
//
//        NOMAD::Double diffProjFrameCenter = proj[i] - frameCenter[i];
//        // Value which will be used in verifyPointIsOnMesh
//        proj[i] = diffProjFrameCenter / delta[i];
//        proj[i] = proj[i].round();
//        proj[i] = frameCenter[i] + proj[i] * delta[i];
//    }

    for (size_t i = 0; i < point.size(); ++i)
    {
        const NOMAD::Double deltaI = delta[i];
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
            // Some values are just ill-conditionned.
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


void NOMAD::GMesh::getLargerMantExp(NOMAD::Double &frameSizeMant, NOMAD::Double &frameSizeExp) const
{
    if (frameSizeMant == 1)
    {
        frameSizeMant = 2;
    }
    else if (frameSizeMant == 2)
    {
        frameSizeMant = 5;
    }
    else
    {
        frameSizeMant = 1;
        ++frameSizeExp;
    }
}
