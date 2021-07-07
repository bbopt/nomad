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

//#include "../../Algos/Mads/GMesh.hpp" // Code using GMesh is commented out
#include "../../Algos/Mads/Mads.hpp"
#include "../../Algos/Mads/MadsIteration.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Mads/MadsUpdate.hpp"
#include "../../Cache/CacheBase.hpp"
#include "../../Output/OutputQueue.hpp"


void NOMAD::MadsMegaIteration::init()
{
    setStepType(NOMAD::StepType::MEGA_ITERATION);
}


NOMAD::ArrayOfPoint NOMAD::MadsMegaIteration::suggest()
{
    // Create a single MadsIteration
    NOMAD::MadsIteration madsIteration (this, _k, _mainMesh);

    OUTPUT_DEBUG_START
    AddOutputDebug("Iteration generated:");
    AddOutputDebug(madsIteration.getName());
    NOMAD::ArrayOfDouble meshSize  = madsIteration.getMesh()->getdeltaMeshSize();
    NOMAD::ArrayOfDouble frameSize = madsIteration.getMesh()->getDeltaFrameSize();
    AddOutputDebug("Mesh size:  " + meshSize.display());
    AddOutputDebug("Frame size: " + frameSize.display());
    OUTPUT_DEBUG_END

    return madsIteration.suggest();
}


void NOMAD::MadsMegaIteration::observe(const std::vector<NOMAD::EvalPoint>& evalPointList)
{
    // Update cache with new points.
    NOMAD::EvalPoint evalPointFound;
    for (auto evalPoint : evalPointList)
    {
        if (NOMAD::CacheBase::getInstance()->find(evalPoint, evalPointFound))
        {
            // New eval for point already in cache
            NOMAD::CacheBase::getInstance()->update(evalPoint, NOMAD::EvalType::BB);
        }
        else
        {
            // Point is not in cache yet
            NOMAD::CacheBase::getInstance()->smartInsert(evalPoint);
        }
    }

    // Update barrier with new points.
    _barrier->updateRefBests();
    _barrier->updateWithPoints(evalPointList, NOMAD::EvalType::BB, NOMAD::ComputeType::STANDARD, _runParams->getAttributeValue<bool>("FRAME_CENTER_USE_CACHE"));

    // Update main mesh
    NOMAD::MadsUpdate update(this);
    update.start();
    update.run();
    update.end();

    OUTPUT_DEBUG_START
    AddOutputDebug("MegaIteration generated: " + getName());
    NOMAD::ArrayOfDouble meshSize  = _mainMesh->getdeltaMeshSize();
    NOMAD::ArrayOfDouble frameSize = _mainMesh->getDeltaFrameSize();
    AddOutputDebug("Mesh size:  " + meshSize.display());
    AddOutputDebug("Frame size: " + frameSize.display());
    OUTPUT_DEBUG_END
}


void NOMAD::MadsMegaIteration::startImp()
{
    // Update main mesh and barrier.
    NOMAD::MadsUpdate update( this );
    update.start();
    update.run();
    update.end();

    // Now that update has used the previous MegaIteration success type, reset it
    setSuccessType(NOMAD::SuccessType::NOT_EVALUATED);

    // Verify mesh stop conditions.
    _mainMesh->checkMeshForStopping( _stopReasons );

    OUTPUT_DEBUG_START
    AddOutputDebug("Mesh Stop Reason: " + _stopReasons->getStopReasonAsString());
    OUTPUT_DEBUG_END
}


bool NOMAD::MadsMegaIteration::runImp()
{
    NOMAD::SuccessType bestSuccessYet = NOMAD::SuccessType::NOT_EVALUATED;

    std::string s;

    if ( _stopReasons->checkTerminate() )
    {
        OUTPUT_DEBUG_START
        s = "MegaIteration: stopReason = " + _stopReasons->getStopReasonAsString() ;
        AddOutputDebug(s);
        OUTPUT_DEBUG_END
        return false;
    }

    // Get Mads ancestor to call terminate(k)
    NOMAD::Mads* mads = getParentOfType<NOMAD::Mads*>();
    if (nullptr == mads)
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "Mads MegaIteration without Mads ancestor");
    }
    if (!mads->terminate(_k))
    {
        // Create a single MadsIteration
        std::shared_ptr<NOMAD::MadsIteration> madsIteration = std::make_shared<NOMAD::MadsIteration>(this, _k, _mainMesh);

        OUTPUT_DEBUG_START
        AddOutputDebug("Iteration generated:");
        AddOutputDebug(madsIteration->getName());
        NOMAD::ArrayOfDouble meshSize  = madsIteration->getMesh()->getdeltaMeshSize();
        NOMAD::ArrayOfDouble frameSize = madsIteration->getMesh()->getDeltaFrameSize();
        AddOutputDebug("Mesh size:  " + meshSize.display());
        AddOutputDebug("Frame size: " + frameSize.display());
        OUTPUT_DEBUG_END

        madsIteration->start();

        bool iterSuccessful = madsIteration->run();
        // Compute MegaIteration success
        NOMAD::SuccessType iterSuccess = madsIteration->getSuccessType();
        if (iterSuccess > bestSuccessYet)
        {
            bestSuccessYet = iterSuccess;
        }

        madsIteration->end();

        if (iterSuccessful)
        {
            OUTPUT_DEBUG_START
            s = getName() + ": new success " + NOMAD::enumStr(iterSuccess);
            AddOutputDebug(s);
            OUTPUT_DEBUG_END
        }

        // Update MegaIteration's stop reason
        if (_stopReasons->checkTerminate())
        {
            OUTPUT_DEBUG_START
            s = getName() + " stop reason set to: " + _stopReasons->getStopReasonAsString();
            AddOutputDebug(s);
            OUTPUT_DEBUG_END
        }

        // Note: Delta (frame size) will be updated in the Update step next time it is called.

        if (_userInterrupt)
        {
            hotRestartOnUserInterrupt();
        }
    }

    // MegaIteration is a success if either a better xFeas or
    // a dominating or partial success for xInf was found.
    // See Algorithm 12.2 from DFBO.
    setSuccessType(bestSuccessYet);

    // return true if we have a partial or full success.
    return (bestSuccessYet >= NOMAD::SuccessType::PARTIAL_SUCCESS);
}


void NOMAD::MadsMegaIteration::display(std::ostream& os) const
{
    os << "MAIN_MESH " << std::endl;
    os << *_mainMesh ;
    NOMAD::MegaIteration::display(os);
}


void NOMAD::MadsMegaIteration::read(std::istream& is)
{
    // Set up structures to gather member info
    size_t k = 0;
    // Read line by line
    std::string name;
    while (is >> name && is.good() && !is.eof())
    {
        if ("MAIN_MESH" == name)
        {
            if (nullptr != _mainMesh)
            {
                is >> *_mainMesh;
            }
            else
            {
                std::string err = "Error: Reading a mesh onto a NULL pointer";
                std::cerr << err;
            }
        }
        else if ("ITERATION_COUNT" == name)
        {
            is >> k;
        }
        else if ("BARRIER" == name)
        {
            if (nullptr != _barrier)
            {
                is >> *_barrier;
            }
            else
            {
                std::string err = "Error: Reading a Barrier onto a NULL pointer";
                std::cerr << err;
            }
        }
        else
        {
            for (size_t i = 0; i < name.size(); i++)
            {
                is.unget();
            }
            break;
        }
    }

    setK(k);
}


std::ostream& NOMAD::operator<<(std::ostream& os, const NOMAD::MadsMegaIteration& megaIteration)
{
    megaIteration.display ( os );
    return os;
}


std::istream& NOMAD::operator>>(std::istream& is, NOMAD::MadsMegaIteration& megaIteration)
{

    megaIteration.read(is);
    return is;

}
