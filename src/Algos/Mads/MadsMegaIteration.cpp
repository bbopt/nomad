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

#include <sstream>

#include "../../Algos/Mads/GMesh.hpp"
#include "../../Algos/Mads/MadsMegaIteration.hpp"
#include "../../Algos/Mads/MadsUpdate.hpp"
#include "../../Algos/Mads/MegaSearchPoll.hpp"

#include "../../Algos/EvcInterface.hpp"


void NOMAD::MadsMegaIteration::init()
{
    _name = getAlgoName() + NOMAD::MegaIteration::getName();
}

void NOMAD::MadsMegaIteration::startImp()
{
    // Create a NOMAD::MadsIteration for each frame center and each desired mesh size.
    // Use all xFeas and xInf available.
    // For now, not using other frame centers.
    size_t k = _k;  // Main iteration counter

    // Update main mesh and barrier.
    NOMAD::MadsUpdate update( this );
    update.start();
    update.run();
    update.end();

    // Now that update has used the previous MegaIteration success type, reset it
    setSuccessType(NOMAD::SuccessType::NOT_EVALUATED);

    // Verify mesh stop conditions.
    _mainMesh->checkMeshForStopping( _stopReasons );

    AddOutputDebug("Mesh Stop Reason: " + _stopReasons->getStopReasonAsString() );
    if ( ! _stopReasons->checkTerminate() )
    {
        // MegaIteration's barrier member is already in sub dimension.
        auto allXFeasPtr = _barrier->getAllXFeas();
        auto allXInfPtr  = _barrier->getAllXInf();

        // Create local copies of xFeas and xInf points
        std::vector<NOMAD::EvalPoint> allXFeas, allXInf;
        std::transform(allXFeasPtr.begin(), allXFeasPtr.end(), std::back_inserter(allXFeas),
                       [](NOMAD::EvalPointPtr evalPointPtr) -> NOMAD::EvalPoint { return *evalPointPtr; });
        std::transform(allXInfPtr.begin(), allXInfPtr.end(), std::back_inserter(allXInf),
                       [](NOMAD::EvalPointPtr evalPointPtr) -> NOMAD::EvalPoint { return *evalPointPtr; });

        // Compute the number of xFeas and xInf points we want to use, to get at
        // most MAX_ITERATION_PER_MEGAITERATION iterations.
        auto maxXFeas = allXFeas.size();
        auto maxXInf = allXInf.size();
        computeMaxXFeasXInf(maxXFeas, maxXInf);

        size_t nbPoints = 0;
        std::vector<NOMAD::EvalPoint>::const_iterator it;
        for (it = allXFeas.begin(), nbPoints = 0; it != allXFeas.end() && nbPoints < maxXFeas; ++it, nbPoints++)
        {
            std::shared_ptr<NOMAD::MadsIteration> madsIteration = std::make_shared<NOMAD::MadsIteration>(this , std::make_shared<NOMAD::EvalPoint>(*it), k, _mainMesh);
            _iterList.push_back(madsIteration);
            k++;
        }
        for (it = allXInf.begin(), nbPoints = 0; it != allXInf.end() && nbPoints < maxXInf; ++it, nbPoints++)
        {
            std::shared_ptr<NOMAD::MadsIteration> madsIteration = std::make_shared<NOMAD::MadsIteration>(this, std::make_shared<NOMAD::EvalPoint>(*it), k, _mainMesh);
            _iterList.push_back(madsIteration);
            k++;
        }

        // Add iterations for larger meshes
        // NOTE: A preliminary test gave not so good results
        // when using larger meshes. Skip for now.
        // TODO: Control if we have the bandwith to generate
        // a lot of Iterations, thus a lot of trial points,
        // or not.
        // TODO: Work on adding finer meshes too. There is more work
        // to do because a solution on a larger mesh is on the
        // main mesh, but not a solution on a finer mesh.
        /*
        if (xFeasDefined)
        {
            addIterationsForLargerMeshes(xFeas, k);
        }
        if (xInfDefined)
        {
            addIterationsForLargerMeshes(xInf, k);
        }
        */

        size_t nbIter = _iterList.size();

        AddOutputInfo(_name + " has " + NOMAD::itos(nbIter) + " iteration" + ((nbIter > 1)? "s" : "") + ".");

        AddOutputDebug("Iterations generated:");
        for (size_t i = 0; i < nbIter; i++)
        {
            // downcast from Iteration to MadsIteration
            std::shared_ptr<NOMAD::MadsIteration> madsIteration = std::dynamic_pointer_cast<NOMAD::MadsIteration>( _iterList[i] );

            if ( madsIteration == nullptr )
            {
                throw NOMAD::Exception(__FILE__, __LINE__, "Invalid shared pointer cast");
            }

            AddOutputDebug( _iterList[i]->getName());
            NOMAD::ArrayOfDouble meshSize  = madsIteration->getMesh()->getdeltaMeshSize();
            NOMAD::ArrayOfDouble frameSize = madsIteration->getMesh()->getDeltaFrameSize();
            auto frameCenter = madsIteration->getFrameCenter();
            AddOutputDebug("Frame center: " + frameCenter->display());
            auto previousFrameCenter = frameCenter->getPointFrom();
            AddOutputDebug("Previous frame center: " + (previousFrameCenter ? previousFrameCenter->display() : "NULL"));
            AddOutputDebug("Mesh size:  " + meshSize.display());
            AddOutputDebug("Frame size: " + frameSize.display());
        }
    }
}


// Commenting this out. Currently not used.
/*
bool NOMAD::MadsMegaIteration::addIterationsForLargerMeshes(const NOMAD::EvalPoint& x0, size_t &k)
{
    bool newMesh = false;

    // Compute new directions for x0, to enlarge mesh.
    size_t dim = _mainMesh->getSize();
    for (size_t i = 0; i < 2 * dim; i++)
    {
        // i = 0..dim-1 : going up
        // i = dim.. 2*dim-1 : going down
        size_t ii = (i < dim) ? i : i - dim;
        NOMAD::EvalPoint dirPoint(x0);
        NOMAD::Double delta = _mainMesh->getDeltaFrameSize(ii);
        dirPoint[ii] += delta;
        NOMAD::Direction dir = NOMAD::Point::vectorize(x0, dirPoint);

        // Create larger mesh.
        const auto largeMeshRef = std::shared_ptr<NOMAD::MeshBase>
            (new NOMAD::GMesh(*(dynamic_cast<NOMAD::GMesh*>(_mainMesh.get()))));
        auto largeMesh = std::shared_ptr<NOMAD::MeshBase>
            (new NOMAD::GMesh(*(dynamic_cast<NOMAD::GMesh*>(largeMeshRef.get()))));
        auto anisotropyFactor = _runParams->getAttributeValue<NOMAD::Double>("ANISOTROPY_FACTOR");
        bool anisotropicMesh = _runParams->getAttributeValue<bool>("ANISOTROPIC_MESH");

        if (largeMesh->enlargeDeltaFrameSize(dir, anisotropyFactor, anisotropicMesh))
        {
            // At least one new mesh was generated.
            newMesh = true;
            //std::string dirStr = "New direction " + dir.display();
            //AddOutputInfo(dirStr);
            std::shared_ptr<NOMAD::MadsIteration> madsIteration = std::make_shared<NOMAD::MadsIteration>(this, std::make_shared<NOMAD::EvalPoint>(dirPoint), k, largeMesh);
            _iterList.push_back(madsIteration);
            k++;

            // Reset largeMesh - Actually create a new one.
            largeMesh = std::shared_ptr<NOMAD::MeshBase>
                (new NOMAD::GMesh(*(dynamic_cast<NOMAD::GMesh*>(largeMeshRef.get()))));
        }

    }

    return newMesh;
}
*/


bool NOMAD::MadsMegaIteration::runImp()
{
    NOMAD::SuccessType bestSuccessYet = NOMAD::SuccessType::NOT_EVALUATED;

    std::string s;

    if ( _stopReasons->checkTerminate() )
    {
        s = "MegaIteration: stopReason = " + _stopReasons->getStopReasonAsString() ;
        AddOutputDebug(s);
        return false;
    }

    if (_iterList.empty())
    {
        throw NOMAD::Exception(__FILE__, __LINE__, "No iterations to run");
    }

    if (_runParams->getAttributeValue<bool>("GENERATE_ALL_POINTS_BEFORE_EVAL"))
    {
        MegaSearchPoll megaStep( this );
        megaStep.start();

        bool successful = megaStep.run();

        megaStep.end();

        if (successful)
        {
            bestSuccessYet = _megaIterationSuccess;
            s = "MadsMegaIteration: new success " + NOMAD::enumStr(bestSuccessYet);
            s += " stopReason = " + _stopReasons->getStopReasonAsString() ;
            AddOutputDebug(s);
        }

        // Note: Delta (frame size) will be updated in the Update step next time it is called.

        // End of running all the iterations
        // Update number of iterations - note: _k is atomic
        _k += _iterList.size();

    }
    else
    {
        for (size_t i = 0; i < _iterList.size(); i++)
        {
            // Get Mads ancestor to call terminate(k)
            NOMAD::Mads* mads = const_cast<NOMAD::Mads*>(dynamic_cast<const NOMAD::Mads*>(getParentOfType<NOMAD::Mads*>()));
            if (nullptr == mads)
            {
                throw NOMAD::Exception(__FILE__, __LINE__, "Mads MegaIteration without Mads ancestor");
            }
            if (_stopReasons->checkTerminate()
                || _stopReasons->testIf(NOMAD::EvalStopType::OPPORTUNISTIC_SUCCESS)
                || mads->terminate(_iterList[i]->getK()))
            {
                break;
            }

            // downcast from Iteration to MadsIteration
            std::shared_ptr<NOMAD::MadsIteration> madsIteration = std::dynamic_pointer_cast<NOMAD::MadsIteration>(_iterList[i]);

            if (madsIteration == nullptr)
            {
                throw NOMAD::Exception(__FILE__, __LINE__, "Invalid shared pointer cast");
            }

            madsIteration->start();

            bool iterSuccessful = madsIteration->run();          // Is this iteration successful
            // Compute MegaIteration success
            NOMAD::SuccessType iterSuccess = madsIteration->getSuccessType();
            if (iterSuccess > bestSuccessYet)
            {
                bestSuccessYet = iterSuccess;
            }

            madsIteration->end();

            if (iterSuccessful)
            {
                s = "MadsMegaIteration: new success " + NOMAD::enumStr(iterSuccess);
                AddOutputDebug(s);
            }

            // Update MegaIteration's stop reason
            if (_stopReasons->checkTerminate())
            {
                s = "MadsMegaIteration stop reason set to: " + _stopReasons->getStopReasonAsString();
                AddOutputDebug(s);
            }

            _nbIterRun++; // Count one more iteration.

            if (_userInterrupt)
            {
                hotRestartOnUserInterrupt();
            }

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
