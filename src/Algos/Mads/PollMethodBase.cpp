
#include "../../Algos/AlgoStopReasons.hpp"
#include "../../Algos/EvcInterface.hpp"
#include "../../Algos/Mads/PollMethodBase.hpp"
#include "../../Output/OutputQueue.hpp"

void NOMAD::PollMethodBase::init()
{
    // A poll method must have a parent
    verifyParentNotNull();

}

void NOMAD::PollMethodBase::startImp()
{
    if ( ! _stopReasons->checkTerminate() )
    {
        // Create EvalPoints and snap to bounds and snap on mesh
        generateTrialPoints();

        // Stopping criterion
        if ( 0 == getTrialPointsCount() )
        {
            auto madsStopReasons = NOMAD::AlgoStopReasons<NOMAD::MadsStopType>::get ( _stopReasons );
            madsStopReasons->set( NOMAD::MadsStopType::MESH_PREC_REACHED );
        }

    }
}

bool NOMAD::PollMethodBase::runImp()
{

    bool foundBetter = false;
    if ( ! _stopReasons->checkTerminate() )
    {
        foundBetter = evalTrialPoints(this);
    }

    return foundBetter;
}

void NOMAD::PollMethodBase::endImp()
{
    // Compute hMax and update Barrier.
    postProcessing(NOMAD::EvcInterface::getEvaluatorControl()->getEvalType());
}


void NOMAD::PollMethodBase::generateTrialPoints()
{
    // Groups of variables.
    auto varGroups = _pbParams->getAttributeValue<NOMAD::ListOfVariableGroup>("VARIABLE_GROUP");;

    auto n = _pbParams->getAttributeValue<size_t>("DIMENSION");

    std::list<NOMAD::Direction> directionsSubSpace, directionsFullSpace;

    if (varGroups.size() == 0)
    {
        // Creation of the poll directions in the full space
        generateUnitPollDirections(directionsFullSpace,n);
    }
    else
    {
        for (auto vg : varGroups)
        {
            size_t nVG = vg.size();

            // Creation of the poll directions in the sub space of the variable group
            generateUnitPollDirections(directionsSubSpace,nVG);

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
                }
                vgIndex++;
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
    scaleAndProjectOnMesh(directionsFullSpace);

    OUTPUT_DEBUG_START
    for (auto dir : directionsFullSpace)
    {
        AddOutputDebug("Scaled and mesh projected poll direction: " + dir.display());
    }
    OUTPUT_DEBUG_END

    OUTPUT_INFO_START
    AddOutputInfo("Generate points for " + _name, true, false);
    OUTPUT_INFO_END

    // We need a frame center to start with.
    auto frameCenter = getIterationFrameCenter();
    if (nullptr == frameCenter || !frameCenter->ArrayOfDouble::isDefined() || frameCenter->size() != n)
    {
        std::string err("Invalid frame center: ");
        if (nullptr != frameCenter)
        {
            err += frameCenter->display();
        }
        throw NOMAD::Exception(__FILE__, __LINE__, err);
    }

    OUTPUT_DEBUG_START
    AddOutputDebug("Frame center: " + frameCenter->display());
    OUTPUT_DEBUG_END

    for (std::list<NOMAD::Direction>::iterator it = directionsFullSpace.begin(); it != directionsFullSpace.end() ; ++it)
    {
        NOMAD::Point pt(n);

        // pt = frame center + direction
        for (size_t i = 0 ; i < n ; ++i )
        {
            pt[i] = (*frameCenter)[i] + (*it)[i];
        }

        // Snap the points and the corresponding direction to the bounds
        if (snapPointToBoundsAndProjectOnMesh(pt, _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("LOWER_BOUND"),
                                              _pbParams->getAttributeValue<NOMAD::ArrayOfDouble>("UPPER_BOUND")))
        {
            if (pt != *frameCenter->getX())
            {
                // New EvalPoint to be evaluated.
                // Add it to the list.
                bool inserted = insertTrialPoint(NOMAD::EvalPoint(pt));

                OUTPUT_INFO_START
                std::string s = "Generated point";
                s += (inserted) ? ": " : " not inserted: ";
                s += pt.display();
                AddOutputInfo(s);
                OUTPUT_INFO_END
            }
        }

    }

    verifyPointsAreOnMesh(getName());
    updatePointsWithFrameCenter();

    OUTPUT_INFO_START
    AddOutputInfo("Generated " + NOMAD::itos(getTrialPointsCount()) + " points");
    AddOutputInfo("Generate points for " + _name, false, true);
    OUTPUT_INFO_END

}

void NOMAD::PollMethodBase::scaleAndProjectOnMesh(std::list<Direction> & dirs)
{
    // Scale the directions and project on the mesh

    std::shared_ptr<NOMAD::MeshBase> mesh = getIterationMesh();
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
